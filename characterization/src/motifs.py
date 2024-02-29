import re
import json
import subprocess
import numpy as np
import polars as pl
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from itertools import product
from Levenshtein import distance
from sklearn.manifold import MDS
from sklearn import metrics
from collections import Counter
from joblib import Parallel, delayed
from polars.exceptions import ComputeError


def define_IUPAC():
    """define IUPAC nucleotide code"""
    letters = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "R": ["A", "G"],
        "Y": ["C", "T"],
        "S": ["G", "C"],
        "W": ["A", "T"],
        "K": ["G", "T"],
        "M": ["A", "C"],
        "B": ["C", "G", "T"],
        "D": ["A", "G", "T"],
        "H": ["A", "C", "T"],
        "V": ["A", "C", "G"],
        "N": ["A", "C", "G", "T"],
    }

    ambiguous = {
        "R": ["A", "G"],
        "Y": ["C", "T"],
        "S": ["G", "C"],
        "W": ["A", "T"],
        "K": ["G", "T"],
        "M": ["A", "C"],
        "B": ["C", "G", "T"],
        "D": ["A", "G", "T"],
        "H": ["A", "C", "T"],
        "V": ["A", "C", "G"],
        "N": ["A", "C", "G", "T"],
    }

    encoding = {
        "A": [1, 0, 0, 0],
        "C": [0, 1, 0, 0],
        "G": [0, 0, 1, 0],
        "T": [0, 0, 0, 1],
        "R": [1, 0, 1, 0],
        "Y": [0, 1, 0, 1],
        "S": [0, 1, 1, 0],
        "W": [1, 0, 0, 1],
        "K": [0, 0, 1, 1],
        "M": [1, 1, 0, 0],
        "B": [0, 1, 1, 1],
        "D": [1, 0, 1, 1],
        "H": [1, 1, 0, 1],
        "V": [1, 1, 1, 0],
        "N": [1, 1, 1, 1],
    }

    encoding_extended = {
        "A": [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "C": [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "G": [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "T": [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "M": [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "R": [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "W": [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "S": [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "Y": [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "K": [0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "V": [1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "H": [1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "D": [1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "B": [0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "N": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    }

    IUPAC = {
        "letters": letters,
        "ambiguous": ambiguous,
        "encoding": encoding,
        "encoding_extended": encoding_extended,
    }

    with open(Path("benchmark/IUPAC.json"), "w") as f:
        json.dump(IUPAC, f)


def get_names(filename):
    seperator = "_peaks" if "_peaks" in filename else "_sel4.0"
    return filename.split(seperator)[0].split("_", 1)


class Colors:
    RED = "\033[31m"
    GREEN = "\033[32m"
    BLUE = "\033[34m"
    RESET = "\033[0m"


class Motifs:
    def __init__(self, result: str):
        self.file = result
        self.name = Path(self.file).stem
        self.motifs = (
            pl.scan_csv(Path(result).resolve()).select("motif").unique().collect()
        )
        self.motifs = self.motifs["motif"].to_list()

    def get_distance_matrix(self) -> np.ndarray:
        """get the levenshtein distance between two strings"""
        matrix = np.zeros((len(self.motifs), len(self.motifs)))
        for i in range(len(self.motifs)):
            for j in range(len(self.motifs)):
                matrix[i, j] = distance(self.motifs[i], self.motifs[j])

        return matrix

    def encode(self, IUPAC: str) -> np.ndarray:
        """encode motifs using IUPAC nucleotide code"""
        with open(Path(IUPAC).resolve(), "r") as f:
            IUPAC = json.load(f)

        encoding = IUPAC["extended"]
        encodeds = np.zeros(
            (
                len(self.motifs),
                max(len(motif) for motif in self.motifs),
                len(encoding["A"]),
            )
        )

        print("Encoding motifs...")

        for i, motif in enumerate(self.motifs):
            for j, n in enumerate(motif):
                encodeds[i, j] = encoding[n]
        return encodeds

    def record_to_fasta(self):
        """write motifs to fasta file"""
        out = self.file.with_suffix(".fasta")
        with open(out, "w") as f:
            for i, motif in enumerate(self.motifs):
                identifier = f">motif_{i + 1}"
                f.write(f"{identifier}\n{motif}\n")

    def refine_motifs(
        self,
        characterize: bool = False,
        repl: bool = False,
        kmer: bool = False,
        remove_short: bool = True,
        trim: bool = True,
    ):
        """refine motifs from raw outputs"""

        def replace(motif):
            match = re.search(r"(....).*?(...)$", motif)
            if match and len(match.group(0)) >= 9 and len(match.group(0)) <= 15:
                return (
                    match.group(1)
                    + "N"
                    * (len(match.group(0)) - len(match.group(1)) - len(match.group(2)))
                    + match.group(2)
                )
            else:
                return motif

        def is_base_safe(motif: str, IUPAC: dict, letter: list=["A", "C"]):
            
            if any(l in motif[0: 4] for l in letter):
                return motif
            
            for iupac, replacement in IUPAC.items():
                if iupac in motif:
                    if "A" in replacement:
                        return motif.replace(iupac, "A", 1)
                    if "C" in replacement and len(motif) <= 6:
                        return motif.replace(iupac, "C", 1)

            return None

        def get_most_frequent_kmer(kmer: str, freq: dict, iupac: dict):
            candidate = [kmer]
            for letter, options in iupac.items():
                if letter in kmer:
                    candidate.extend(kmer.replace(letter, char) for char in options)
            best = max(candidate, key=lambda x: freq.get(x, 0))
            if best == kmer:
                best = candidate[1]
            return best

        motifs = self.motifs[:]
        
        if trim:
            motifs = [m.lstrip("N").rstrip("N") for m in motifs]

        if remove_short:
            motifs = [m for m in motifs if len(m) >= 4]

        outdir = self.file.parent.as_posix().replace("results", "annotation")

        if repl:
            motifs = [replace(m) for m in motifs]

        with open(Path("benchmark/IUPAC.json").resolve(), "r") as f:
            IUPAC = json.load(f)["ambiguous"]

        if kmer:
            # TODO: use best filter approach for cisFinder
            kmer_freqs = Counter(motif[:3] for motif in motifs), Counter(
                motif[-3:] for motif in motifs
            )

            for i, motif in enumerate(motifs):

                if len(motif) > 2:
                    first_kmer = motif[:3]
                    if any(letter in first_kmer for letter in IUPAC.keys()):
                        motif = (
                            get_most_frequent_kmer(first_kmer, kmer_freqs[0], IUPAC)
                            + motif[3:]
                        )

                if len(motif) > 5:
                    last_kmer = motif[-3:]
                    if any(letter in last_kmer for letter in IUPAC.keys()):
                        motif = motif[:-3] + get_most_frequent_kmer(
                            last_kmer, kmer_freqs[1], IUPAC
                        )

                    motifs[i] = motif

        motifs = sorted(list(dict.fromkeys(motifs)), key=len, reverse=True)

        if characterize:

            outdir = outdir.replace("annotation", "characterize")
            motifs = [m for m in motifs if is_base_safe(m, IUPAC)]

        self.motifs = motifs[:]
        
        if not Path(outdir).exists():
            Path(outdir).mkdir(parents=True, exist_ok=True)

        with open(Path(f"{outdir}/{self.name}.json").resolve(), "w") as f:
            json.dump(self.motifs, f)

    def clustering(self):
        """cluster motifs using mmseqs2"""
        outdir = f"benchmark/cluster/{self.name}/"

        if not Path(outdir).exists():
            Path(outdir).mkdir(parents=True, exist_ok=True)

        self.record_to_fasta()

        # cluster motifs using mmseqs2
        cmd = [
            "mmseqs",
            "easy-cluster",
            f"{outdir}motifs.fasta",
            f"{outdir}clusterRes",
            f"{outdir}tmp",
            "-c",
            "0.8",
            "--target-search-mode",
            "1",
            "--min-seq-id",
            "0.3",
            "--alph-size",
            "15",
        ]

        try:
            with open(Path(outdir + "cluster.log").resolve(), "w") as f:
                subprocess.run(cmd, stdout=f, check=True)
        except Exception as e:
            print(f"Error occurred: {e} when running {cmd}")
            return None

    def global_clustering(
        self, rebase: str = "Methylase recognition sequences (plain text)..."
    ):
        """globally cluster current motifs with all motifs in ReBase"""

        outdir = f"benchmark/global_cluster/{self.name}/"

        if not Path(outdir).exists():
            Path(outdir).mkdir(parents=True, exist_ok=True)
            Path(f"{outdir}tmp/").mkdir(parents=True, exist_ok=True)

        self.record_to_fasta()

        record = self.hashing(rebase)

        # convert motifs and record into a single fasta file
        fasta = ""
        for motif in self.motifs:
            fasta += f">{motif}({self.name})\n{motif}\n"

        for motif in record.keys():
            fasta += f">{motif}(ReBase)\n{motif}\n"

        with open(Path(outdir + "cluster.fasta").resolve(), "w") as f:
            f.write(fasta)

        cluster = [
            "mmseqs",
            "easy-cluster",
            f"{outdir}cluster.fasta",
            f"{outdir}result",
            f"{outdir}tmp",
            "-c",
            "0.5",
            "--target-search-mode",
            "1",
            "--similarity-type",
            "2",
            "--alignment-mode",
            "1",
            "--min-seq-id",
            "0.3",
            "--alph-size",
            "15",
            "--adjust-kmer-len",
            "1",
        ]

        try:
            with open(Path(outdir + "cluster.log").resolve(), "w") as f:
                subprocess.run(cluster, stdout=f, check=True)
        except Exception as e:
            print(f"Error occurred: {e}")
            return None

        df = pl.scan_csv(
            f"{outdir}result_cluster.tsv",
            separator="\t",
            has_header=False,
            new_columns=["centroid", "member"],
        ).filter(
            (pl.col("member") != pl.col("centroid"))
            & ~(
                (pl.col("member").str.ends_with("(ReBase)"))
                & (pl.col("centroid").str.ends_with("(ReBase)"))
            )
        )
        df.collect().write_csv(f"{outdir}cluster.csv", has_header=True, separator=",")

    def search(self, rebase: str = "Methylase recognition sequences (plain text)..."):
        """perform a many to many search using mmseqs2 to search given motifs in ReBase
        :param rebase: str, rebase annotation file name
        """

        record = self.hashing(rebase)

        # write motifs to fasta file
        if not Path(f"benchmark/rebase/{rebase}.fasta").exists():
            with open(Path(f"benchmark/rebase/{rebase}.fasta").resolve(), "w") as f:
                for i, motif in enumerate(list(record.keys())):
                    f.write(f">rec_seq_{i + 1}\n{motif}\n")

        # create db and index
        if not Path(f"benchmark/rebase/mmseqs/").exists():
            Path(f"benchmark/rebase/mmseqs/").mkdir(parents=True, exist_ok=True)
            Path(f"benchmark/rebase/mmseqs/tmp/").mkdir(parents=True, exist_ok=True)

        create_db = [
            "mmseqs",
            "createdb",
            f"benchmark/rebase/{rebase}.fasta",
            f"benchmark/rebase/mmseqs/{rebase}",
            "--dbtype",
            "2",
            "-v",
            "3",
        ]

        create_index = [
            "mmseqs",
            "createindex",
            f"benchmark/rebase/mmseqs/{rebase}",
            "benchmark/rebase/mmseqs/tmp",
            "-k",
            "3",
            "--search-type",
            "3",
            "--alph-size",
            "15",
        ]

        subprocess.run(create_db, check=True)
        subprocess.run(create_index, check=True)

        # search motifs in ReBase

        query = f"benchmark/search/{self.name}"

        if not Path(query).exists():
            Path(query).mkdir(parents=True, exist_ok=True)
            Path(f"{query}/tmp/").mkdir(parents=True, exist_ok=True)

        if not self.file.with_suffix(".fasta").exists():
            self.record_to_fasta()

        create_db = [
            "mmseqs",
            "createdb",
            f"{self.file.with_suffix('.fasta')}",
            f"{query}/{self.name}",
            "--dbtype",
            "2",
            "-v",
            "3",
        ]

        search = [
            "mmseqs",
            "search",
            f"{query}/{self.name}",
            f"benchmark/rebase/mmseqs/{rebase}",
            f"{query}/result",
            f"{query}/tmp",
            "-e",
            "1",
            "--target-search-mode",
            "1",
            "--alph-size",
            "15",
            "--min-seq-id",
            "0.2",
            "--search-type",
            "0",
            "--alignment-mode",
            "3",
            "--gap-open",
            "1",
        ]

        convert = [
            "mmseqs",
            "convertalis",
            f"{query}/{self.name}",
            f"benchmark/rebase/mmseqs/{rebase}",
            f"{query}/result",
            f"{query}/result.m8",
            "--format-output",
            "query,target",
        ]

        subprocess.run(create_db, check=True)
        subprocess.run(search, check=True)
        subprocess.run(convert, check=True)

    def plot_motifs_distance(self, dissimilarity: str = "precomputed", out: str = None):
        """plot the distance matrix
        :param dissimilarity: str, precomputed (default) or euclidean
        """
        mds = MDS(
            n_components=2,
            dissimilarity=dissimilarity,
            random_state=1,
            normalized_stress=False,
        )

        if dissimilarity == "precomputed":
            results = mds.fit_transform(self.get_distance_matrix())
        else:
            # flatten the encoded motifs
            encodes = self.encode("benchmark/IUPAC.json").reshape(len(self.motifs), -1)
            results = mds.fit_transform(encodes)

        plt.figure(figsize=(8, 6))
        plt.scatter(results[:, 0], results[:, 1], edgecolors="k", c="orange", s=10)
        for i, (x, y) in enumerate(zip(results[:, 0], results[:, 1])):
            plt.text(x, y, self.motifs[i], fontsize=6)
        plt.xlabel("First MDS dimension")
        plt.ylabel("Second MDS dimension")
        plt.title(f"MDS plot of motif distances with {dissimilarity} distance")

        if out:
            plt.savefig(Path(out).resolve(), format="png", dpi=300)

        plt.show()

    def build_trie(self):
        """build a trie for quick search"""
        tries = []
        for motif in self.motifs:
            trie = {}
            trie[motif] = motif
            for i in range(1, len(motif)):
                suffix = motif[i:]
                prefix = motif[:i]
                if len(prefix) > 1:
                    trie[prefix] = motif
                if len(suffix) > 1:
                    trie[suffix] = motif
            tries.append(trie)
        return tries

    def annotate(self, outdir: str = None):
        """find the longest common substring"""

        def is_match(s1: str, s2: str, letters: dict):
            # check if two motifs are matched
            candidates = ["".join(p) for p in product(letters[s1])]
            for c in candidates:
                if c == s2:
                    return True

        def longest_common_substring(str1: str, str2: str):
            m = len(str1)
            n = len(str2)
            dp = [[0] * (n + 1) for _ in range(m + 1)]
            max_len = 0
            end = 0
            for i in range(1, m + 1):
                for j in range(1, n + 1):
                    if str1[i - 1] == str2[j - 1]:
                        dp[i][j] = dp[i - 1][j - 1] + 1
                        if dp[i][j] > max_len:
                            max_len = dp[i][j]
                            end = i
            return str1[end - max_len : end], str1, str2

        print(f"processing {self.file}...")

        outdir = (
            self.file.parent.as_posix().replace("results", "annotation")
            if not outdir
            else outdir
        )

        self.refine_motifs()

        if not Path(outdir).exists():
            Path(outdir).mkdir(parents=True, exist_ok=True)

        with open("benchmark/IUPAC.json", "r") as f:
            letters = json.load(f)["letters"]

        record = self.hashing()

        with open(Path(f"{outdir}/{self.name}.tsv").resolve(), "w") as f:
            f.write(
                "motif\tmatched\trebase\toverlap\tmeth_base\tmeth_type\tcomp_meth_base\tcomp_meth_type\tgenuine\tpacbio_only\n"
            )
            for m1 in self.motifs:
                results = Parallel(n_jobs=12)(
                    delayed(longest_common_substring)(m1, m2) for m2 in record.keys()
                )
                results = sorted(results, key=lambda x: len(x[0]), reverse=True)
                maximum = len(results[0][0])
                results = [r for r in results if len(r[0]) == maximum]
                for r in results:
                    annotations = "\t".join(v for v in record[r[2]].values())
                    f.write(
                        f"{r[1]}\t{r[0]}\t{r[2]}\t{(len(r[0]) / (len(r[2]))):.2f}\t{annotations}\n"
                    )
                    
    def refine(self):
        """refine motif from raw outputs using nucleotide code substitution"""
        
        outdir = self.file.parent.as_posix().replace("results", "characterize")
        dataset = get_names(self.name)[1]
        
        with open(Path(f"{outdir}/{self.name}.json").resolve(), "r") as f:
            motifs = json.load(f)

            print(f"{Colors.GREEN}refine motifs from {self.name}...{Colors.RESET}")

            nb_threads = "12"
            base_name = self.name
            path_diff_data = Path(
                f"data/samples/{dataset}/diff/{dataset}_difference_FN.RDS"
            ).as_posix()
            path_output = Path(f"{outdir.replace("characterize", "refinement")}/{self.name}").as_posix()
            list_motif = ",".join(motifs)
            genome = Path(f"data/samples/{dataset}/ref_genome/{dataset}.fasta")

            cmd = [
                "Rscript",
                Path("benchmark/src/nanodisco/code/refine.R").as_posix(),
                "-p",
                nb_threads,
                "-b",
                base_name,
                "-d",
                Path(path_diff_data).as_posix(),
                "-o",
                Path(path_output).as_posix(),
                "-m",
                list_motif,
                "-M",
                "all",
                "-r",
                genome,
            ]

            try:
                subprocess.run(cmd, check=True)
            except Exception as e:
                print(f"Error occurred: {e} when running {cmd}")
                return None

    def characterize(self, classifier: str):

        outdir = self.file.parent.as_posix().replace("results", "characterize")
        tool, dataset = get_names(self.name)

        if tool == "CoMos":
            self.refine_motifs()
        else:
            self.refine_motifs(characterize=True, repl=True)

        with open(Path(f"{outdir}/{self.name}.json").resolve(), "r") as f:
            motifs = json.load(f)

            print(f"{Colors.GREEN}characterizing {self.name}...{Colors.RESET}")

            nb_threads = "12"
            base_name = self.name
            path_diff_data = Path(
                f"data/samples/{dataset}/diff/{dataset}_difference_FN.RDS"
            ).as_posix()
            path_output = Path(f"{outdir}/{self.name}").as_posix()
            list_motif = ",".join(motifs)
            genome = Path(f"data/samples/{dataset}/ref_genome/{dataset}.fasta")

            if not Path(
                f"{path_output}/Motifs_classification_{base_name}_{classifier}_model.tsv"
            ).exists():

                cmd = [
                    "Rscript",
                    Path("benchmark/src/nanodisco/code/characterize.R").as_posix(),
                    "-p",
                    nb_threads,
                    "-b",
                    base_name,
                    "-d",
                    Path(path_diff_data).as_posix(),
                    "-o",
                    Path(path_output).as_posix(),
                    "-m",
                    list_motif,
                    "-t",
                    classifier,
                    "-r",
                    genome,
                ]

                try:
                    subprocess.run(cmd, check=True)
                except Exception as e:
                    print(f"Error occurred: {e} when running {cmd}")
                    return None
                
    def get_characterized_motifs(self, classifier: str = "nn"):
        outdir = Path(f"{self.file.parent.as_posix().replace("results", "characterize")}/{self.name}")
        
        try: 
            df = pl.scan_csv(Path(outdir / f"Motifs_classification_{self.name}_{classifier}_model.tsv"), separator="\t")
            return df
        except ComputeError:
            pass
        
    def write_best_classification_results(self, classification_results, base_name, path_output):
        classification_results = pd.read_csv(classification_results, sep='\t')
        
        classification_results['is_within'] = np.where((classification_results['mod_pos'] >= 0) & (classification_results['mod_pos'] <= classification_results['motif'].str.len() - 1), True, False)
        classification_results['replacing'] = classification_results['mod_type'].apply(lambda x: re.split('', x)[2])
        classification_results = classification_results[(classification_results['mod_pos'] >= 0) & (classification_results['mod_pos'] <= 3)]
        classification_results['safe_mod_pos'] = np.where(classification_results['is_within'], classification_results['mod_pos'] + 2, 1)
        classification_results['safe_base'] = np.where(classification_results['motif'].str.len() > 0, classification_results.apply(lambda row: re.split('', row['motif'])[row['safe_mod_pos']], axis=1), 'Z')
        classification_results['is_consistent'] = np.where(classification_results['is_within'] & (classification_results['replacing'] == classification_results['safe_base']), True, False)
        classification_results = classification_results[classification_results['is_consistent']]
        classification_results = classification_results.groupby('motif').apply(lambda x: x.loc[x['score'].idxmax()]).reset_index(drop=True)
        classification_results['clean_motif'] = classification_results.apply(lambda row: row['motif'][:row['mod_pos']] + row['mod_type'] + row['motif'][row['mod_pos'] + 2:], axis=1)
        classification_results['clear_mod_pos'] = classification_results['mod_pos'] + 1
        classification_results['clear_score'] = classification_results['score'].round(2)
        best_prediction = classification_results[['clean_motif', 'motif', 'clear_mod_pos', 'mod_type', 'clear_score']]
        best_prediction.columns = ['Characterized_motif', 'Motif', 'Predicted_position', 'Predicted_type', 'Prediction_score']

        file_name = f"Motifs_classification_{base_name}.tsv"
        output_file_name = path_output + file_name

        best_prediction.to_csv(output_file_name, sep='\t', index=False)

    @staticmethod
    def hashing(annotation: str = "Methylase recognition sequences (plain text)..."):
        """hashing all annotations for a quick search"""

        with open(f"benchmark/rebase/{annotation}.txt", "r") as f:
            lines = f.readlines()

        record = {}

        for line in lines:
            pattern = [
                "<rec_seq>",
                "<meth_base>",
                "<meth_type>",
                "<comp_meth_base>",
                "<comp_meth_type>",
                "<genuine>",
                "<pacbio_only>",
            ]
            for p in pattern:
                if line.startswith(p):
                    v = line.replace(p, "").strip("\n")
                    if p == "<rec_seq>":
                        seq = v
                        record[seq] = {}
                    else:
                        record[seq].update({p.strip("<>"): v})

        return record

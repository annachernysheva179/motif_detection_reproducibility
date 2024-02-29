import yaml
import seaborn as sns
import pandas as pd
import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt
from postprocess import get_pacbio_annotation
from Bio import SeqIO


def plot_motif_counts():
    """plot for every dataset the count of motifs for every tool before and after postprocessing
    :param motifs: motifs.yaml
    """

    motifs = yaml.load(
        open(Path("benchmark/motifs.yaml").resolve(), "r"), Loader=yaml.FullLoader
    )

    df = pd.DataFrame(columns=["program", "dataset", "motif", "suffix", "source"])

    with open("missing_characterize.txt", "a") as f:
        f.write("program\tdataset\tfilter\n")

    for program, v in motifs.items():
        for dataset, suffix in v.items():
            for s, detected in suffix.items():
                try:
                    postprocessed = pd.read_csv(
                        Path(
                            f"benchmark/postprocess/{program}/{dataset}/{s}/Motifs_classification_{program}_{dataset}_nn_model.tsv"
                        ).resolve(),
                        sep="\t",
                    )["Motif"].tolist()

                    assessment = pd.read_csv(
                        Path(
                            f"benchmark/postprocess/{program}/{dataset}/{s}/all_motifs_nn_model_classification.tsv"
                        ).resolve(),
                        sep="\t",
                    )
                except FileNotFoundError:
                    with open("missing_characterize.txt", "a") as f:
                        f.write(f"{program}\t{dataset}\t{s}\n")
                    continue

                counts = assessment.groupby("motif")["count"].sum()

                res = pd.DataFrame(
                    {
                        "program": [program] * len(detected),
                        "dataset": [dataset] * len(detected),
                        "motif": detected,
                        "suffix": [s] * len(detected),
                        "source": ["yaml"] * len(detected),
                    }
                )

                res["counts"] = res["motif"].map(counts)
                df = pd.concat([df, res], ignore_index=True)

                res_postprocessed = pd.DataFrame(
                    {
                        "program": [program] * len(postprocessed),
                        "dataset": [dataset] * len(postprocessed),
                        "motif": postprocessed,
                        "suffix": [s] * len(postprocessed),
                        "source": ["postprocess"] * len(postprocessed),
                    }
                )

                res_postprocessed["counts"] = res_postprocessed["motif"].map(counts)
                df = pd.concat([df, res_postprocessed], ignore_index=True)

    datasets = df["dataset"].unique()
    n_datasets = len(datasets)
    n_rows = n_datasets // 4 + (n_datasets % 4 > 0)
    n_cols = min(n_datasets, 4)
    fig, axs = plt.subplots(
        nrows=n_rows, ncols=n_cols, figsize=(10 * n_cols, 8 * n_rows), sharey=True
    )
    if isinstance(axs, np.ndarray):
        axs = axs.flatten()
    else:
        axs = [axs]

    for ax, dataset in zip(axs, datasets):
        df_dataset = df[df["dataset"] == dataset]
        sns.barplot(x="program", y="counts", hue="source", data=df_dataset, ax=ax)
        ax.set_title(dataset)

    plt.tight_layout()
    plt.savefig("benchmark/plot/motif_counts/motif_counts.png")


class Result:
    def __init__(self, tool: str, file: str) -> None:
        self.name = tool
        self.file = Path(file).resolve()

    def get_motifs(self):
        with open(self.file, "r") as f:
            raw = yaml.load(f, Loader=yaml.FullLoader)
            return raw[self.name]

    def evaluate(
        self,
        organism: str,
        genome_legnth: int,
        motifs: list,
        data_dir: str = "/data/samples/",
    ):
        """evaluate a result on a single dataset for a given tool"""
        organism = " ".join(organism.split(" ")[0:2])
        ground_truth = get_pacbio_annotation(organism)
        print(ground_truth)

    def benchmarking(self, data_dir: str = "data/samples/"):
        """benchmark all results on all datasets for a given tool"""
        data = self.get_motifs()
        all_motifs = set([motif for v in data.values() for motif in v])
        mapped = ""  # TODO
        for identifer, motif in data.items():
            dataset = mapped[identifer]
            genome = SeqIO.read(
                Path(data_dir, dataset, f"ref_genome/{dataset}.fasta").resolve(),
                "fasta",
            )
            length = len(genome)
            organism = (
                genome.description.replace(f"{genome.id}", "").strip().split(",")[0]
            )
            self.evaluate(organism, motif, length, data_dir)


if __name__ == "__main__":
    plot_motif_counts()

import re
import yaml
import polars as pl
from pathlib import Path
from postprocess import characterize, score, get_pacbio_annotation
from glob import glob
from Bio import SeqIO
from joblib import Parallel, delayed

# true posive: motifs macht rebase
# false positive: reported by tool not found in rebase + nanodisco
# false negative: not found in tools but reported in rebase
# true negative: not reported, no match in rebase + nanodisco (currently can not apply)


def map_name(datasets: list[str]):
    mapped = {}
    namespaces = [f.name for f in Path("/data/samples").iterdir() if f.is_dir()]
    for filename in datasets:
        combined = "".join(re.split(r"_|\.", filename))
        for name in namespaces:
            if name.replace("_", "").lower() in combined.lower():
                mapped[filename] = name
                break
    return mapped


def read_motifs(file: str):
    """read motifs from file"""
    df = pl.scan_csv(Path(file).resolve())
    return df.collect()["motif"].to_list()


def postprocess(func: str = "characterize"):
    """postprocess results using nanodisco
    :param func: Nanodisco function [characterize|score] to run
    """

    functions = {"characterize": characterize, "score": score}
    function = functions[func]

    def get_names(file: str) -> tuple[str, str]:
        names = file.split("/")[-1].rsplit(".", 1)[0].split("_")
        program = file.rsplit("/", 2)[1]
        dataset = "_".join(names[1:])
        return program, dataset

    motif, motifs, datasets = {}, {}, []

    # iterate over all result in benchmark/results
    for result in sorted(glob("/data/output/results/**/*.csv", recursive=True)):
        program, dataset_with_suffix = get_names(result)
        key, value = map_name([dataset_with_suffix]).popitem()
        dataset, suffix = value, key.replace(value, "")[1:]
        datasets.append(dataset)

        motif_with_suffix = {suffix: read_motifs(result)}

        (
            motif[dataset].update(motif_with_suffix)
            if dataset in motif
            else motif.update({dataset: motif_with_suffix})
        )

        if program in motifs:
            motifs[program].update(motif)
        else:
            motif.clear()
            motifs.update({program: {dataset: motif_with_suffix}})

    with open(Path("benchmark/motifs.yaml").resolve(), "w") as f:
        yaml.dump(motifs, f, width=float("inf"))

    inputs = []
    for program, v in motifs.items():
        for dataset, suffix in v.items():
            for s, m in suffix.items():
                if Path(
                    f"/data/output/benchmark/postprocess/{program}/{dataset}/{s}/all_motifs_nn_model_classification.tsv"
                ).exists():
                    continue
                curr_motifs = ",".join(m)
                if callable(function) and function.__name__ == "characterize":
                    inputs.append(
                        {
                            "nb_threads": "8",
                            "base_name": program + "_" + dataset,
                            "path_diff_data": Path(
                                f"/data/samples/{dataset}/diff/{dataset}_difference.RDS"
                            ).as_posix(),
                            "path_output": Path(
                                f"benchmark/postprocess/{program}/{dataset}/{s}"
                            ).as_posix(),
                            "list_motif": curr_motifs,
                            "type_model": "nn",
                            "genome": Path(
                                f"/data/samples/{dataset}/ref_genome/{dataset}.fasta"
                            ).as_posix(),
                        }
                    )
                elif callable(function) and function.__name__ == "score":
                    inputs.append(
                        {
                            "nb_threads": "8",
                            "path_diff_data": Path(
                                f"/data/samples/{dataset}/diff/{dataset}_difference.RDS"
                            ).as_posix(),
                            "path_reference": Path(
                                f"/data/samples/{dataset}/ref_genome/{dataset}.fasta"
                            ).as_posix(),
                            "base_name": program + "_" + dataset,
                            "path_output": Path(
                                f"benchmark/postprocess/{program}/{dataset}/{s}"
                            ).as_posix(),
                            "list_motif": motifs,
                        }
                    )
                else:
                    raise ValueError("Invalid function name")

    Parallel(n_jobs=1, verbose=10)(delayed(function)(**input) for input in inputs)


if __name__ == "__main__":
    postprocess()

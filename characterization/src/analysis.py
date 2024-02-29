import re
import json
import subprocess
import polars as pl
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from polars.exceptions import ComputeError
from glob import glob
from motifs import get_names


def nanodisco_motif(
    nb_threads: str,
    base_name: str,
    path_diff_data: str,
    path_output: str,
    genome: str,
    list_config: str = None,
    contigs_file: str = None,
    auto: bool = False,
):
    """adapted from motif.R"""

    cmd = [
        "Rscript",
        Path("benchmark/src/nanodisco/code/motif.R").as_posix(),
        "-p",
        nb_threads,
        "-b",
        base_name,
        "-d",
        Path(path_diff_data).as_posix(),
        "-o",
        Path(path_output).as_posix(),
        "-r",
        genome,
    ]

    cmd += ["--list_contig", list_config] if list_config is not None else []
    cmd += (
        ["--contigs_file", Path(contigs_file).as_posix()]
        if contigs_file is not None
        else []
    )
    cmd += ["-a"] if auto else []

    try:
        subprocess.run(cmd, check=True)
    except Exception as e:
        print(f"Error occurred: {e} when running {cmd}")
        return None


def find_default_motifs(result_dir: str):

    dfs = []

    for file in glob(
        f"{result_dir}/nanodisco/**/Motifs_classification_*_model.tsv",
        recursive=True,
    ):
        dfs.append(
            pl.scan_csv(file, has_header=True, separator="\t")
            .with_columns(
                pl.lit(file.split("_")[-2]).alias("classifier"),
                pl.lit(Path(file).parent.name).alias("dataset"),
                pl.lit("nanodisco").alias("name"),
            )
            .filter(pl.col("Prediction_score") > 0)
        )

    df = pl.concat(dfs, how="vertical").sort("dataset").unique().collect()
    df.write_csv(
        f"{result_dir}/nanodisco_classification.tsv",
        has_header=True,
        separator="\t",
    )

    df = (
        (
            (df.drop(["classifier", "name"]).unique())
            .rename({"Motif": "motif"})
            .rename({"Characterized_motif": "characterized_motif"})
            .rename({"Predicted_type": "meth_type"})
            .sort("motif")
        )
        .groupby("motif", "dataset", "characterized_motif", "meth_type")
        .agg(pl.col("Prediction_score").max().alias("score"))
    )

    df.write_csv(
        f"{result_dir}/nanodisco_motifs.tsv",
        has_header=True,
        separator="\t",
    )


def is_motif_characterized(result_dir: str, renew: bool = True, plot: bool = True):
    """Check if the motif is characterized"""

    classifer = result_dir.split("_")[-1]

    if renew:
        result = {
            "name": [],
            "dataset": [],
            "postprocess": [],
            "motif": [],
            "accepted": [],
            "occurred": [],
            "characterized": [],
            "characterized_motif": [],
            "meth_type": [],
            "meth_pos": [],
            "score": [],
        }

        for file in glob(f"{result_dir}/**/*.json", recursive=True):
            separator = "_peaks_|_sel4.0_"
            filename = Path(file).name
            tool, dataset = get_names(filename)
            postprocess = re.split(separator, filename)[1].replace(".json", "")

            with open(file, "r") as f:
                motifs = json.load(f)

            try:
                classification = pl.scan_csv(
                    Path(
                        f"{Path(file).with_suffix('')}/all_motifs_{classifer}_model_classification.tsv"
                    ),
                    has_header=True,
                    separator="\t",
                )
                classification = classification.filter(pl.col("count") > 0).collect()

            except ComputeError:
                classification = None

            try:
                characterization = pl.scan_csv(
                    Path(
                        f"{Path(file).with_suffix('')}/Motifs_classification_{Path(file).stem}_{classifer}_model.tsv"
                    ),
                    has_header=True,
                    separator="\t",
                )
                characterization = characterization.filter(
                    pl.col("Prediction_score") > 0
                ).collect()

            except ComputeError:
                characterization = None

            for motif in motifs:
                result["name"].append(tool)
                result["dataset"].append(dataset)
                result["postprocess"].append(postprocess)
                result["motif"].append(motif)
                signature_filename = f"Signature_center_detection_{Path(file).stem}_example_{motif}_v1.pdf"
                if Path(f"{Path(file).with_suffix('')}/{signature_filename}").exists():
                    result["accepted"].append(True)
                else:
                    result["accepted"].append(False)

                if classification is not None:

                    result["occurred"].append(
                        True if motif in classification["motif"].to_list() else False
                    )

                else:
                    result["occurred"].append(False)

                if characterization is not None:
                    try:
                        matched_types = characterization.filter(
                            pl.col("Motif") == motif
                        )["Predicted_type"].to_list()
                        matched_motifs = characterization.filter(
                            pl.col("Motif") == motif
                        )["Characterized_motif"].to_list()
                        matched_pos = characterization.filter(pl.col("Motif") == motif)[
                            "Predicted_position"
                        ].to_list()
                        matched_score = characterization.filter(
                            pl.col("Motif") == motif
                        )["Prediction_score"].to_list()

                        if (
                            not matched_types
                            or not matched_motifs
                            or not matched_pos
                            or not matched_score
                        ):
                            raise ValueError

                        result["meth_type"].append(matched_types)
                        result["characterized_motif"].append(matched_motifs)
                        result["meth_pos"].append(matched_pos)
                        result["score"].append(matched_score)

                        result["characterized"].append(True)

                    except ValueError:
                        result["meth_type"].append(None)
                        result["characterized"].append(False)
                        result["characterized_motif"].append(None)
                        result["meth_pos"].append(None)
                        result["score"].append(None)
                else:
                    result["characterized"].append(False)
                    result["meth_type"].append(None)
                    result["characterized_motif"].append(None)
                    result["meth_pos"].append(None)
                    result["score"].append(None)

        df = (
            pl.DataFrame(result)
            .explode("meth_type")
            .explode("characterized_motif")
            .explode("meth_pos")
            .explode("score")
            .unique(maintain_order=True)
        )

        df.write_csv(
            f"{result_dir}/motif_characterized.csv",
            has_header=True,
            separator=",",
        )

    # plot

    if plot:
        df = pd.read_csv(f"{result_dir}/motif_characterized.csv").drop(
            columns=[
                "characterized_motif",
                "meth_type",
                "meth_pos",
            ]
        )
        # Group by 'process' and 'tool', count 'characterized' and 'not characterized' motifs
        deduplicated = df.drop_duplicates(subset=["name", "dataset", "motif"])

        # Reshape the DataFrame to a long format
        df_tool = deduplicated.melt(
            id_vars=["name", "dataset", "motif"],
            value_vars=["accepted", "occurred", "characterized"],
            var_name="level",
            value_name="count",
        )

        # Separate the DataFrame into 'True' and 'False' DataFrames
        df_tool_true = df_tool[df_tool["count"] == True]
        df_tool_false = df_tool[df_tool["count"] == False]

        # Create a figure with two subplots
        fig, axs = plt.subplots(1, 2, figsize=(10, 5))

        # Create a bar plot for 'True' values
        plot = sns.countplot(
            data=df_tool_true,
            x="name",
            hue="level",
            ax=axs[0],
        )

        for p in plot.patches:
            height = p.get_height()
            if height > 0:
                plot.text(
                    p.get_x() + p.get_width() / 2.0,
                    height,
                    "%d" % int(p.get_height()),
                    fontsize=6,
                    color="black",
                    ha="center",
                    va="bottom",
                )

        plot.legend_.remove()

        axs[0].set_xlabel("")  # Remove x-axis label
        axs[0].set_ylabel("Number of accepted motifs")  # Set y-axis label

        # Create a bar plot for 'False' values
        plot = sns.countplot(
            data=df_tool_false,
            x="name",
            hue="level",
            ax=axs[1],
        )

        for p in plot.patches:
            height = p.get_height()
            if height > 0:
                plot.text(
                    p.get_x() + p.get_width() / 2.0,
                    height,
                    "%d" % int(p.get_height()),
                    fontsize=6,
                    color="black",
                    ha="center",
                    va="bottom",
                )

        axs[1].set_xlabel("")  # Remove x-axis label
        axs[1].set_ylabel("Number of rejected motifs")  # Set y-axis label

        fig.suptitle(
            "Number of motifs accepted, occurred and characterized by each tool"
        )

        plt.tight_layout()
        plt.savefig(f"{result_dir}/tagged_motifs_tools.png", dpi=300)

        fig, axs = plt.subplots(1, 3, figsize=(15, 7), sharey=True)
        for i, level in enumerate(["accepted", "occurred", "characterized"]):
            # Filter the DataFrame for the current level
            df_level = df_tool[df_tool["level"] == level]
            # Pivot the DataFrame to have datasets as columns and tools as index
            df_pivot = df_level.pivot_table(
                index="dataset",
                columns="name",
                values="count",
                aggfunc="sum",
            )

            # Create a stacked bar plot
            df_pivot.plot.bar(stacked=True, ax=axs[i], legend=False)
            axs[i].set_xlabel("")  # Remove x-axis label
            axs[i].set_title(level)  # Set the title

        # Create a legend for the figure
        handles, labels = axs[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper right")

        fig.suptitle(
            "Number of motifs accepted, occurred and characterized by each dataset"
        )

        # Adjust the spacing between the subplots
        fig.subplots_adjust(top=0.85)

        plt.tight_layout()
        plt.savefig(f"{result_dir}/tagged_motifs_datasets.png", dpi=300)
        plt.show()


def find_all_characterized_motifs(result_dir: str):
    """Find all characterized motifs"""

    df = (
        pl.scan_csv(
            f"{result_dir}/motif_characterized.csv",
            has_header=True,
            separator=",",
        )
        .drop_nulls(subset=["characterized_motif"])
        .select(["name", "dataset", "motif"])
        .unique()
        .collect()
    )

    df = df.group_by(["motif", "name"]).agg(pl.col("dataset"))
    df = df.pivot(index="motif", columns="name", values="dataset")

    df = df.with_columns(
        pl.col("meme").list.join(","),
        pl.col("memechip").list.join(","),
        pl.col("CoMoS").list.join(","),
        pl.col("vCNN").list.join(","),
        pl.col("cisFinder").list.join(","),
    )

    """
    meme = df.select(["dataset", "meme"]).explode("meme")
    memechip = df.select(["dataset", "memechip"]).explode("memechip")
    CoMoS = df.select(["dataset", "CoMoS"]).explode("CoMoS")
    cisFinder = df.select(["dataset", "cisFinder"]).explode("cisFinder")
    vCNN = df.select(["dataset", "vCNN"]).explode("vCNN")

    res = (
        vCNN.join(cisFinder, on=["dataset"], how="left")
        .join(meme, on=["dataset"], how="left")
        .join(memechip, on=["dataset"], how="left")
        .join(CoMoS, on=["dataset"], how="left")
    )

    df = (
        (
            df.explode("meme")
            .explode("memechip")
            .explode("CoMoS")
            .explode("cisFinder")
            .explode("vCNN")
        )
        .unique()
        .group_by(["meme", "memechip", "CoMoS"])
        .agg(pl.col("vCNN"), pl.col("cisFinder"))
    )
    
    res = df.with_columns(
        pl.when(df["CoMoS"].is_first_distinct()).then(df["CoMoS"]).otherwise(None),
        pl.when(df["meme"].is_first_distinct()).then(df["meme"]).otherwise(None),
        pl.when(df["memechip"].is_first_distinct())
        .then(df["memechip"])
        .otherwise(None),
        pl.when(df["vCNN"].is_first_distinct()).then(df["vCNN"]).otherwise(None),
    )
    """

    df.write_csv(
        f"{result_dir}/characterized_motifs_summary.tsv",
        has_header=True,
        separator="\t",
    )


def evaluate_plot(result_dir: str):

    def heatmap(data, **kwargs):
        pivot = data.pivot_table(
            index="dataset",
            columns="name",
            values="score",
        )
        sns.heatmap(pivot, **kwargs)

    Path.mkdir(Path(f"{result_dir}/plot"), exist_ok=True)

    df = pd.read_csv(f"{result_dir}/motif_characterized.csv")

    summary = (
        df.groupby(["name", "dataset", "postprocess"])
        .agg(
            is_characterized=("characterized", "mean"),
            not_occurred=("occurred", lambda x: 1 - x.mean()),
        )
        .reset_index()
    )

    reshaped = pd.melt(
        summary,
        id_vars=["name", "dataset", "postprocess"],
        value_vars=["is_characterized", "not_occurred"],
        var_name="metric",
        value_name="percentage",
    )

    # Create the FacetGrid, plot the boxplots for each metric
    g = sns.FacetGrid(
        reshaped,
        row="metric",
        height=4,
        aspect=1.5,
        margin_titles=True,
    )

    g.map_dataframe(
        sns.boxplot,
        x="name",
        y="percentage",
        hue="postprocess",
        palette="Set2",
    )

    g.set(xlabel=None)
    g.tick_params(axis="x", labelrotation=45, labelsize=8)

    legend_data = {
        "amb2.0": "amb2.0",
        "uBH_0.001": "naive",
        "uBH_0.001_peak_dist_2_min_cov_20_min_dist_20": "advanced",
        "uBH_0.001_peak_dist_2_min_cov_20_min_dist_20_k_3_kmer_quantile_0.25": "advanced_with_kmer",
    }

    handles, _ = g.axes[0][0].get_legend_handles_labels()
    labels = [legend_data[label] for label in _]
    g.figure.legend(
        handles=handles, labels=labels, loc="lower right", fontsize=6, ncol=4
    )

    g.tight_layout()

    g.savefig(f"{result_dir}/plot/rate_plot.png", dpi=300)

    # plot heatmap
    data = df[
        (df["postprocess"] == "uBH_0.001_peak_dist_2_min_cov_20_min_dist_20")
        | (df["postprocess"] == "amb2.0")
    ].drop(columns=["postprocess"])

    pivot = data.pivot_table(
        index="dataset",
        columns="name",
        values="score",
    )

    plt.figure(figsize=(8, 7))
    g = sns.heatmap(pivot, annot=True, cmap="YlGnBu", fmt=".2f")
    g.set(xlabel=None, ylabel=None)

    plt.savefig(f"{result_dir}/plot/heatmap.png", dpi=300)

    nanodisco = (
        pl.scan_csv(
            "benchmark/characterize_old/nanodisco_motifs.tsv",
            has_header=True,
            separator="\t",
        )
        .with_columns(pl.lit("nanodisco").alias("name"))
        .collect()
    )

    published = (
        pl.scan_csv(
            "benchmark/results/published_motif.csv", has_header=True, separator=","
        )
        .with_columns(
            pl.lit("published").alias("name"),
            pl.lit(100.00).alias("score"),
        )
        .collect()
    )

    df = (
        (
            pl.scan_csv(
                f"{result_dir}/motif_characterized.csv",
                has_header=True,
                separator=",",
            )
            .filter(pl.col("characterized") == True)
            .drop(["postprocess", "accepted", "occurred", "characterized", "meth_pos"])
        )
        .unique()
        .collect()
    )

    data = df.join(
        nanodisco,
        on=["name", "dataset", "motif", "score", "characterized_motif", "meth_type"],
        how="outer",
    )

    data = data.join(
        published,
        on=["name", "dataset", "motif", "score", "characterized_motif", "meth_type"],
        how="outer",
    )

    # change meth_type with multi values to unknown
    data = (
        data.group_by(["name", "dataset", "motif"]).agg(
            pl.col("meth_type"), pl.col("score").mean()
        )
    ).with_columns(
        pl.when(pl.col("meth_type").list.len() > 1)
        .then(pl.lit(None).alias("meth_type"))
        .otherwise(pl.col("meth_type").list.first().alias("meth_type"))
    )

    # add new column level
    data = (
        data.with_columns(pl.col("name").alias("level"))
        .group_by(["dataset", "motif"])
        .agg(pl.col("level"), pl.col("meth_type"), pl.col("score").mean())
    )

    data = data.with_columns(
        pl.when(pl.col("level").list.contains("published"))
        .then("published")
        .when(
            pl.col("level").list.contains("nanodisco"),
        )
        .then("nanodisco")
        .when(pl.col("level").list.len() > 1)
        .then("multiple tools")
        .otherwise(data["level"])
        .alias("detected by"),
    ).with_columns(pl.col("detected by").list.join(","))

    data = data.with_columns(
        pl.col("meth_type").list.drop_nulls().list.unique(),
    )

    data = data.with_columns(
        pl.when(pl.col("meth_type").list.len() != 0)
        .then(pl.col("meth_type").list.first())
        .otherwise(pl.lit("unknown").alias("meth_type")),
    )

    motifs = (
        data.filter(
            (pl.col("score").is_between(90, 100))
            & (pl.col("detected by") != "published")
        )
        .select("motif", "score")
        .to_pandas()
    )

    # plot scatter plot
    plt.figure(figsize=(10, 7))
    axes = sns.scatterplot(
        data=data.to_pandas(),
        x="motif",
        y="score",
        hue="dataset",
        style="detected by",
        size="meth_type",
        sizes={"unknown": 10, "6mA": 30, "4mC": 40, "5mC": 50},
        legend="full",
        palette="Set2",
    )

    for _, row in motifs.iterrows():
        text = row["motif"]
        t = re.sub(r"N+", lambda m: f"N({len(m.group(0))})", text)
        axes.annotate(t, (row["motif"], row["score"]), fontsize=6)

    axes.set(xlabel="Motif", ylabel="Score")
    axes.set_xticklabels([])
    axes.legend(loc="upper center", bbox_to_anchor=(0.5, -0.05), ncol=10, fontsize=6)

    plt.tight_layout()
    plt.savefig(f"{result_dir}/plot/scatter_plot.png", dpi=300)

    """

    g = sns.relplot(
        data=data.to_pandas(),
        x="motif",
        y="score",
        col="meth_type",
        hue="dataset",
        style="detected by",
        kind="scatter",
    )

    g._legend.remove()
    g.set(xlabel=None)
    g.set_xticklabels([])
    g.figure.legend(loc="upper center", ncol=12, fontsize=6)

    plt.tight_layout()
    g.savefig(f"{result_dir}/plot/scatter_plot_type.png", dpi=300)
    """

    # Display the plot
    plt.show()


if __name__ == "__main__":
    evaluate_plot("benchmark/characterize_knn")
    # is_motif_characterized("benchmark/characterize_knn", renew=True, plot=True)
    # find_all_characterized_motifs("benchmark/characterize_rf")

import yaml
import subprocess
import requests
import pandas as pd
import polars as pl
from io import StringIO
from pathlib import Path
from glob import glob
from bs4 import BeautifulSoup
from joblib import Parallel, delayed
from motifs import Motifs


HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3"
}


def score(
    nb_threads: str,
    path_diff_data: str,
    path_reference: str,
    base_name: str,
    path_output: str,
    list_motif: str,
    list_config: str = None,
    contigs_file: str = None,
):
    "adapted from score.R"

    cmd = [
        "Rscript",
        Path("benchmark/src/nanodisco/code/score.R").as_posix(),
        "-p",
        nb_threads,
        "-d",
        Path(path_diff_data).as_posix(),
        "-r",
        Path(path_reference).as_posix(),
        "-b",
        base_name,
        "-o",
        Path(path_output).as_posix(),
        "-m",
        list_motif,
    ]

    cmd += ["--list_contig", list_config] if list_config is not None else []
    cmd += (
        ["--contigs_file", Path(contigs_file).as_posix()]
        if contigs_file is not None
        else []
    )

    print(cmd)

    subprocess.run(cmd, check=True)

def nanodisco_characterize(
    nb_threads: str,
    base_name: str,
    path_diff_data: str,
    path_output: str,
    list_motif: str,
    type_model: str,
    genome: str,
    list_config: str = None,
    contigs_file: str = None,
):
    """adapted from characterize.R"""

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
        type_model,
        "-r",
        genome,
    ]

    cmd += ["--list_contig", list_config] if list_config is not None else []
    cmd += (
        ["--contigs_file", Path(contigs_file).as_posix()]
        if contigs_file is not None
        else []
    )

    print(cmd)

    try:
        subprocess.run(cmd, check=True)
    except Exception as e:
        print(f"Error occurred: {e} when running {cmd}")
        return None
    
    
def annotate(url: str = "http://rebase.neb.com/rebase/rebase.charts.html"):
    """"annotate motifs with ReBase annotations from 
    http://rebase.neb.com/rebase/rebase.charts.html
    :param motif: motif to annotate
    :param url: url to ReBase annotation lists
    """
    
    def get_all_annotations():
        
        if not Path("benchmark/rebase/rebase_list.html").exists():
            response = requests.get(url, headers=HEADERS)
            
            with open("benchmark/rebase/rebase_list.html", "w") as f:
                f.write(response.text)
                
            links = {}
            
            soup = BeautifulSoup(response.text, "html.parser")
            tables = soup.find_all('table', attrs={'bgcolor': '#FFFFFF', 'cellpadding': '6'})
            
            for link in tables[3].find_all('a'):
                links[link.text] = link.get('href')

            with open("benchmark/rebase/rebase_annotations.yaml", "w") as f:
                yaml.dump(links, f, width=float("inf"))
            
            
    def download_annotations(name: str):
        """download annotations from ReBase
        :param name: list name of the annotations to download
        """
        if not Path("benchmark/rebase/rebase_annotations.yaml").exists():           
            get_all_annotations()
        with open("benchmark/rebase/rebase_annotations.yaml", "r") as f:
            links = yaml.load(f, Loader=yaml.FullLoader)
            
        url = links[name]
        base = "http://rebase.neb.com/" if url.startswith("/cgi-bin") else "http://rebase.neb.com/rebase/"
        url = base + url
        response = requests.get(url, headers=HEADERS)
        if 'plain text' in name:
            with open(f"benchmark/rebase/{name}.txt", "wb") as f:
                f.write(response.content)
        else:
            with open(f"benchmark/rebase/{name}.html", "w") as f:
                f.write(response.text)
                
    

def get_pacbio_metadata(url: str = "http://rebase.neb.com/cgi-bin/pblist"):
    """get the list of all PacBio data on ReBase"""

    response = requests.get(url, headers=HEADERS)
    with open("benchmark/rebase/rebase.html", "w") as f:
        f.write(response.text)
    df = pd.read_html(
        Path("benchmark/rebase.html").resolve(),
        attrs={"bgcolor": "beige", "cellpadding": "4"},
    )[0]

    # remove useless columns
    df.columns = df.iloc[0]
    df = df.iloc[1:]
    df.set_index("Org#", inplace=True)

    df.to_csv("benchmark/rebase/rebase_pacbio_annot.csv")


def get_pacbio_annotation(org: str) -> pd.DataFrame:
    """search true positive motifs in PacBio data on ReBase"""

    df = pl.scan_csv(Path("benchmark/rebase/rebase_pacbio_annot.csv").resolve())
    selected = df.filter(pl.col("Organism").str.contains(org)).collect()
    value = selected.select("Org#").item()

    url = "http://rebase.neb.com/cgi-bin/pacbioget?{}+w".format(value)
    response = requests.get(url, headers=HEADERS)

    content = StringIO(response.text)

    df = pd.read_html(
        content,
        attrs={"cellpadding": "2", "bgcolor": "beige", "border": "0"},
    )

    table = df[0]
    table.dropna(axis=1, inplace=True)
    table.columns = table.iloc[1]
    motifs = table.iloc[1:][
        ["Motif", "Count", "Type", "Unique", "Genuine", "% Detected"]
    ]

    table = df[2]
    table.dropna(axis=1, inplace=True)
    table.columns = table.iloc[1]
    table = table.iloc[1:][
        ["Motif", "Count", "Type", "Unique", "Genuine", "% Detected"]
    ]

    motifs = pd.concat([motifs, table[1:]], axis=0)[1:]

    motifs.to_csv(f"benchmark/rebase/{org.replace(" ", "_")}_motifs.csv")

    return motifs

def process_characterization(output, classifier):
    if "nanodisco" not in output and "published" not in output and "cisFinder" in output:
        motifs = Motifs(Path(output))
        motifs.characterize(classifier=classifier)


if __name__ == "__main__":
    outputs = glob("benchmark/results/**/*.csv", recursive=True)
    Parallel(n_jobs=3)(delayed(process_characterization)(output, "rf") for output in outputs)
    

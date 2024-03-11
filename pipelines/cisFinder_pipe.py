import os
from datetime import datetime
import argparse
import numpy as np
import pandas as pd

start_time = datetime.now()

parser = argparse.ArgumentParser(description='Execute cisFinder')
parser.add_argument('--fasta', type=str, required=True,
                    help='An input fasta file in upper case')
parser.add_argument('--nmotifs', type=int, required=True,
                    help='Maximal number of motifs to produce by cisFinder')
parser.add_argument('--out', type=str, required=True,
                    help='Output directory for files')
parser.add_argument('--improve', type=bool, help='True/False if execution of patternTest needed; Default False ',
                     default=False)
parser.add_argument('--top', type=int, help='Number of top sequnces to be in the output; Default 20',
                     default=20)
#parser.add_argument('--p_val', type=str, required=True, help='p value used for filtering peaks, instead of . use -, for example for p_value 0.01 write 0_01')


#args = parser.parse_args()
args = vars(parser.parse_args())

#out_dir = "/data/output/cisF_output/"
input_path = args["fasta"]
out_dir = args["out"]
top_n = args["top"]
nmotifs = args["nmotifs"]
#p_val = args["p_val"]

base = os.path.basename(input_path).split(".fasta")[0]

output_path_find = out_dir +  "cisFinder_" + base + "_motifs.tsv"
output_path_cluster = out_dir + "cisFinder_" + base + "_motifs_clustered.tsv"
output_path_imporved = out_dir +  "cisFinder_"+ base + "_motifs_improved.tsv"


os.system(r"awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' "+input_path + "  > input.fasta")

input_fasta = "input.fasta"
print("Executing cisFinder function patternFind: ")
os.system("./patternFind -i {} -o {} -n {}".format(input_fasta,output_path_find, nmotifs))
print()
print("Executing cisFinder function patternClust: ")
os.system("./patternCluster -i {} -o {}".format(output_path_find,output_path_cluster))

final_file = ""

if args["improve"]:
    print()
    print("Executing cisFinder function patternTest: ")
    os.system("./patternTest -i {} -f {} -o {}".format(output_path_cluster,input_fasta,output_path_imporved))
    final_file = output_path_imporved
else:
    final_file = output_path_cluster

def reverse_complement(input_sequence):
    reverse_reg = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "V": "B", "H": "D",
                   "D": "H", "B": "V", "M": "K",
                   "R": "Y", "W": "W", "S": "S",
                   "Y": "R", "K": "M"}

    reverse_complement = ""

    for index in range(len(input_sequence)-1, -1, -1):
        base = input_sequence[index]
        complement = reverse_reg[base]
        reverse_complement += complement

    return reverse_complement


def get_file_data(path):

    file_data = {}
    table_data = []
    first_lines = ""
    with open(path, "r") as file:
        previous_line = ""
        for line in file:
            if not (line.startswith("Parameters") or line.startswith("Headers")):
                if line.startswith(">"):
                    motif_line = line
                elif line[0].isdigit() and previous_line != "Zvalue\n":
                    table_data.append(line.rstrip().split("\t"))
                elif previous_line.rstrip().split("\t") in table_data:
                    file_data[motif_line] = table_data
                    table_data = []
                else:
                    pass
                previous_line = line
            else:
                first_lines += line
    return file_data

# Define the PPM function
def ppm(C):
    return C / sum(C)

# Define the S function, B is a background
def s(C, B):
    return np.log2(ppm(C) / B)

# Define the ICtotal function
def ic_total(C):
    return np.log2(len(C))

# Define the U function
def u(C):
    return -np.sum(ppm(C) * np.log2(ppm(C)), axis=0)

# Define the ICfinal function
def ic_final(C):
    return ic_total(C) - u(C)

# Define the IC function
def ic(C):
    return ppm(C) * ic_final(C)


def calc_ic_df(pfm_df):
    B = 0.25
    ppm_df = pfm_df.apply(ppm, axis=1)
    s_df = ppm_df.apply(lambda x: s(x, B))
    pfm_df
    ic_df = pfm_df.apply(ic, axis=1)
    ic_df = ic_df.transpose()
    ic_t = ic_df.sum(axis=0)
    ic_t = ic_t.round(2)
    return ic_t.to_list()

def trim_motif (motif, ic_t, reverse):

    began = False
    new_motif = ""
    ran = range(0)
    if reverse:
        ran = reversed(range(len(motif)))
    else:
        ran = range(0, len(motif))
    new_ic = []
    for i in ran:
        if ic_t[i]>1.5 :
            began = True
            new_motif += motif[i]
            new_ic.append(ic_t[i])
            began = True
        elif ic_t[i]<=1.5 and began:
            new_motif += motif[i]
            new_ic.append(ic_t[i])

    if reverse:
        return new_motif[::-1], list(reversed(new_ic))

    return trim_motif(new_motif, new_ic, True)


def process_motifs(file_data,base,out_dir):
    new_data = {}
    output_path =out_dir +  "cisFinder_{}.csv".format(base)

    with open(output_path, "w") as file:
        file.write("motif,reverse,freq,fdr\n")
        for line, table in file_data.items():
            att = line.split("\t")

            motif = att[1].rstrip()
            pfm_df = pd.DataFrame(table, columns=["index", "A", "C", "G", "T"])
            pfm_df.set_index('index', inplace=True)
            pfm_df = pfm_df.astype(int)
            ic_t = calc_ic_df(pfm_df)
            new_motif, new_ic = trim_motif(motif, ic_t, False)
            reverse = reverse_complement(new_motif)
            freq = att[3]
            fdr = att[8]
            if len(new_motif) >= 3:
                file.write("{},{},{},{}\n".format(new_motif, reverse, freq, fdr))


file_data = get_file_data(final_file)
process_motifs(file_data,base,out_dir)
end_time = datetime.now()
print("Done postprocessig. The result is saved into file {} .".format(out_dir+"cisFinder_{}.csv".format(base)))
print('Duration: {}'.format(end_time - start_time))
os.system("rm input.fasta")
os.system("rm {}".format(output_path_find))
os.system("rm {}".format(output_path_cluster))
if args["improve"]:
    os.system("rm {}".format(output_path_imporved))



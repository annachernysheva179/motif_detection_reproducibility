import os
import sys
import json
import pandas as pd
import numpy as np

from datetime import datetime
import argparse
from collections import OrderedDict

start_time = datetime.now()

parser = argparse.ArgumentParser(description='Run MEME-ChIP')
parser.add_argument('--fasta', metavar='-f', type=str, required=True,
                    help='An input fasta file')
parser.add_argument('--nmotifs', metavar='-n', type=int, required=True,
                    help='Number of motifs')
parser.add_argument('--out_dir', metavar='-o', type=str, required=True,
                    help='Output directory for meme files')
parser.add_argument('--out_file', metavar='-of', type=str, required=True,
                    help='Output CSV file')


args = vars(parser.parse_args())

input_path = args["fasta"]
out_dir = args["out_dir"]
out_file = args["out_file"]
motifs = args["nmotifs"]
current_path = os.environ.get('PATH', '')
new_path = '/data/meme/bin'
new_path_value = '{}:{}'.format(current_path, new_path)
os.environ['PATH'] = new_path_value


command = "meme-chip -o {} -order 2 -streme-nmotifs {} -meme-nmotifs {} {}".format(out_dir, motifs, motifs, input_path)
os.system(command)

with open(out_dir + '/meme-chip.html', 'r') as fcc_file:
    lines = fcc_file.readlines()

ans = ""
i = 7
line = lines[i]
while "};" not in line:
    line = lines[i].lstrip().rstrip()
    ans += line
    i += 1
ans1 = ans[11:-1]


json_object = json.loads(ans1)

def ic_total(C):
    return np.log2(len(C))

def u(C):
    return -np.sum(C * np.log2(C), axis=0)

def ic_final(C):
    return ic_total(C) - u(C)

def ic_unique(C):
    return ic_final(C) * C

def calc_ic_df(pfm_df):
    ic_df = pfm_df.apply(ic_unique, axis=1)
    ic_df = ic_df.transpose()
    ic_t = ic_df.sum(axis=0)
    ic_t = ic_t.round(2)
    return ic_t.to_list()

def trim_motif(motif, ic_t, reverse):
    began = False
    new_motif = ""
    ran = range(0)
    if reverse:
        ran = reversed(range(len(motif)))
    else:
        ran = range(0, len(motif))
    new_ic = []
    for i in ran:
        if ic_t[i]>1:
            began = True
            new_motif += motif[i]
            new_ic.append(ic_t[i])
            began = True
        elif ic_t[i]<=1 and began:
            new_motif += motif[i]
            new_ic.append(ic_t[i])
    if reverse:
        return new_motif[::-1], list(reversed(new_ic))
    return trim_motif(new_motif, new_ic, True)


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


motifs = []
visited=[]
for dict in json_object['motifs']:
    motif = dict["consensus"]
    ind = motif.find("-")
    if ind != -1:
        motif = motif[ind+1:]
    if "centrimo_sites" in dict.keys():
        pwm = dict["pwm"]
        df_pwm = pd.DataFrame(pwm)
        ic_list = calc_ic_df(df_pwm)
        new_motif, new_ic = trim_motif(motif, ic_list, False)
        if new_motif not in visited and len(new_motif) >= 3:
            motifs.append([new_motif, dict["sig"]])
            reverse = reverse_complement(new_motif)
            visited.append(reverse)
            visited.append(new_motif)


df = pd.DataFrame(motifs, columns=["motif", "E-val"])
end_time = datetime.now()
print('Duration: {}'.format(end_time - start_time))
df.to_csv(out_file)



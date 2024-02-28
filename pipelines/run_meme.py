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


command = "meme -minw 4 -maxw 14 -o {} -nmotifs {} {}".format(out_dir, motifs, input_path)
os.system(command)

with open(out_dir + '/meme.html', 'r') as fcc_file:
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
    motif = dict["id"]
    ind = motif.find("-")
    if ind != -1:
        motif = motif[ind+1:]
    motifs.append([motif, dict["evalue"]])


df = pd.DataFrame(motifs, columns=["motif", "E-val"])
end_time = datetime.now()
print('Duration: {}'.format(end_time - start_time))
df.to_csv(out_file)



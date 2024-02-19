import os
import sys
import numpy as np
from Bio import SeqIO
import pandas as pd


def calculate_mutations(fasta_file):
    with open(fasta_file, "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))

    M = len(records)
    N = len(records[0].seq)

    count = np.zeros((N, 3))
    ave_seg = np.zeros(1)
    L_count = np.zeros((N + 1, 2))

    ref_seq = records[0].seq

    for j in range(1, M):
        seq = records[j].seq
        for k in range(N):
            count[k, 0] = k

            if seq[k] != ref_seq[k]:
                count[k, 1] += 1

    mu = np.mean(count[:, 1])
    count[:, 2] = count[:, 1] / mu
    ave_seg[0] = np.mean(count[:, 2])
    max_mut = np.max(count[:, 1])
    L_count[0:int(max_mut) + 1, 1] = np.linspace(0, max_mut, int(max_mut) + 1)

    for j in range(N):
        mutation = count[j, 1]
        L_count[int(mutation), 0] += 1

    L_count[:, 0] /= mu
    L_count[:, 1] /= N

    header = ["Position", "c/u", "MCount"]
    header2 = ["Prob", "c/u"]

    return count, ave_seg, L_count, header, header2


# Function to write results to CSV
def write_results(output_file1, output_file2, output_file3, count, ave_seg, L_count, header, header2):
    count_df = pd.DataFrame(count, columns=header)
    count_df.to_csv(output_file1, index=False)

    ave_seg_df = pd.DataFrame({"AveSeg": ave_seg})
    ave_seg_df.to_csv(output_file2, index=False)

    L_count_df = pd.DataFrame(L_count, columns=header2)
    L_count_df.to_csv(output_file3, index=False)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <fasta_file> <output_file1> <output_file2> <output_file3>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    output_file1 = sys.argv[2]
    output_file2 = sys.argv[3]
    output_file3 = sys.argv[4]

    count, ave_seg, L_count, header, header2 = calculate_mutations(fasta_file)
    write_results(output_file1, output_file2, output_file3, count, ave_seg, L_count, header, header2)

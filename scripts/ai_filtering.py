import sys
from Bio import SeqIO
from joblib import load
import numpy as np

def compute_features(seq):
    seq = str(seq)
    length = len(seq)
    gc = (seq.count('G') + seq.count('C')) / length if length > 0 else 0
    at_skew = abs(seq.count('A') - seq.count('T')) / length if length > 0 else 0
    entropy = -sum(seq.count(n)/length * np.log2(seq.count(n)/length) 
                   for n in "ACGT" if seq.count(n) > 0)
    max_homopolymer = max([len(s) for s in seq.replace("C"," ").replace("G"," ").replace("T"," ").split()])
    return [length, gc, at_skew, entropy, max_homopolymer]

def filter_reads(input_r1, input_r2, output_r1, output_r2, model_path):
    clf = load(model_path)
    with open(output_r1, "w") as out1, open(output_r2, "w") as out2:
        for rec1, rec2 in zip(SeqIO.parse(input_r1, "fastq"), SeqIO.parse(input_r2, "fastq")):
            X = compute_features(rec1.seq)
            prediction = clf.predict([X])[0]
            if prediction == 1:
                SeqIO.write(rec1, out1, "fastq")
                SeqIO.write(rec2, out2, "fastq")

if __name__ == "__main__":
    input_r1 = sys.argv[1]
    input_r2 = sys.argv[2]
    output_r1 = sys.argv[3]
    output_r2 = sys.argv[4]
    model_path = sys.argv[5]
    filter_reads(input_r1, input_r2, output_r1, output_r2, model_path)

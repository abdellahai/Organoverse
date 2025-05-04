import sys
from Bio import SeqIO
from Bio.Seq import Seq

def find_overlap(seq, min_overlap=100):
    """Find overlap between start and end of the sequence."""
    for i in range(min_overlap, len(seq) // 2):
        if seq[:i] == seq[-i:]:
            return i
    return None

def circularize(input_fasta, output_fasta):
    record = list(SeqIO.parse(input_fasta, "fasta"))[0]
    seq = str(record.seq)
    overlap = find_overlap(seq, min_overlap=100)
    if overlap:
        new_seq = seq[:-overlap]
        print(f"[INFO] Circularization performed. Overlap of {overlap} bp removed.")
    else:
        new_seq = seq
        print("[INFO] No overlap found. Sequence kept linear.")
    
    record.seq = Seq(new_seq)
    record.description += " [circularized]" if overlap else " [linear]"
    SeqIO.write(record, output_fasta, "fasta")

if __name__ == "__main__":
    circularize(sys.argv[1], sys.argv[2])

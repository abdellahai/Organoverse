import argparse
import os
import subprocess
import glob
from Bio import SeqIO
import sys

def merge_references(ref_dir, merged_fasta):
    """Concatenate all FASTA in ref_dir into one merged_fasta."""
    fasta_paths = []
    for ext in ("*.fasta", "*.fa", "*.fna"):
        fasta_paths.extend(glob.glob(os.path.join(ref_dir, ext)))
    if not fasta_paths:
        raise ValueError(f"No FASTA files found in reference directory: {ref_dir}")
    if len(fasta_paths) == 1:
        # Single reference, no merge needed
        return fasta_paths[0]
    # Merge multiple references
    with open(merged_fasta, 'w') as w:
        for path in fasta_paths:
            with open(path) as r:
                w.write(r.read())
    print(f"[INFO] Merged {len(fasta_paths)} references into {merged_fasta}")
    return merged_fasta

def run_fastp(r1, r2, out_r1, out_r2, html, json, threads=4,
              qualified_quality_phred=20, length_required=50):
    cmd = [
        "fastp",
        "-i", r1, "-I", r2,
        "-o", out_r1, "-O", out_r2,
        "--html", html, "--json", json,
        "--qualified_quality_phred", str(qualified_quality_phred),
        "--length_required", str(length_required),
        "--detect_adapter_for_pe",
        "--thread", str(threads)
    ]
    print("[INFO] Running fastp:", " ".join(cmd))
    subprocess.run(cmd, check=True)

def extract_with_bbduk(r1, r2, ref_fasta, out_dir, k=31, hdist=1, threads=4):
    os.makedirs(out_dir, exist_ok=True)
    out1 = os.path.join(out_dir, "extracted_R1.fq.gz")
    out2 = os.path.join(out_dir, "extracted_R2.fq.gz")
    cmd = [
        "bbduk.sh",
        f"in1={r1}", f"in2={r2}",
        f"out1={out1}", f"out2={out2}",
        f"ref={ref_fasta}",
        f"k={k}", f"hdist={hdist}",
        f"threads={threads}", "minlen=50"
    ]
    print("[INFO] Running BBduk:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return out1, out2

def remove_chimeras(fasta_in, fasta_out, ref_fasta, threads=4):
    cmd = [
        "vsearch",
        "--uchime_ref", fasta_in,
        "--db", ref_fasta,
        "--nonchimeras", fasta_out,
        "--threads", str(threads)
    ]
    print("[INFO] Removing chimeras with vsearch:", " ".join(cmd))
    subprocess.run(cmd, check=True)

def remove_duplicates(fasta_in, fasta_out):
    cmd = [
        "seqkit", "rmdup",
        "-s", fasta_in,
        "-o", fasta_out
    ]
    print("[INFO] Removing duplicates with seqkit:", " ".join(cmd))
    subprocess.run(cmd, check=True)

def remove_contamination(fastq1, fastq2, out1, out2, kraken_db, threads=4):
    cmd = [
        "kraken2",
        "--db", kraken_db,
        "--paired",
        "--output", "/dev/null",
        "--unclassified-out", f"{out1},{out2}",
        "--threads", str(threads),
        fastq1, fastq2
    ]
    print("[INFO] Removing contamination with kraken2:", " ".join(cmd))
    subprocess.run(cmd, check=True)

def subsample_reads(fastq1, fastq2, out1, out2, fraction, seed=42):
    cmd1 = f"seqtk sample -s{seed} {fastq1} {fraction} > {out1}"
    cmd2 = f"seqtk sample -s{seed} {fastq2} {fraction} > {out2}"
    print("[INFO] Subsampling reads:", cmd1)
    subprocess.run(cmd1, shell=True, check=True)
    print("[INFO] Subsampling reads:", cmd2)
    subprocess.run(cmd2, shell=True, check=True)

def ai_filter(fastq1, fastq2, out1, out2, model_path):
    cmd = [
        "python", "scripts/ai_filtering.py",
        fastq1, fastq2, out1, out2, model_path
    ]
    print("[INFO] Running AI filtering:", " ".join(cmd))
    subprocess.run(cmd, check=True)

def main():
    parser = argparse.ArgumentParser(
        description="Preprocess WGS reads for organelle assembly"
    )
    parser.add_argument("--r1", required=True, help="Input WGS R1 FASTQ")
    parser.add_argument("--r2", required=True, help="Input WGS R2 FASTQ")
    parser.add_argument("--refdir", required=True,
                        help="Directory containing one or more reference FASTA(s)")
    parser.add_argument("--output", required=True,
                        help="Output directory")
    parser.add_argument("--kraken_db", default=None,
                        help="Path to Kraken2 database for contamination removal")
    parser.add_argument("--model", default=None,
                        help="Path to AI model (.joblib) for final filtering")
    parser.add_argument("--subsample", type=float, default=1.0,
                        help="Fraction to subsample reads (0 < fraction ≤ 1)")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    # Prepare directories
    os.makedirs(args.output, exist_ok=True)
    work = os.path.join(args.output, "work")
    os.makedirs(work, exist_ok=True)

    # Merge or select reference FASTA(s)
    merged_ref = merge_references(args.refdir, os.path.join(work, "merged_refs.fasta"))

    # Step 1: Extract reads matching ref organelles
    ext_dir = os.path.join(work, "extracted")
    ext_r1, ext_r2 = extract_with_bbduk(args.r1, args.r2, merged_ref, ext_dir,
                                        threads=args.threads)

    # Step 2: Quality trimming and filtering
    qc_r1 = os.path.join(work, "qc_R1.fq.gz")
    qc_r2 = os.path.join(work, "qc_R2.fq.gz")
    qc_html = os.path.join(work, "fastp.html")
    qc_json = os.path.join(work, "fastp.json")
    run_fastp(ext_r1, ext_r2, qc_r1, qc_r2, qc_html, qc_json,
              threads=args.threads)

    # Step 3: Chimera removal
    chimera_free = os.path.join(work, "nochim.fasta")
    remove_chimeras(qc_r1, chimera_free, merged_ref, threads=args.threads)

    # Step 4: Duplicate removal
    dedup_fasta = os.path.join(work, "dedup.fasta")
    remove_duplicates(chimera_free, dedup_fasta)

    # Step 5: Contamination removal
    if args.kraken_db:
        clean_r1 = os.path.join(work, "clean_R1.fq")
        clean_r2 = os.path.join(work, "clean_R2.fq")
        remove_contamination(chimera_free, qc_r2, clean_r1, clean_r2,
                             args.kraken_db, threads=args.threads)
    else:
        clean_r1, clean_r2 = dedup_fasta, qc_r2

    # Step 6: Subsample reads (optional)
    if args.subsample < 1.0:
        sub_r1 = os.path.join(work, "sub_R1.fq")
        sub_r2 = os.path.join(work, "sub_R2.fq")
        subsample_reads(clean_r1, clean_r2, sub_r1, sub_r2,
                        args.subsample)
        proc_r1, proc_r2 = sub_r1, sub_r2
    else:
        proc_r1, proc_r2 = clean_r1, clean_r2

    # Step 7: AI-based final filter
    if args.model:
        final_r1 = os.path.join(args.output, "final_R1.fq")
        final_r2 = os.path.join(args.output, "final_R2.fq")
        ai_filter(proc_r1, proc_r2, final_r1, final_r2, args.model)
    else:
        final_r1, final_r2 = proc_r1, proc_r2

    print("[DONE] Preprocessing complete.")
    print("Final reads:", final_r1, final_r2)

if __name__ == "__main__":
    main()


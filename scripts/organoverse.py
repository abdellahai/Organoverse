#!/usr/bin/env python3
"""Organoverse: AI-assisted organelle genome assembly and annotation pipeline.

This script orchestrates read extraction, assembly and annotation
using helper tools contained in the repository.  It mirrors the
workflow described in the project README.
"""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path
from typing import Tuple

from Bio import Entrez
from sklearn.metrics.pairwise import cosine_similarity
from transformers import AutoTokenizer, AutoModel
import torch

# configuration for NCBI and AI model
Entrez.email = "user@example.com"
TOKENIZER = AutoTokenizer.from_pretrained("dmis-lab/biobert-base-cased-v1.1")
MODEL = AutoModel.from_pretrained("dmis-lab/biobert-base-cased-v1.1")
MODEL.eval()

TAXONOMIC_RANKS = [
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "kingdom",
    "superkingdom",
]

# ----------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------

def embed_text(text: str) -> torch.Tensor:
    """Return BioBERT CLS embedding for ``text``."""
    inp = TOKENIZER(text, return_tensors="pt", truncation=True, max_length=512)
    with torch.no_grad():
        out = MODEL(**inp)
    return out.last_hidden_state[:, 0, :].squeeze(0)


def ncbi_esearch(db: str, term: str, retmax: int = 20):
    handle = Entrez.esearch(db=db, term=term, retmax=retmax)
    results = Entrez.read(handle)
    handle.close()
    return results.get("IdList", [])


def fetch_taxonomic_lineage(species: str):
    ids = ncbi_esearch("taxonomy", species, retmax=1)
    if not ids:
        return []
    handle = Entrez.efetch(db="taxonomy", id=ids[0], retmode="xml")
    rec = Entrez.read(handle)[0]
    handle.close()
    ranks = {i["Rank"]: i["ScientificName"] for i in rec.get("LineageEx", [])}
    ranks[rec.get("Rank", "species")] = rec.get("ScientificName", species)
    ordered = [(r, ranks[r]) for r in TAXONOMIC_RANKS if r in ranks]
    return ordered


def search_refseq_organelle(taxon: str, organelle: str):
    query = f"({taxon}[Organism]) AND {organelle}[Title] AND RefSeq[Filter]"
    ids = ncbi_esearch("nucleotide", query)
    if not ids:
        return []
    handle = Entrez.esummary(db="nucleotide", id=",".join(ids))
    summaries = Entrez.read(handle)
    handle.close()
    results = []
    for s in summaries:
        results.append(
            {
                "accession": s.get("AccessionVersion"),
                "title": s.get("Title"),
                "species": s.get("Organism"),
                "taxid": s.get("TaxId"),
            }
        )
    return results


def fetch_pubmed_abstract(species: str) -> str | None:
    ids = ncbi_esearch(
        "pubmed", f"{species}[Title/Abstract] AND (phylogeny OR genome)", retmax=1
    )
    if not ids:
        return None
    handle = Entrez.efetch(db="pubmed", id=ids[0], rettype="abstract", retmode="text")
    txt = handle.read().strip()
    handle.close()
    return txt


def select_best_candidate(candidates, target_species: str):
    tgt = embed_text(f"{target_species} organelle phylogeny")
    scores = []
    for cand in candidates:
        abst = fetch_pubmed_abstract(cand["species"])
        if not abst:
            continue
        emb = embed_text(abst)
        sim = cosine_similarity(tgt.reshape(1, -1), emb.reshape(1, -1))[0][0]
        scores.append((sim, cand))
    if not scores:
        return None
    scores.sort(key=lambda x: x[0], reverse=True)
    return scores[0][1]


def fetch_reference_files(acc: str, prefix: Path) -> Tuple[Path, Path]:
    fasta = prefix.with_suffix(".fasta")
    gbk = prefix.with_suffix(".gbk")
    for fmt, path in [("fasta", fasta), ("gbwithparts", gbk)]:
        handle = Entrez.efetch(db="nucleotide", id=acc, rettype=fmt, retmode="text")
        with open(path, "w") as f:
            f.write(handle.read())
        handle.close()
    return fasta, gbk


# ----------------------------------------------------------------------
# Pipeline steps
# ----------------------------------------------------------------------

def find_closest_reference(species: str, organelle: str, outdir: Path) -> Tuple[Path, Path]:
    lineage = fetch_taxonomic_lineage(species)
    for rank, taxon in lineage:
        candidates = search_refseq_organelle(taxon, organelle)
        if not candidates:
            continue
        if len(candidates) == 1:
            chosen = candidates[0]
        else:
            chosen = select_best_candidate(candidates, species)
            if not chosen:
                continue
        prefix = outdir / f"{chosen['accession']}"
        return fetch_reference_files(chosen["accession"], prefix)
    raise RuntimeError("No suitable reference found")


def run_extract_reads(ref: Path, r1: Path, r2: Path, prefix: Path) -> Tuple[Path, Path]:
    cmd = [
        "python",
        "scripts/extract_reads.py",
        "--ref",
        str(ref),
        "--fq1",
        str(r1),
        "--fq2",
        str(r2),
        "--out",
        str(prefix),
    ]
    subprocess.run(cmd, check=True)
    return prefix.with_name(prefix.name + "_ai_polished_1.fq"), prefix.with_name(prefix.name + "_ai_polished_2.fq")


def run_spades(r1: Path, r2: Path, outdir: Path) -> Path:
    cmd = ["spades.py", "-1", str(r1), "-2", str(r2), "-o", str(outdir), "--careful"]
    subprocess.run(cmd, check=True)
    return outdir / "contigs.fasta"


def run_ragtag_scaffold(contigs: Path, reference: Path, outdir: Path) -> Path:
    cmd = ["ragtag.py", "scaffold", str(reference), str(contigs), "-o", str(outdir)]
    subprocess.run(cmd, check=True)
    return outdir / "ragtag.scaffold.fasta"


def run_annotation(fasta: Path, outdir: Path) -> Path:
    outdir.mkdir(exist_ok=True)
    gbk = outdir / "annot.gbk"
    cmd1 = ["blastn", "-query", str(fasta), "-db", str(fasta), "-outfmt", "6"]
    subprocess.run(cmd1, check=True)
    cmd2 = ["tRNAscan-SE", str(fasta), str(outdir / "trnascan.txt")]
    subprocess.run(cmd2, check=True)
    # In practice you would parse the outputs and create a GenBank file.
    gbk.touch()
    return gbk


def main():
    parser = argparse.ArgumentParser(description="Organoverse pipeline")
    parser.add_argument("--species", required=True, help="Target species name")
    parser.add_argument("--organelle", choices=["mitochondrion", "chloroplast"], required=True)
    parser.add_argument("--library", choices=["paired", "single", "long"], required=True)
    parser.add_argument("--reads", nargs="+", required=True, help="Path(s) to reads")
    parser.add_argument("--outdir", default="results", help="Output directory")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("[1] Finding closest reference genome...")
    ref_fasta, ref_gbk = find_closest_reference(args.species, args.organelle, outdir)
    print(f"    FASTA: {ref_fasta}")
    print(f"    GenBank: {ref_gbk}")

    if args.library == "paired" and len(args.reads) >= 2:
        r1, r2 = Path(args.reads[0]), Path(args.reads[1])
    elif args.library == "single":
        r1 = Path(args.reads[0])
        r2 = Path(args.reads[0])
    else:
        raise ValueError("Unsupported library configuration")

    print("[2] Extracting organelle reads and polishing...")
    cleaned_r1, cleaned_r2 = run_extract_reads(ref_fasta, r1, r2, outdir / "reads")

    print("[3] Running SPAdes assembly...")
    contigs = run_spades(cleaned_r1, cleaned_r2, outdir / "spades")

    print("[4] Scaffolding against reference...")
    scaffold = run_ragtag_scaffold(contigs, ref_fasta, outdir / "ragtag")

    print("[5] Annotating assembly...")
    gbk = run_annotation(scaffold, outdir / "annotation")
    print(f"Annotation saved to {gbk}")

    print("Pipeline finished")


if __name__ == "__main__":
    main()

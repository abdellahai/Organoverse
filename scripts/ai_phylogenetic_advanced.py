import argparse
import os
import time
from Bio import Entrez, SeqIO
from transformers import AutoTokenizer, AutoModel
import torch
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

# Set your email to comply with NCBI policies
Entrez.email = "idrissi.azami.abdellah@gmail.com"

# Load SciBERT model
print("[INFO] Loading SciBERT model...")
tokenizer = AutoTokenizer.from_pretrained("allenai/scibert_scivocab_uncased")
model = AutoModel.from_pretrained("allenai/scibert_scivocab_uncased")

def embed_text(text):
    """Return the SciBERT embedding for a given text (taxon summary)."""
    inputs = tokenizer(text, return_tensors="pt", truncation=True, max_length=128)
    with torch.no_grad():
        outputs = model(**inputs)
    return outputs.last_hidden_state[:, 0, :].numpy()[0]  # CLS token

def fetch_taxon_summary(taxid):
    """Fetch the scientific name and brief lineage summary for a taxid."""
    handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
    rec = Entrez.read(handle)[0]
    name = rec["ScientificName"]
    lineage = rec.get("Lineage", "")
    summary = f"{name}. Lineage: {lineage}"
    return summary

def search_refseq(name, organelle):
    """Search for RefSeq organelle genomes for a given taxon name."""
    query = f"{name}[Organism] AND {organelle}[Filter] AND srcdb_refseq[PROP]"
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=10)
    return Entrez.read(handle)["IdList"]

def get_taxid(genus):
    """Get NCBI Taxonomy ID for a given genus."""
    handle = Entrez.esearch(db="taxonomy", term=genus)
    rec = Entrez.read(handle)
    return rec["IdList"][0] if rec["IdList"] else None

def get_higher_taxa(taxid):
    """Get lineage (excluding root) for the taxon."""
    handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
    rec = Entrez.read(handle)[0]
    return [t for t in rec["LineageEx"] if t["Rank"] in ("family", "order", "class", "phylum", "kingdom")]

def ai_phylogenetic_search(genus, organelle):
    """AI-based selection of closest higher taxon using SciBERT embeddings."""
    taxid = get_taxid(genus)
    if not taxid:
        print("[ERROR] TaxID not found for genus.")
        return []

    target_summary = fetch_taxon_summary(taxid)
    target_emb = embed_text(target_summary)

    lineage = get_higher_taxa(taxid)
    candidates = []

    for node in lineage:
        name = node["ScientificName"]
        ids = search_refseq(name, organelle)
        if ids:
            desc = fetch_taxon_summary(node["TaxId"])
            emb = embed_text(desc)
            sim = cosine_similarity([target_emb], [emb])[0][0]
            candidates.append((name, sim, ids))

    if not candidates:
        print("[INFO] No RefSeq genomes found in higher taxa.")
        return []

    candidates.sort(key=lambda x: x[1], reverse=True)
    print(f"[AI] Closest match: {candidates[0][0]} (similarity={candidates[0][1]:.3f})")
    return candidates[0][2]

def download_genomes(ids, outdir):
    os.makedirs(outdir, exist_ok=True)
    for seq_id in ids:
        handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        path = os.path.join(outdir, f"{record.id}.fasta")
        SeqIO.write(record, path, "fasta")
        print(f"[+] Downloaded {path}")
        time.sleep(0.5)

def main():
    parser = argparse.ArgumentParser(description="AI-enhanced RefSeq organelle genome retriever")
    parser.add_argument("--genus", required=True, help="Target genus name (e.g. Camellia)")
    parser.add_argument("--organelle", choices=["chloroplast", "mitochondrion"], required=True)
    parser.add_argument("--output", default="refseq_genomes", help="Output directory")

    args = parser.parse_args()

    ids = search_refseq(args.genus, args.organelle)
    if ids:
        print(f"[INFO] Found RefSeq organelle genomes for genus {args.genus}")
        download_genomes(ids, args.output)
    else:
        print(f"[INFO] No genus-level match found. Initiating AI-based search...")
        ids = ai_phylogenetic_search(args.genus, args.organelle)
        if ids:
            download_genomes(ids, args.output)
        else:
            print("[WARNING] No suitable genomes found even with AI assistance.")

if __name__ == "__main__":
    main()
    print("[INFO] Script completed.")
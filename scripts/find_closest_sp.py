# Organoverse AI-powered Reference Genome Selector
# Version: v2.3 with BioBERT + taxonomy fallback + CLI + FASTA + CSV + accurate species display

import os
import time
import argparse
import csv
from Bio import Entrez
from sklearn.metrics.pairwise import cosine_similarity
from transformers import AutoTokenizer, AutoModel
import torch

# --- CONFIGURATION ---
Entrez.email = "doctor.abdellah.idrissi.azami@um6ss.ma"
NCBI_API_DELAY = 0.4
MAX_RETRY_ATTEMPTS = 3
MIN_SIMILARITY_THRESHOLD = 0.4
TAXONOMIC_RANKS = ["species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"]

print("[AI] Loading BioBERT model (dmis-lab/biobert-base-cased-v1.1)...")
tokenizer = AutoTokenizer.from_pretrained("dmis-lab/biobert-base-cased-v1.1")
model = AutoModel.from_pretrained("dmis-lab/biobert-base-cased-v1.1")
model.eval()

def ncbi_request_with_retry(func, *args, **kwargs):
    for attempt in range(MAX_RETRY_ATTEMPTS):
        try:
            time.sleep(NCBI_API_DELAY)
            return func(*args, **kwargs)
        except Exception as e:
            print(f"[NCBI Retry {attempt + 1}] Error: {e}")
            time.sleep(1.5 * (attempt + 1))
    return None

def embed_text(text):
    try:
        inputs = tokenizer(text, return_tensors="pt", truncation=True, max_length=512)
        with torch.no_grad():
            outputs = model(**inputs)
        return outputs.last_hidden_state[:, 0, :].squeeze().numpy()
    except Exception as e:
        print(f"  [AI Error] Embedding failed: {e}")
        return None

def fetch_taxonomic_lineage(species_name):
    handle = ncbi_request_with_retry(Entrez.esearch, db="taxonomy", term=species_name, retmode="xml")
    if not handle:
        return []
    search = Entrez.read(handle)
    handle.close()
    if not search['IdList']:
        return []
    taxid = search['IdList'][0]
    handle = ncbi_request_with_retry(Entrez.efetch, db="taxonomy", id=taxid, retmode="xml")
    lineage_data = Entrez.read(handle)
    handle.close()

    ranks = {item['Rank']: item['ScientificName'] for item in lineage_data[0].get('LineageEx', []) if item['Rank'] in TAXONOMIC_RANKS}
    current_rank = lineage_data[0].get('Rank', 'species')
    current_name = lineage_data[0].get('ScientificName', species_name)
    if current_rank in TAXONOMIC_RANKS:
        ranks[current_rank] = current_name
    ordered = [(rank, ranks[rank]) for rank in TAXONOMIC_RANKS if rank in ranks]
    return ordered

def fetch_pubmed_abstract(species_name):
    try:
        query = f"{species_name}[Title/Abstract] AND (phylogeny OR genome OR organelle)"
        handle = ncbi_request_with_retry(Entrez.esearch, db="pubmed", term=query, retmax=1)
        if not handle:
            return None
        result = Entrez.read(handle)
        handle.close()
        if not result['IdList']:
            return None
        pubmed_id = result['IdList'][0]
        handle = ncbi_request_with_retry(Entrez.efetch, db="pubmed", id=pubmed_id, rettype="abstract", retmode="text")
        if not handle:
            return None
        abstract = handle.read()
        handle.close()
        return abstract.strip()
    except Exception as e:
        print(f"  [PubMed Error] Failed fetching abstract for {species_name}: {e}")
        return None

def fetch_fasta(accession):
    try:
        handle = ncbi_request_with_retry(Entrez.efetch, db="nucleotide", id=accession, rettype="fasta", retmode="text")
        if not handle:
            return None
        fasta = handle.read()
        handle.close()
        return fasta
    except Exception as e:
        print(f"  [FASTA Error] Could not fetch FASTA for {accession}: {e}")
        return None

def save_summary_to_csv(species, organelle, result, filename="organoverse_selection_summary.csv"):
    try:
        with open(filename, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["species", "accession", "title", "taxid", "organelle"])
            writer.writeheader()
            writer.writerow({
                "species": result["species"],
                "accession": result["accession"],
                "title": result["title"],
                "taxid": result["taxid"],
                "organelle": organelle
            })
        print(f"[Output] Summary saved to: {filename}")
    except Exception as e:
        print(f"  [CSV Error] Could not save summary: {e}")

def search_refseq_organelle(taxon_name, organelle_type):
    query = f"({taxon_name}[Organism]) AND RefSeq[Filter] AND {organelle_type}[Title] AND complete genome[Title]"
    print(f"[NCBI] Query: {query}")
    handle = ncbi_request_with_retry(Entrez.esearch, db="nucleotide", term=query, retmax=20)
    if not handle:
        return []
    search_results = Entrez.read(handle)
    handle.close()
    ids = search_results.get("IdList", [])
    if not ids:
        return []
    handle = ncbi_request_with_retry(Entrez.esummary, db="nucleotide", id=','.join(ids))
    summaries = Entrez.read(handle)
    handle.close()
    results = []
    for s in summaries:
        taxid = s.get("TaxId")
        acc = s.get("AccessionVersion", "N/A")
        species = s.get("Organism", s.get("Title", "Unknown species").split(" ")[0])
        results.append({
            "accession": acc,
            "species": species,
            "title": s.get("Title", "N/A"),
            "taxid": str(taxid) if taxid else "0"
        })
    return results

def select_best_candidate_with_ai(candidates, target_species):
    print("\n[AI] Starting semantic relevance analysis using BioBERT...")
    target_embedding = embed_text(f"{target_species} phylogeny organelle genome")
    if target_embedding is None:
        print("  [AI Fallback] Could not embed target species. Returning None.")
        return None
    candidate_scores = []
    for cand in candidates:
        abstract = fetch_pubmed_abstract(cand["species"])
        if not abstract:
            continue
        emb = embed_text(abstract)
        if emb is None:
            continue
        score = cosine_similarity([target_embedding], [emb])[0][0]
        print(f"    {cand['species']} (Acc: {cand['accession']}): Similarity = {score:.3f}")
        candidate_scores.append((cand, score))
    if not candidate_scores:
        return None
    best = max(candidate_scores, key=lambda x: x[1])
    if best[1] < MIN_SIMILARITY_THRESHOLD:
        return None
    print(f"\n[AI] Selected: {best[0]['species']} (Accession: {best[0]['accession']}) with score {best[1]:.3f}")
    return best[0]

def main():
    parser = argparse.ArgumentParser(description="Organoverse AI Genome Selector")
    parser.add_argument("--species", type=str, required=True, help="Target species name")
    parser.add_argument("--organelle", type=str, required=True, choices=["chloroplast", "mitochondrion"], help="Organelle type")
    args = parser.parse_args()
    print("\n=== Organoverse AI Genome Selector (v2.3) ===")
    species = args.species.strip()
    organelle = args.organelle.strip().lower()

    lineage = fetch_taxonomic_lineage(species)
    if not lineage:
        print("[Error] Failed to fetch taxonomic lineage.")
        return

    selected = None
    for rank, taxon in lineage:
        print(f"\n--- Trying at {rank.upper()}: {taxon} ---")
        candidates = search_refseq_organelle(taxon, organelle)
        if candidates:
            selected = select_best_candidate_with_ai(candidates, species)
            if selected:
                break
            else:
                print("[Info] AI did not find a suitable candidate. Trying next rank...")

    if not selected:
        print("[Failure] No suitable organelle genome found.")
        return

    print("\n--- FINAL SELECTION ---")
    print(f"Species: {selected['species']}")
    print(f"Accession: {selected['accession']}")
    print(f"Title: {selected['title']}")
    print(f"TaxID: {selected['taxid']}")

    fasta = fetch_fasta(selected['accession'])
    if fasta:
        filename = f"{species.replace(' ', '_')}_{organelle}_{selected['accession']}.fasta"
        with open(filename, "w") as f:
            f.write(fasta)
        print(f"[Output] FASTA saved to: {filename}")
    else:
        print("[Warning] FASTA download failed.")

    save_summary_to_csv(species, organelle, selected)

if __name__ == "__main__":
    main()
# 🌍 Organoverse

**Organoverse** is a fully automated, AI-enhanced pipeline for the assembly, annotation, and phylogenetic analysis of organelle genomes (chloroplast and mitochondrial) directly from whole-genome sequencing (WGS) data.  
It integrates machine learning models, curated databases, and modern assembly/annotation tools to ensure accurate and reproducible results.

---

## 🧬 Key Features

- 🔬 **K-mer & ML-Based Organelle Read Extraction**  
  Combines taxon-aware k-mer analysis and Random Forest models to enrich for true organelle reads and remove nuclear contaminants (NUMTs/NUPTs).

- 🧠 **AI-Enhanced Annotation & Validation**  
  Uses Prokka and BLAST, followed by AI-based checks for codon completeness, internal stops, and gene frame integrity.

- 🧪 **De Novo Assembly & Circularization**  
  Assembles organelle genomes using SPAdes, then applies overlap-based circularization logic.

- 🌿 **Phylogenomic Tree Construction**  
  Detects conserved genes across taxa, aligns them with MAFFT, and infers phylogeny using IQ-TREE.

- 📊 **Interactive Reports (HTML + PDF)**  
  Generates complete, user-friendly reports including genome stats, feature tables, and trees.

- ⚙️ **Workflow Automation**  
  Available in both [Snakemake](https://snakemake.readthedocs.io) and [Nextflow](https://www.nextflow.io) with full Conda/Bioconda support.

---

## 🏗️ Installation

### 1. Clone the repository

```bash
git clone https://github.com/your-username/organoverse.git
cd organoverse
```
### 2.  Create the Conda environment
```bash
conda env create -f envs/organoverse.yaml
conda activate organoverse
```
### 3. Download test data
```bash
make test-data
```
## 🚀 Usage
### 1. Using Snakemake
```bash
snakemake --cores 4
```
### 2. Using Nextflow
```bash
nextflow run workflows/main.nf
```
Configure species and paths in config/config.yaml.

### 3. Direct CLI
An all-in-one script is provided for smaller jobs and testing:
```bash
python scripts/organoverse.py --species "Arabidopsis thaliana" \
    --organelle chloroplast --library paired \
    --reads R1.fq.gz R2.fq.gz --outdir results
```
## 🧪 Test Datasets
Sample mitochondrial and chloroplast genomes from:

- Arabidopsis thaliana (cp, mt)

- Olea europaea (mt)

- Homo sapiens (mt)

Provided in /test/ after running make test-data.
## 📂 Output
results/circularized/*.fasta → final genome

results/annotation/*.gbk → annotated GenBank file

results/report/*.html and *.pdf → report
## 🧠 Citation
If you use Organoverse in your research, please cite:

Idrissi Azami, A., et al. (2025). Organoverse: An AI-powered pipeline for organelle genome assembly and annotation. Bioinformatics Tools and Applications. (in prep)
## 🤝 Contributors
Dr. Abdellah Idrissi Azami – LinkedIn | GitHub
Special thanks to contributors from Organoverse Project Team

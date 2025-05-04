
import sys
from jinja2 import Environment, FileSystemLoader
from Bio import SeqIO

def generate_html_report(fasta_file, annotation_file, output_html):
    # Charger séquences et infos d’annotation
    fasta_record = next(SeqIO.parse(fasta_file, "fasta"))
    annotation_record = next(SeqIO.parse(annotation_file, "genbank"))
    
    genome_length = len(fasta_record.seq)
    gene_count = sum(1 for feature in annotation_record.features if feature.type == "CDS")
    rna_count = sum(1 for feature in annotation_record.features if "RNA" in feature.type.upper())
    
    # Environnement Jinja2
    env = Environment(loader=FileSystemLoader("templates"))
    template = env.get_template("report_template.html")

    html_content = template.render(
        fasta_file=fasta_file,
        annotation_file=annotation_file,
        genome_length=genome_length,
        gene_count=gene_count,
        rna_count=rna_count,
    )

    with open(output_html, "w") as f:
        f.write(html_content)
    print(f"[✔] Report generated: {output_html}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python generate_report.py <FASTA> <GBK> <OUTPUT_HTML>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    annotation_file = sys.argv[2]
    output_html = sys.argv[3]

    generate_html_report(fasta_file, annotation_file, output_html)

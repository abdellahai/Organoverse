nextflow.enable.dsl=2

params.config = file("config/config.yaml")
workflow {

    reads = Channel.fromPath(params.config.input_reads.R1)
    reads2 = Channel.fromPath(params.config.input_reads.R2)

    cleaned_reads = quality_control(reads, reads2)
    filtered_reads = ai_filtering(cleaned_reads[0], cleaned_reads[1])
    contigs = assemble(filtered_reads[0], filtered_reads[1])
    circular = circularize(contigs)
    annotated = annotate(circular)
    report(annotated, circular)
}

process quality_control {

    input:
    path r1
    path r2

    output:
    path "results/cleaned_reads/R1_trimmed.fq.gz"
    path "results/cleaned_reads/R2_trimmed.fq.gz"

    script:
    """
    fastp -i ${r1} -I ${r2} \
          -o results/cleaned_reads/R1_trimmed.fq.gz \
          -O results/cleaned_reads/R2_trimmed.fq.gz \
          --html results/qc_reports/fastp.html
    """
}

process ai_filtering {

    input:
    path r1
    path r2

    output:
    path "results/verified_reads/R1.fq"
    path "results/verified_reads/R2.fq"

    script:
    """
    python scripts/ai_filtering.py ${r1} ${r2} \
           results/verified_reads/R1.fq \
           results/verified_reads/R2.fq \
           ${params.config.ml_model}
    """
}

process assemble {

    input:
    path r1
    path r2

    output:
    path "results/assembly/contigs.fasta"

    script:
    """
    spades.py -1 ${r1} -2 ${r2} -o results/assembly --${params.config.spades_mode}
    cp results/assembly/contigs.fasta results/assembly/contigs.fasta
    """
}

process circularize {

    input:
    path contigs

    output:
    path "results/circularized/organelle_circular.fasta"

    script:
    """
    python scripts/circularize.py ${contigs} results/circularized/organelle_circular.fasta
    """
}

process annotate {

    input:
    path fasta

    output:
    path "results/annotation/organelle.gbk"

    script:
    """
    prokka --outdir results/annotation \
           --prefix organelle \
           --kingdom Bacteria \
           --genus ${params.config.species.split(" ")[0]} \
           --species "${params.config.species}" \
           --gcode ${params.config.translation_table} \
           ${fasta}
    """
}

process report {

    input:
    path annotation
    path fasta

    output:
    path "results/report/final_report.html"
    path "results/report/final_report.pdf"

    script:
    """
    python scripts/generate_report.py ${fasta} ${annotation} results/report/final_report.html
    python scripts/generate_pdf_report.py results/report/final_report.html results/report/final_report.pdf
    """
}

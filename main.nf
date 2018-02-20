#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    ================================================================
     guidemapper-nf
    ================================================================
    DESCRIPTION

    Usage:
    nextflow run zuberlab/guidemapper-nf --genome <genome_name>

        <genome_name>       Name of genome to align to. Corresponding GTF file,
                            bowtie2 index, bioconductor organism package, and
                            type of 'gene_id' in GTF file. Must be specified
                            in config file.

    Options:
        --inputDir          Directory containing input files (default: 'input').
                            Input files must contain at least the columns 'id',
                            'group', and 'sequence'.
        --resultsDir        Directory to write output to (default: 'results').
        --genomeDir         Base directory containing genome annotation files
                            referenced in <genome_name>.
        --agnostic          Annotate guides without using the 'group' columns
                            in input file (default: false).
        --legacy_feature    Type of annotation in 'group' column in input files.
                            Must be one of 'SYMBOL', 'ENSEMBL', or 'ENTREZID'
                            (default: 'SYMBOL').
        --max_alignments    Maximum number of perfect alignment to report
                            (default: 10). Increasing this number could
                            considerably slow down the alignment process in
                            the presence of unspecific guides and large genomes.
    Profiles:
        standard            local execution
        singularity         local execution with singularity
        ii2                 SLURM execution with singularity on IMPIMBA2

    Docker:
    zuberlab/guidemapper-nf:latest

    Author:
    Jesse J. Lipp (jesse.lipp@imp.ac.at)
    """.stripIndent()
}

params.help = false

if (params.help){
    helpMessage()
    exit 0
}

params.bt2 = params.genome ? params.genomes[ params.genome ].bt2 ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.org_package = params.genome ? params.genomes[ params.genome ].org_package ?: false : false
params.gtf_feature = params.genome ? params.genomes[ params.genome ].gtf_feature ?: false : false

log.info ""
log.info " parameters "
log.info " ======================"
log.info " input directory         : ${params.inputDir}"
log.info " result directory        : ${params.resultsDir}"
log.info " genome base directory   : ${params.genomeDir}"
log.info " genome                  : ${params.genome}"
log.info " GTF file                : ${params.gtf}"
log.info " Bowtie2 index           : ${params.bt2}"
log.info " agnostic annotation     : ${params.agnostic}"
log.info " legacy feature type     : ${params.legacy_feature}"
log.info " GTF feature type        : ${params.gtf_feature}"
log.info " max reported alignments : ${params.max_alignments}"
log.info " ======================"
log.info ""

Channel
    .fromPath( "${params.inputDir}/*.txt" )
    .map { file -> tuple( file.baseName, file ) }
    .into { libInputFiles; libAnnotationFiles }

Channel
    .fromPath( params.gtf )
    .ifEmpty { exit 1, "GTF file not found: ${params.gtf}" }
    .set { gtfFile }

Channel
    .fromPath( "${params.bt2}*.bt2" )
    .ifEmpty { exit 1, "Bowtie2 index not found: ${params.bt2}" }
    .toList()
    .map { files -> tuple( file(params.bt2).baseName, files) }
    .set { bt2Index }

process input2fasta {

    tag { name }

    input:
    set val(name), file(library) from libInputFiles

    output:
    set val(name), file("${name}.fa") into libFastaFiles

    script:
    """
    input2fasta.R ${library} ${name}
    """
}

process align {

    tag { name }

    publishDir path: "${params.resultsDir}/${name}",
               mode: 'copy',
               overwrite: 'true',
               saveAs: {filename ->
                   if (filename.indexOf(".log") > 0) "${filename}"
                   else if (filename.indexOf(".bam") > 0) "${filename}"
                   else null}

    input:
    set val(name), file(fasta) from libFastaFiles
    set val(prefix), file(index) from bt2Index

    output:
    set val(name), file("${name}.bam") into alignedBamFiles
    file "${name}_bt2.log" into alignLogs

    script:
    """
    bowtie2 \
        --threads 1 \
        -x ${prefix} \
        -f ${fasta} \
        -L 18 \
        --score-min 'C,0,-1' \
        -k ${params.max_alignments} \
        -S ${name}.sam 2> ${name}_bt2.log
    samtools view -bS ${name}.sam > ${name}.bam
    """
}

libAnnotationFiles
    .concat(alignedBamFiles)
    .groupTuple()
    .set{ annotationFiles }

process annotate_features {

    tag { name }

    publishDir path: "${params.resultsDir}/${name}",
               mode: 'copy',
               overwrite: 'true'

    input:
    set val(name), file(files) from annotationFiles
    each file(gtf) from gtfFile

    output:
    set val(name), file("${name}_annotation.txt") into resultsFiles

    script:
    """
    annotate_features.R \
        ${files[0]} \
        ${files[1]} \
        ${gtf} \
        ${params.org_package} \
        ${params.agnostic} \
        ${params.legacy_feature} \
        ${params.gtf_feature} \
        ${name}
    """
}

workflow.onComplete {
	println ( workflow.success ? "COMPLETED!" : "FAILED" )
}

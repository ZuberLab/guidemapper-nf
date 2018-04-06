#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    ================================================================
     guidemapper-nf
    ================================================================
    DESCRIPTION

    Usage:
    nextflow run zuberlab/guidemapper-nf

        <genome_name>       Name of genome to align to. Corresponding GTF file,
                            bowtie2 index, bioconductor organism package, and
                            bioconductor genome package. Must be specified
                            in config file.

    Options:
        --parameters        Tab-separated file (default: 'parameters.txt')
                            containing with two columns:
                            library_id (must match input file) and genome
                            (must be defined in config in 'genomes')
        --inputDir          Directory containing input files (default: 'input').
                            Input files must contain at least the columns 'id',
                            'group', and 'sequence'.
        --resultsDir        Directory to write output to (default: 'results').
        --genomeDir         Base directory containing genome annotation files
                            referenced in <genome_name>.
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

log.info ""
log.info " parameters "
log.info " ======================"
log.info " parameter file          : ${params.parameters}"
log.info " input directory         : ${params.inputDir}"
log.info " result directory        : ${params.resultsDir}"
log.info " genome base directory   : ${params.genomeDir}"
log.info " max reported alignments : ${params.max_alignments}"
log.info " ======================"
log.info ""

Channel
    .fromPath( params.parameters )
    .splitCsv(sep: '\t', header: true)
    .map { it -> tuple(
        it.library_id, [
            name : it.library_id,
            bt2 : params.genomes[ it.genome ].bt2,
            org_package : params.genomes[ it.genome ].org_package,
            genome_package : params.genomes[ it.genome ].genome_package,
            gtf : file(params.genomes[ it.genome ].gtf) ] ) }
    .set { libraryParameters }

Channel
    .fromPath( "${params.inputDir}/*.txt" )
    .map { file -> tuple( file.baseName, file ) }
    .into { libraryFiles; libraryGrouping }

libraryParameters
    .join( libraryFiles )
    .set { libraryInputs }

process input2fasta {

    tag { name }

    input:
    set val(name), val(parameters), file(seqs) from libraryInputs

    output:
    set val(name), val(parameters), file("${name}.fa") into libraryFasta

    script:
    """
    input2fasta.R ${seqs} ${name}
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
    set val(name), val(parameters), file(fasta) from libraryFasta

    output:
    set val(name), val(parameters), file("${name}.bam") into libraryAligned
    file "${name}_bt2.log" into alignLogs

    script:
    """
    bowtie2 \
        --threads 1 \
        -x ${parameters.bt2} \
        -f ${fasta} \
        -L 18 \
        --score-min 'C,0,-1' \
        -k ${params.max_alignments} \
        -S ${name}.sam 2> ${name}_bt2.log
    samtools view -bS ${name}.sam > ${name}.bam
    """
}

libraryAligned
    .join( libraryGrouping )
    .set{ libraryAnnotation }


process annotate_features {

    tag { name }

    publishDir path: "${params.resultsDir}/${name}",
               mode: 'copy',
               overwrite: 'true'

    input:
    set val(name), val(parameters), file(bam), file(seqs) from libraryAnnotation

    output:
    set val(name), file("${name}_annotation.txt") into resultsFiles

    script:
    """
    annotate_features.R \
        ${seqs} \
        ${bam} \
        ${parameters.gtf} \
        ${parameters.org_package} \
        ${parameters.genome_package} > ${name}_annotation.txt
    """
}

workflow.onComplete {
	println ( workflow.success ? "COMPLETED!" : "FAILED" )
}

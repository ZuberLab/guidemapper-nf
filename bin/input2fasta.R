#!/usr/bin/env Rscript

################################################################################
# convert library input file to fasta for alignment
# part of guidemapper-nf pipeline
# github: zuberlab/guidemapper-nf
# 
# Jesse J. Lipp
# Institute of Molecular Pathology (IMP), Vienna, Austria
# started 2018/02/16
################################################################################

### functions
`%>%` <- dplyr::`%>%`

### command line parameters
args     <- commandArgs(trailingOnly = TRUE)
file_lib <- args[1]
name_lib <- args[2]

### convert to fasta
readr::read_tsv(file_lib) %>%
  dplyr::select(id, sequence) %>%
  tibble::deframe() %>%
  Biostrings::DNAStringSet() %>%
  Biostrings::writeXStringSet(paste0(name_lib, ".fa"))
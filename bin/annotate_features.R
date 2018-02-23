#!/usr/bin/env Rscript

################################################################################
# annotate library using original input file, GTF file, and aligned BAM file
# part of guidemapper-nf pipeline
# github: zuberlab/guidemapper-nf
# 
# Jesse J. Lipp
# Institute of Molecular Pathology (IMP), Vienna, Austria
# started 2018/02/16
################################################################################

### packages
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(GenomicAlignments)
library(GenomicFeatures)
library(Biostrings)
library(tidyverse)
library(stringr)

### command line parameters
args <- commandArgs(trailingOnly = TRUE)
file_lib <- args[1]
file_bam <- args[2]
file_gtf <- args[3]
org_db <- args[4]
annotate_agnostic <- as.logical(args[5])
legacy_feature_type <- args[6]
feature_type <- args[7]
name_lib <- args[8]

stopifnot(legacy_feature_type %in% c("ENTREZID", "SYMBOL", "ENSEMBL"))
stopifnot(feature_type %in% c("ENTREZID", "SYMBOL", "ENSEMBL"))
stopifnot(org_db %in% c("org.Hs.eg.db", "org.Mm.eg.db"))

### import
seq <- read_tsv(file_lib) %>%
  rename(legacy_feature = group) %>%
  mutate(legacy_feature = as.character(legacy_feature))

# set agnostic to TRUE if all legacy features are NA
if (all(is.na(seq$legacy_feature))) {
  annotate_agnostic <- TRUE
}

bam <- readGAlignments(file = file_bam, use.names = TRUE)
txdb <- makeTxDbFromGFF(file = file_gtf)

cds <- cdsBy(txdb, by = "gene")
txs <- transcriptsBy(txdb, by = "gene")

cds_ov <- findOverlaps(bam, cds, ignore.strand = TRUE)
txs_ov <- findOverlaps(bam, txs, ignore.strand = TRUE)

### alignment
genome <- tibble(id = names(bam), 
                 chromosome = as.character(seqnames(bam)), 
                 strand = as.character(strand(bam)), 
                 start = start(bam), 
                 end = end(bam)) %>%
  group_by(id) %>%
  mutate(pm_genome = n(),
         chromosome = if_else(pm_genome > 1, NA_character_, chromosome),
         strand = if_else(pm_genome > 1, NA_character_, strand),
         start = if_else(pm_genome > 1, NA_integer_, start),
         end = if_else(pm_genome > 1, NA_integer_, end)) %>%
  ungroup %>%
  distinct

### target features
features_cds <- tibble(id = names(bam)[queryHits(cds_ov)], 
                       target_cds = names(cds)[subjectHits(cds_ov)]) %>%
  group_by(id) %>%
  summarize(target_cds = paste(sort(target_cds), collapse = "|"), 
            pm_cds = n()) %>%
  ungroup

features_txs <- tibble(id = names(bam)[queryHits(txs_ov)], 
                       target_txs = names(txs)[subjectHits(txs_ov)]) %>%
  group_by(id) %>%
  summarize(target_txs = paste(sort(target_txs), collapse = "|"),
            pm_txs = n()) %>%
  ungroup

temp <- seq %>%
  left_join(genome) %>%
  left_join(features_cds) %>%
  left_join(features_txs)

if (!annotate_agnostic & !feature_type == legacy_feature_type) {
  features <- temp %>%
    mutate(feature = coalesce(target_cds, target_txs)) %>%
    select(id, legacy_feature, feature) %>%
    mutate(feature = str_split(feature, pattern = "[|]")) %>%
    unnest(feature) %>%
    mutate(legacy_feature_mapped = mapIds(eval(parse(text = org_db)),
                                          keys = legacy_feature,
                                          keytype = legacy_feature_type,
                                          column = feature_type)) %>%
    filter(feature == legacy_feature_mapped) %>%
    distinct(id, feature)
} else if (!annotate_agnostic & feature_type == legacy_feature_type) {
  features <- temp %>%
    mutate(feature = coalesce(target_cds, target_txs)) %>%
    select(id, legacy_feature, feature) %>%
    mutate(feature = str_split(feature, pattern = "[|]")) %>%
    unnest(feature) %>%
    filter(feature == legacy_feature) %>%
    distinct(id, feature)
} else {
  features <- temp %>%
    mutate(feature = coalesce(target_cds, target_txs)) %>%
    select(id, feature)
}

results <- temp %>%
  left_join(features) %>%
  mutate(feature = coalesce(feature, target_cds, target_txs)) %>%
  mutate(type_guide = case_when(!is.na(target_cds) ~ "coding",
                                !is.na(target_txs) ~ "non_coding",
                                TRUE ~ "unknown")) %>%
  mutate(type_alignment = case_when(pm_genome == 1 ~ "unique",
                                    pm_genome > 1 ~ "ambiguous",
                                    TRUE ~ "not_mapped")) %>%
  mutate(type_feature = case_when(pm_cds == 1 ~ "unique", 
                                  pm_cds > 1 ~ "ambiguous", 
                                  TRUE ~ "unknown")) %>%
  select(id, 
         feature, 
         legacy_feature, 
         sequence, 
         chromosome, 
         strand, 
         start, 
         end, 
         pm_genome, 
         pm_cds, 
         pm_txs, 
         type_guide, 
         type_alignment, 
         type_feature, 
         target_cds, 
         target_txs)

results %>%
  write_tsv(paste0(name_lib, "_annotation.txt"))

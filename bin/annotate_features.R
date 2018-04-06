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
library(BSgenome)
library(GenomicAlignments)
library(GenomicFeatures)
library(Biostrings)
library(tidyverse)
library(stringr)
library(rlang)

### parameters
args <- commandArgs(trailingOnly = TRUE)
file_guides <- args[1]
file_bam <- args[2]
file_gtf <- args[3]
org_db <- args[4]
genome_db <- args[5]

### debug
# file_guides <- "test_input/test_avana.txt"
# file_bam <- "results/test_avana/test_avana.bam"
# file_gtf <- "/Volumes/groups/zuber/zubarchive/GROUPS/Bioinfo/genomes/GRCh38/gtf/Homo_sapiens.GRCh38.91.gtf"
# org_db <- "org.Hs.eg.db"
# genome_db <- "BSgenome.Hsapiens.UCSC.hg38"

### functions
read_guides <- function(file) {
  read_tsv(file, col_types = "ccc") %>%
    dplyr::rename(legacy_id = id, legacy_group = group)
}

get_genome_object <- function(genome_db) {
  object <- genome_db %>%
    str_split(pattern = "[.]") %>%
    unlist %>%
    .[2]
  paste(genome_db, object, sep = "::") %>%
    rlang::parse_quosure()
}

get_org_object <- function(org_db) {
  rlang::parse_quosure(paste(rep(org_db, 2), collapse = "::"))
}

get_context_ranges <- function(gr) {
  gr %>%
    resize(width = 24, fix = "end") %>%
    resize(width = 30, fix = "start")
}

get_contexts <- function(genome_db, gr) {
  Views(eval_tidy(get_genome_object(genome_db)), gr)
}

get_ngg <- function(contexts) {
  str_detect(str_sub(contexts, start = 25, end = 27), "[ATGC]GG")
}

# get_out_of_bounds <- function(gr) {
#   seq(gr) %in% GenomicRanges:::get_out_of_bound_index(gr)
# }

clean_mapping <- function(x) {
  ifelse(is.null(x), NA_character_, as.character(x))
}

map_redundantly <- function(df, org_db) {
  df %>%
    mutate(entrez_id = mapIds(eval_tidy(get_org_object(org_db)),
                                keys = ensembl_id,
                                keytype = "ENSEMBL",
                                column = "ENTREZID", 
                                multiVals = "list"), 
           entrez_id = map_chr(entrez_id, clean_mapping)) %>%
    mutate(hgnc_symbol = mapIds(eval_tidy(get_org_object(org_db)),
                              keys = ensembl_id,
                              keytype = "ENSEMBL",
                              column = "SYMBOL", 
                              multiVals = "list"), 
           hgnc_symbol = map_chr(hgnc_symbol, clean_mapping))
}

get_match_numbers <- function(df) {
  df %>%
    group_by(legacy_id) %>%
    mutate(pm_total = n(), 
           pm_pam = sum(has_ngg), 
           pm_pam_cds = sum(!is.na(ensembl_id) & has_ngg), 
           pm_pam_cds_unique = length(unique(ensembl_id[!is.na(ensembl_id) & has_ngg]))) %>%
    ungroup
}

replace_empty_with_na <- function(x) {
  if_else(x == "", NA_character_, x)
}

summarize_features <- function(df) {
  df %>%
    group_by(legacy_id) %>%
    # remove alignments outside CDS, keep coordinates and annotation for cases
    # where there is a single cds match but other perfect PAM matches in the genome
    filter(!(is.na(ensembl_id) & n() > 1 & pm_pam_cds > 0)) %>%
    summarize(pm_total = unique(pm_total),
              pm_pam = unique(pm_pam),
              pm_pam_cds = unique(pm_pam_cds),
              pm_pam_cds_unique = unique(pm_pam_cds_unique), 
              entrez_id = paste(sort(unique(entrez_id)), collapse = "|"), 
              ensembl_id = paste(sort(unique(ensembl_id)), collapse = "|"), 
              hgnc_symbol = paste(sort(unique(hgnc_symbol)), collapse = "|"), 
              # use ifelse instead of if_else, otherwise NA's are not silently replicated
              chromosome = ifelse(n() > 1, NA_character_, chromosome),
              strand = ifelse(n() > 1, NA_character_, strand),
              start = ifelse(n() > 1, NA_integer_, start),
              end = ifelse(n() > 1, NA_integer_, end), 
              context = ifelse(n() == 1, context, NA_character_)) %>%
    ungroup %>%
    # replace empty with NA
    mutate_at(vars(hgnc_symbol, entrez_id, ensembl_id), replace_empty_with_na)
}

assign_groups <- function(df) {
  df %>%
    mutate(group = case_when(legacy_group %in% c("NONTARGETING", "SAFETARGETING") ~ legacy_group, 
                             is.na(pm_total) ~ "UNMAPPED", 
                             is.na(ensembl_id) ~ "NOFEATURE", 
                             str_detect(ensembl_id, "[|]") ~ "AMBIGUOUS",
                             TRUE ~ ensembl_id))
}

classify_guides <- function(df) {
  df %>%
    mutate(class = case_when(legacy_group == "NONTARGETING" ~ "control_nontargeting", 
                             legacy_group == "SAFETARGETING" ~ "control_safetargeting", 
                             pm_total == 1 & pm_pam == 1 & pm_pam_cds == 1 & pm_pam_cds_unique == 1 ~ "ok", 
                             pm_total > 1 & pm_pam == 1 & pm_pam_cds == 1 & pm_pam_cds_unique == 1 ~ "multi_no_pam",
                             pm_total > 1 & pm_pam > 1 & pm_pam_cds == 1 & pm_pam_cds_unique == 1 ~ "multi_pam",
                             pm_total > 1 & pm_pam > 1 & pm_pam_cds > 1 & pm_pam_cds_unique == 1 ~ "multi_same_cds",
                             pm_total > 1 & pm_pam > 1 & pm_pam_cds > 1 & pm_pam_cds_unique > 1 ~ "multi_cds",
                             is.na(pm_total) ~ "not_mapped", 
                             is.na(ensembl_id) ~ "no_feature", 
                             TRUE ~ "unknown"
    ))
}

generate_guide_id <- function(df) {
  df %>%
    mutate(id = paste(group, sequence, sep = "_"))
}

format_output <- function(df) {
  df %>%
    dplyr::select(id, group, sequence, hgnc_symbol, entrez_id, ensembl_id, 
                  chromosome, strand, start, end, 
                  pm_total, pm_pam, pm_pam_cds, pm_pam_cds_unique, context, 
                  legacy_id, legacy_group, class) %>%
    arrange(ensembl_id)
}

### import
guides <- read_guides(file_guides)
bam <- readGAlignments(file_bam, use.names = TRUE)
txdb <- makeTxDbFromGFF(file_gtf)

### processing

# contexts
bam_ucsc <- bam
seqlevelsStyle(bam_ucsc) <- "UCSC"
gr_context <- get_context_ranges(GRanges(bam_ucsc))
contexts <- get_contexts(genome_db, gr_context)

# features
entrez <- keys(eval_tidy(get_org_object(org_db)), keytype = "ENTREZID")
ensembl <- mapIds(eval_tidy(get_org_object(org_db)), 
                            keys = entrez, 
                            keytype = "ENTREZID", 
                            column = "ENSEMBL")
ensembl <- ensembl[!is.na(ensembl)]

cds <- cdsBy(txdb, by = "gene")
cds <- cds[names(cds) %in% ensembl]
cds_ov <- findOverlaps(bam, cds, ignore.strand = TRUE)

features <- tibble(legacy_id = names(bam), 
                   ensembl_id = NA_character_, 
                   chromosome = as.character(seqnames(bam)),
                   strand = as.character(strand(bam)),
                   start = start(bam),
                   end = end(bam), 
                   context = as.character(contexts), 
                   has_ngg = get_ngg(contexts))
# queryHits(cds_ov) is NOT unique.
# one alignment can match multiple (overlapping) cds
# this information is lost in this assignment in favor of the latest feature
features$ensembl_id[queryHits(cds_ov)] <- names(cds)[subjectHits(cds_ov)]

annotation <- features %>%
  # add perfect matches to features
  get_match_numbers %>%
  # add gene mapping to features
  filter(has_ngg) %>%
  map_redundantly(org_db) %>%
  # summarize features
  summarize_features %>%
  # add original library information
  right_join(guides) %>%
  # assign groups
  assign_groups %>%
  # classify guides
  classify_guides %>%
  # clean up
  generate_guide_id %>%
  format_output

# export
annotation %>%
  format_tsv %>%
  cat
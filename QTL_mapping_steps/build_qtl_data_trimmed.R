library(reticulate)
library(dplyr)
library(stringr)
library(biomaRt)
library(httr2)
library(tidyverse)

base_dir <- "/Users/reillj/Documents/QTL/results"
llm_data_dir <- "/Users/reillj/Documents/LLM/py_version/data"

phenos <- c("LD", "OFA", "HB", "NPP")

qtl_data_trimmed <- list()

for (p in phenos) {
  
  message("Building QTL object for group: ", p)
  
  # ---------------- LOAD TRIMMED COEFFICIENTS ----------------
  coef_file <- file.path(
    base_dir,
    paste0("HPC_filter_output_", p),
    paste0(p, "_filter_combined"),
    paste0(p, "_combined_coef_trim.rds")
  )
  
  genes_file <- file.path(
    base_dir,
    paste0("HPC_filter_output_", p),
    paste0(p, "_filter_combined"),
    paste0(p, "_combined_genes.rds")
  )
  
  snps_file <- file.path(
    base_dir,
    paste0("HPC_filter_output_", p),
    paste0(p, "_filter_combined"),
    paste0(p, "_combined_scan1snps.rds")
  )
  
  scan1_file <- file.path(
    base_dir,
    paste0("HPC_filter_output_", p),
    paste0(p, "_filter_combined"),
    paste0(p, "_filter_combined_scan1_out.rds")
  )
  
  go_file <- file.path(
    base_dir,
    paste0("HPC_filter_output_", p),
    paste0(p, "_filter_combined"),
    paste0(p, "_combined_GO.rds")
  )
  
  if (!file.exists(coef_file)) {
    stop("Missing trimmed coef file for group: ", p)
  }
  
  genes <- readRDS(genes_file)
  coef_trimmed <- readRDS(coef_file)
  snps <- readRDS(snps_file)
  scan1 <- readRDS(scan1_file)
  go <- readRDS(go_file)
  
  message("Cleaning gene table for group: ", p)
  
  # Remove unused columns
  genes <- genes %>%
    dplyr::select(-any_of(c("Parent", "score", "phase")))
  
  # ---------------- REMOVE PREDICTED / PSEUDOGENES ----------------
  
  message("Filtering predicted genes in group: ", p)
  
  predicted_genes <- genes %>%
    filter(str_detect(tolower(description), "^(predicted|pseudogene)")) %>%
    pull(Name)
  genes <- genes %>%
    filter(!Name %in% predicted_genes)
  
  message(length(predicted_genes), " predicted/pseudogenes removed.")
  
  # ---------------- REMOVE RIK GENES ----------------
  
  message("Filtering Rik genes in group: ", p)
  
  rik_genes <- genes %>%
    filter(str_detect(Name, "Rik")) %>%
    pull(Name)
  
  genes <- genes %>%
    filter(!Name %in% rik_genes)
  
  message(length(rik_genes), " Rik genes removed.")
  
  # ---------------- REMOVE MISSING gene_id ----------------
  
  before_na <- nrow(genes)
  genes <- genes %>%
    filter(!is.na(gene_id) & gene_id != "")
  na_removed <- before_na - nrow(genes)
  
  message(na_removed, " entries removed due to missing gene_id.")
  
  # ---------------- FILTER GO TERMS BY Q-VALUE ----------------
  
  message("Filtering GO terms by q-value in group: ", p)
  
  before_go <- nrow(go)
  
  go <- go %>%
    filter(!is.na(qvalue) & qvalue <= 0.01)
  
  after_go <- nrow(go)
  
  message(before_go - after_go, " GO terms removed (q-value > 0.01).")
  
  # --------------- REMOVE snps with LOD < 5 ----------------
  
  before_LOD_filter <- nrow(snps)
  
  snps <- snps %>%
    group_by(index, interval) %>%
    filter(any(lod > 5, na.rm = TRUE)) %>%
    ungroup()
  
  removed <- before_LOD_filter - nrow(snps)
  
  message(removed, " SNPs removed due to low LOD groups.")
  
  # ---------------- BUILD GROUP ENTRY ----------------
  qtl_data_trimmed[[p]] <- list(
    genes  = genes,
    snps   = snps,
    coef   = coef_trimmed,   
    scan1  = scan1,
    go     = go
  )
}

groups <- names(qtl_data_trimmed)

for (g in groups) {
  
  message("Renaming founder columns in group: ", g)
  
  df <- qtl_data_trimmed[[g]]$coef
  
  df_renamed <- df %>%
    rename("A_J" = "A",
           "C57BL_6J" = "B",
           "129S1_SvImJ" = "C",
           "NOD_ShiLtJ" = "D",
           "NZO_HlLtJ" = "E",
           "CAST_EiJ" = "F",
           "PWK_PhJ" = "G",
           "WSB_EiJ" = "H")
  
  qtl_data_trimmed[[g]]$coef <- df_renamed
}

message("✅ Founder strain columns renamed successfully.")

# MAP ENSEMBL IDs to GENE NAMES

mappings <- read_csv("/Users/reillj/Documents/QTL/data/histoUCF_Grcm38_map.csv")
mappings <- dplyr::rename(mappings, "Name" = "name")

groups <- names(qtl_data_trimmed)

for (g in groups) {
  message("Mapping ENSEMBL IDs to gene names in group: ", g)
  
  genes_df <- qtl_data_trimmed[[g]]$genes
  snps_df  <- qtl_data_trimmed[[g]]$snps
  
  # Extract ENSEMBL IDs from Dbxref column
  genes_df <- genes_df %>%
    mutate(ENSEMBL = str_extract(Dbxref, "ENSMUSG[0-9]+")) %>%
    filter(!is.na(ENSEMBL)) %>%
    dplyr::select(ENSEMBL, Name)
  
  # Join to SNP table
  snps_df <- snps_df %>%
    left_join(genes_df, by = c("ensembl_gene" = "ENSEMBL"))
  
  snps_df <- snps_df %>%
    left_join(mappings, by = c("ensembl_gene" = "gene"))
  
  snps_df <- snps_df %>%
    mutate(
      gene_name = if_else(
        !is.na(Name.y) & Name.x == Name.y,  # if they match
        Name.x,                             # keep the matched name
        Name.x                              # otherwise, take Name.x
      )
    ) %>%
    dplyr::select(-Name.x, -Name.y)  # drop the old columns
  
  # Save back
  qtl_data_trimmed[[g]]$snps <- snps_df
}

py_save_object(
  qtl_data_trimmed,
  file.path(llm_data_dir, "qtl_data.pkl")
)

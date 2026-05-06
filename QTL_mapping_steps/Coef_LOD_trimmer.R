library(dplyr)
library(readr)
library(purrr)

# -----------------------------
# LOAD DATA
# -----------------------------
intervals <- read_tsv("/Users/reillj/Documents/QTL/results/all_lod_support_intervals.tsv")

intervals <- intervals %>%
  rename(
    start_mb = start_Mb,
    end_mb   = end_Mb,
    chr      = chr
  ) %>%
  select(group, phenotype, chr, start_mb, end_mb)

intervals <- intervals %>%
  mutate(phenotype = str_remove(phenotype, "^qtl_"))

base_dir <- "/Users/reillj/Documents/QTL/results"

load("/Users/reillj/Documents/QTL/data/gm_DO2816_qc.RData")

# -----------------------------
# FUNCTION TO TRIM ONE PHENOTYPE
# -----------------------------
process_one <- function(group, phenotype, chr, start_mb, end_mb) {
  
  message("\nProcessing: ", group, " | ", phenotype, " | chr", chr)
  
  # Build coefficient file path
  coef_path <- file.path(
    base_dir,
    paste0("HPC_filter_output_", group),
    paste0(group, "_", phenotype),
    paste0("coef_filter_chr", chr, ".rds")
  )
  
  if (!file.exists(coef_path)) {
    message("  ❌ Coef file not found: ", coef_path)
    return(NULL)
  }
  
  coef_obj <- readRDS(coef_path)
  
  # -----------------------------
  # GET MARKER POSITIONS
  # -----------------------------
  chr_map <- gm_DO2816_qc$pmap[[as.character(chr)]]
  
  if (is.null(chr_map)) {
    message("  ❌ Chromosome map not found")
    return(NULL)
  }
  
  marker_df <- tibble(
    marker = names(chr_map),
    Mb = as.numeric(chr_map)
  )
  
  # Markers inside LOD interval
  keep_markers <- marker_df %>%
    filter(Mb >= start_mb & Mb <= end_mb) %>%
    pull(marker)
  
  message("  Markers in interval: ", length(keep_markers))
  
  # -----------------------------
  # FILTER COEFFICIENT OBJECT
  # -----------------------------
  # Case 1: markers are rownames
  if (!is.null(rownames(coef_obj))) {
    coef_trimmed <- coef_obj[rownames(coef_obj) %in% keep_markers, , drop = FALSE]
    
    # Case 2: markers in column
  } else if ("marker" %in% colnames(coef_obj)) {
    coef_trimmed <- coef_obj %>% filter(marker %in% keep_markers)
    
  } else {
    message("  ❌ Could not find marker IDs in coef object")
    return(NULL)
  }
  
  message("  Rows before: ", nrow(coef_obj))
  message("  Rows after:  ", nrow(coef_trimmed))
  
  # -----------------------------
  # SAVE BACK
  # -----------------------------
  out_path <- file.path(
    dirname(coef_path),
    paste0("coef_filter_chr", chr, "_LODtrim.rds")
  )
  
  saveRDS(coef_trimmed, out_path)
  
  message("  ✅ Saved: ", out_path)
}

# -----------------------------
# RUN FOR ALL PHENOTYPES
# -----------------------------
pwalk(
  intervals,
  process_one
)

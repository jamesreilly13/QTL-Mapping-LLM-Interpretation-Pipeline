#!/usr/bin/env Rscript
# convert_to_alleleprobs_singlecore.R

DATA_DIR <- "/projects/csna/colby2025/QTLdata"
INPUT_GENOPROB <- file.path(DATA_DIR, "genotype_probs.rds")
OUTPUT_ALLELEPROB <- file.path(DATA_DIR, "allele_probs.rds")

pkgs <- c("qtl2")
missing_pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing packages: ", paste(missing_pkgs, collapse=", "))
library(qtl2)

message("Loading genotype probabilities...")
genoprobs <- readRDS(INPUT_GENOPROB)

message("Converting to allele probabilities (single-core)...")
t0 <- Sys.time()
# Force single core so we don't fork many large processes
alleleprobs <- genoprob_to_alleleprob(genoprobs, quiet = FALSE, cores = 1)
t1 <- Sys.time()
message("Conversion finished in ", round(as.numeric(difftime(t1,t0,units="secs")),1), " sec")

# Basic sanity checks
if (is.null(alleleprobs) || length(alleleprobs) == 0) {
  stop("Converted alleleprobs is empty; aborting save.")
}

message("Saving allele probs to: ", OUTPUT_ALLELEPROB)
saveRDS(alleleprobs, OUTPUT_ALLELEPROB)
message("Done.")

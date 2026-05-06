#!/usr/bin/env Rscript
# convert_to_alleleprobs_singlecore.R

DATA_DIR <- "/projects/csna/colby2025/QTLdata"
INPUT_ALLELEPROB <- file.path(DATA_DIR, "allele_probs.rds")
OUTPUT_KINSHIP_LOCO <- file.path(DATA_DIR, "kinship_loco.rds")
OUTPUT_KINSHIP_OVL <- file.path(DATA_DIR, "kinship_overall.rds")

pkgs <- c("qtl2")
missing_pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing packages: ", paste(missing_pkgs, collapse=", "))
library(qtl2)

ncores_env <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA)
if (!is.na(ncores_env) && nzchar(ncores_env)) {
  ncores <- suppressWarnings(as.integer(ncores_env))
} else {
  ncores <- parallel::detectCores(logical = TRUE)
}
if (is.na(ncores) || ncores < 1) ncores <- 1L
# cap to detectCores just in case
max_cores <- parallel::detectCores(logical = TRUE)
if (!is.na(max_cores) && ncores > max_cores) ncores <- max_cores
message("Using ", ncores, " cores (detected / requested)")

options(mc.cores = ncores)

message("Loading allele probabilities...")
genoprobs <- readRDS(INPUT_ALLELEPROB)

message("Calculating Kinship…")
t0 <- Sys.time()
message("Calculating kinship matrices (loco and overall)...")
k <- calc_kinship(genoprobs, "loco", cores=ncores)
k.overall <- calc_kinship(genoprobs, cores=ncores)
t1 <- Sys.time()
message("Kinship calculated in ", round(as.numeric(difftime(t1,t0,units="secs")),1), " sec")

saveRDS(k, OUTPUT_KINSHIP_LOCO)
saveRDS(k.overall, OUTPUT_KINSHIP_OVL)
message("Done.")

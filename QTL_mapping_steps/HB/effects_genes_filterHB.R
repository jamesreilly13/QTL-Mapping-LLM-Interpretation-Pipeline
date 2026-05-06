#!/usr/bin/env Rscript
# qtl_effects_filter.R
# For a single phenotype:
#  - rebuild phenotype + covariates
#  - determine chr from scan1 + permutations
#  - run scan1blup
#  - query genes + variants in LOD interval
#  - run scan1snps
#  - plot CC allele effects + SNP association

suppressPackageStartupMessages({
  library(qtl2)
  library(tidyverse)
  library(qtl2convert)
})

## ----------------------------
## PHENOTYPE SELECTION
## ----------------------------

default_pheno <- "total.time.in.corner"
args <- commandArgs(trailingOnly = TRUE)
PHENO_COL <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  Sys.getenv("PHENO_COL", unset = default_pheno)
}

message("PHENO_COL set to: ", PHENO_COL)

## ----------------------------
## PHENOTYPE GROUP + SAFE NAME
## ----------------------------

safe_name <- function(x) gsub("[^A-Za-z0-9_\\-]", "_", x)

PHENO_GROUP <- "HB"
safe_pheno  <- safe_name(PHENO_COL)

pheno_dirname <- paste0(PHENO_GROUP, "_", safe_pheno)
message("Phenotype directory name:", pheno_dirname)

## ----------------------------
## PATHS
## ----------------------------

DATA_DIR    <- "/projects/csna/colby2025/QTLdata"
RESULTS_DIR <- "/projects/csna/colby2025/results"

pheno_dir <- file.path(RESULTS_DIR, pheno_dirname)
outdir_ph <- file.path(pheno_dir, "qtl_effects")

dir.create(outdir_ph, showWarnings = FALSE, recursive = TRUE)

message("Phenotype results directory: ", pheno_dir)
message("QTL effects output directory: ", outdir_ph)

genes_db    <- file.path(DATA_DIR, "mouse_genes.sqlite")
variants_db <- file.path(DATA_DIR, "cc_variants.sqlite")

## ----------------------------
## CORES
## ----------------------------

ncores_env <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA)
ncores <- if (!is.na(ncores_env) && nzchar(ncores_env)) {
  as.integer(ncores_env)
} else {
  parallel::detectCores(logical = TRUE)
}
if (is.na(ncores) || ncores < 1) ncores <- 1L
options(mc.cores = ncores)
message("Using ", ncores, " cores")

## ----------------------------
## HELPERS
## ----------------------------

norm_rank_transform <- function(x, c = 0) {
  x_noNA <- which(!is.na(x))
  N <- length(x_noNA)
  x[x_noNA] <- qnorm((rank(x[x_noNA], ties.method = "average") - c) /
                       (N - (2 * c) + 1))
  x
}

## ----------------------------
## LOAD DATA
## ----------------------------

message("Loading QC + covariates...")
load(file.path(DATA_DIR, "gm_13batches_newid_qc.RData"))
stopifnot(exists("gm_after_qc"))

message("Loading genotype probabilities...")
apr.revlearn <- readRDS(file.path(DATA_DIR, "allele_probs.rds"))

message("Loading kinship...")
k <- readRDS(file.path(DATA_DIR, "kinship_loco.rds"))

message("Reading phenotype CSV...")
pheno.do <- read.csv(
  file.path(DATA_DIR, "HB_residuals_07012020.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

message("Reading marker grid...")
grid <- read.delim(
  file.path(DATA_DIR, "marker_grid_0.02cM_plus.txt"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
grid.map.list <- map_df_to_list(
  grid,
  pos_column = "pos",
  chr_column = "chr",
  marker_column = "marker"
)

pmap <- gm_after_qc$pmap

## ----------------------------
## BUILD PHENOTYPE OBJECT
## ----------------------------

message("Merging covariates and phenotype data...")
covar <- gm_after_qc$covar

pheno.obj <- merge(
  covar,
  pheno.do,
  by.x = "name",
  by.y = "Mouse.ID",
  all = FALSE
)

if (!(PHENO_COL %in% colnames(pheno.obj))) {
  stop("Phenotype not found: ", PHENO_COL)
}

message("Normal-rank transforming phenotype...")
pheno_raw <- pheno.obj[[PHENO_COL]]
pheno1 <- norm_rank_transform(pheno_raw)
names(pheno1) <- pheno.obj$name

sex <- (pheno.obj$sex == "M") * 1
gen <- pheno.obj$ngen
addcovar <- model.matrix(~ sex + gen)[, -1, drop = FALSE]
rownames(addcovar) <- pheno.obj$name

## ----------------------------
## LOAD SCAN1 + PERMUTATIONS
## ----------------------------

scan1_file <- file.path(pheno_dir, paste0(safe_pheno, "_scan1_out.rds"))
perm_file  <- "/projects/csna/colby2025/results/permutation_total.time.in.corner_n1000.rds"

message("Loading scan1 output from: ", scan1_file)
if (!file.exists(scan1_file)) stop("Missing scan1 file: ", scan1_file)
out <- readRDS(scan1_file)

LOD_cutoff <- 6
message("Applying LOD cutoff: ", LOD_cutoff)

message(
  "scan1 dim BEFORE filter: ",
  paste(dim(out), collapse = " x ")
)

out_filt <- out
out_filt[out_filt < LOD_cutoff] <- NA

message(
  "scan1 dim AFTER filter:  ",
  paste(dim(out_filt), collapse = " x ")
)

message("Loading permutation results from: ", perm_file)
if (!file.exists(perm_file)) stop("Missing permutation file: ", perm_file)
operm <- readRDS(perm_file)

## ----------------------------
## DETERMINE CHROMOSOME
## ----------------------------

message("Computing significance threshold...")
thr <- summary(operm, alpha = 0.10) # significance threshold

message("Finding peaks...")
peaks <- find_peaks(
  scan1_output = out_filt,
  map = pmap,
  threshold = thr,
  prob = 0.95,
  expand2markers = FALSE
)

print(peaks)

if (nrow(peaks) == 0) {
  stop("No significant peaks detected for ", PHENO_COL)
}

chr <- as.character(peaks$chr[1])
message(">>> Chromosome of interest: ", chr)

## ----------------------------
## scan1blup
## ----------------------------

message("Running scan1blup on chr ", chr)

coef <- scan1blup(
  genoprobs = apr.revlearn[, chr],
  pheno     = pheno1,
  kinship   = k[[chr]],
  addcovar  = addcovar,
  cores     = ncores
)

saveRDS(coef, file = file.path(outdir_ph, paste0("coef_filter_chr", chr, ".rds")))

## ----------------------------
## GENE + VARIANT QUERY
## ----------------------------

message("Creating gene query function (MGI)...")
query_genes <- create_gene_query_func(
  dbfile = genes_db,
  filter = "source='MGI'"
)

message("Creating variant query function...")
query_variants <- create_variant_query_func(variants_db)

message("Computing LOD support interval...")
lod_int_chr <- lod_int(out_filt, map = pmap, chr = chr)

lowest  <- lod_int_chr[1]
peakpos <- lod_int_chr[2]
highest <- lod_int_chr[3]

message(sprintf(
  "LOD interval chr %s: %.3f – %.3f (peak %.3f)",
  chr, lowest, highest, peakpos
))

genes <- query_genes(chr, lowest, highest)
variants <- query_variants(chr, lowest, highest)

saveRDS(genes,    file = file.path(outdir_ph, paste0("genes_filter_chr", chr, ".rds")))
saveRDS(variants, file = file.path(outdir_ph, paste0("variants_filter_chr", chr, ".rds")))

## ----------------------------
## SNP ASSOCIATION (scan1snps)
## ----------------------------

message("Running scan1snps on chr ", chr)

out_snps <- scan1snps(
  genoprobs = apr.revlearn[, chr],
  map       = pmap,
  pheno     = pheno1,
  kinship   = k[[chr]],
  addcovar  = addcovar,
  query_func = query_variants,
  chr       = chr,
  start     = lowest,
  end       = highest,
  keep_all_snps = TRUE,
  cores = ncores
)

saveRDS(out_snps, file = file.path(outdir_ph, paste0("scan1snps_filter_chr", chr, ".rds")))

## ----------------------------
## PLOTS
## ----------------------------

# SNP association plot
pdf(
  file = file.path(outdir_ph, paste0("snpasso_filter_chr", chr, ".pdf")),
  width = 10,
  height = 6
)

plot_snpasso(
  out_snps$lod,
  out_snps$snpinfo,
  drop_hilit = 1.5,
  genes = genes
)

dev.off()

# CC effects plot
pdf(
  file = file.path(outdir_ph, paste0("coef_filter_chr", chr, ".pdf")),
  width = 10,
  height = 6
)

plot_coefCC(
  x = coef,
  map = pmap,
  scan1_output = out,
  main = paste0("QTL Effects – Chr ", chr, "\n", PHENO_COL),
  legend = "bottomleft"
)

dev.off()

message("Done with phenotype: ", PHENO_COL)

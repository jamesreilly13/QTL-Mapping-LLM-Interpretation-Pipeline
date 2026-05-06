#!/usr/bin/env Rscript
# run_qtl2_onepheno.R
# Single-phenotype scan1 test with verbose diagnostics (no permutations).

# --- Allow PHENO_COL to come from env var or command-line arg (fallback to default)
default_pheno <- "total.time.in.corner"
PHENO_COL_arg <- if (length(commandArgs(trailingOnly = TRUE)) >= 1) commandArgs(trailingOnly = TRUE)[1] else NA
PHENO_COL_env <- Sys.getenv("PHENO_COL", unset = NA)

if (!is.na(PHENO_COL_arg) && nzchar(PHENO_COL_arg)) {
  PHENO_COL <- PHENO_COL_arg
} else if (!is.na(PHENO_COL_env) && nzchar(PHENO_COL_env)) {
  PHENO_COL <- PHENO_COL_env
} else {
  PHENO_COL <- default_pheno
}
message("PHENO_COL set to: ", PHENO_COL)

DATA_DIR <- "/projects/csna/colby2025/QTLdata"
OUT_DIR  <- "/projects/csna/colby2025/results"

# basic setup
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
message("Output will be written to: ", OUT_DIR)

pkgs <- c("qtl2", "tidyverse", "qtl2convert", "parallel")
missing_pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "),
       "\nPlease ensure these are installed.")
}
library(qtl2); library(tidyverse); library(qtl2convert)

# cores
ncores_env <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA)
if (!is.na(ncores_env) && nzchar(ncores_env)) {
  ncores <- suppressWarnings(as.integer(ncores_env))
} else {
  ncores <- parallel::detectCores(logical = TRUE)
}
if (is.na(ncores) || ncores < 1) ncores <- 1L
message("Using ", ncores, " cores")
options(mc.cores = ncores)

# safe filename helper
safe_name <- function(x) gsub("[^A-Za-z0-9_\\-]", "_", x)

# normal rank transformation function
norm_rank_transform <- function(x, c = 0) {
  stopifnot(is.numeric(x) & is.vector(x))
  x_noNA = which(!is.na(x))
  N = length(x_noNA)
  x[x_noNA] = qnorm((rank(x[x_noNA], ties.method = "average") - c) / (N - (2 * c) + 1))
  return(x)
}

# -- paths
rdata_gm    <- file.path(DATA_DIR, "gm_13batches_newid_qc.RData")
rds_apr     <- file.path(DATA_DIR, "allele_probs.rds")
rdata_apr   <- file.path(DATA_DIR, "LD_residuals_07012020.csv")
grid_file   <- file.path(DATA_DIR, "marker_grid_0.02cM_plus.txt")
ovr_kin     <- file.path(DATA_DIR, "kinship_overall.rds")
loco_kin    <- file.path(DATA_DIR, "kinship_loco.rds")
sqlite_vari <- file.path(DATA_DIR, "cc_variants.sqlite")
genes_db    <- file.path(DATA_DIR, "mouse_genes.sqlite")

files_needed <- c(rdata_gm, rds_apr, rdata_apr, grid_file, ovr_kin, loco_kin)
for (f in files_needed) {
  if (!file.exists(f)) stop("Required file not found: ", f)
}
message("All required files present.")

# load / read
message("Loading RData: ", rdata_gm)
load(rdata_gm)
if (!exists("gm_after_qc")) stop("Expected object 'gm_after_qc' not found in ", rdata_gm)
message("Loading genotype probs from: ", rds_apr)
genoprobs <- readRDS(rds_apr)
message("Loading kinships...")
k <- readRDS(loco_kin)
k.overall <- readRDS(ovr_kin)

message("Reading phenotype CSV: ", rdata_apr)
pheno.do <- read.csv(rdata_apr, stringsAsFactors = FALSE, check.names = FALSE)

# map / pmap
message("Reading grid file: ", grid_file)
grid <- read.delim(grid_file, stringsAsFactors = FALSE, check.names = FALSE)
grid.map.list <- map_df_to_list(grid, pos_column = "pos", chr_column = "chr", marker_column = "marker")
pmap <- gm_after_qc$pmap

# Merge covariates & phenotypes
message("Merging covariates (gm_after_qc$covar) with phenotype CSV by covar$name <-> pheno.Mouse.ID")
covar <- gm_after_qc$covar
if (!("name" %in% colnames(covar))) stop("covar missing 'name' column")
if (!("Mouse.ID" %in% colnames(pheno.do))) stop("pheno CSV missing 'Mouse.ID' column")
pheno.obj <- merge.data.frame(covar, pheno.do, by.x = "name", by.y = "Mouse.ID", all.x = FALSE, all.y = FALSE)
message("Merged data has ", nrow(pheno.obj), " rows and ", ncol(pheno.obj), " columns")

# covariates
sex <- (pheno.obj$sex == "M") * 1; names(sex) <- pheno.obj$name
gen <- pheno.obj$ngen; names(gen) <- pheno.obj$name
addcovar <- model.matrix(~sex + gen, data = pheno.obj)[, -1, drop = FALSE]
row.names(addcovar) <- pheno.obj$name

# prepare phenotype vector (plain named numeric)
if (!PHENO_COL %in% colnames(pheno.obj)) stop("Phenotype column not found: ", PHENO_COL)
pheno1 <- suppressWarnings(as.numeric(pheno.obj[[PHENO_COL]]))
names(pheno1) <- pheno.obj$name

# helper to get genotype sample IDs
get_geno_ids <- function(gp) {
  if (is.list(gp)) {
    first <- gp[[1]]
    return(rownames(first))
  } else if (is.array(gp)) {
    return(dimnames(gp)[[1]])
  } else if (is.matrix(gp) || is.data.frame(gp)) {
    return(rownames(gp))
  } else stop("Unrecognized genoprobs structure")
}
geno_ids <- get_geno_ids(genoprobs)
message("Genoprobs sample count: ", length(geno_ids))
message("Phenotype vector sample count: ", length(pheno1))

# sample overlap diagnostics
common_ids <- intersect(names(pheno1), geno_ids)
message("Samples overlapping between phenotype & genotype: ", length(common_ids), " / ", length(geno_ids))
if (length(common_ids) == 0) {
  message("Example phenotype names (first 10): ", paste(head(names(pheno1), 10), collapse = ", "))
  message("Example genotype IDs (first 10): ", paste(head(geno_ids, 10), collapse = ", "))
  stop("No overlapping sample IDs between phenotype and genotype. Fix sample ids before scanning.")
}
# print few examples of missing/extra ids
missing_from_pheno <- setdiff(geno_ids, names(pheno1))
missing_from_geno  <- setdiff(names(pheno1), geno_ids)
if (length(missing_from_pheno) > 0) {
  message("Number of genotype IDs missing phenotype: ", length(missing_from_pheno))
  message("Examples (up to 10) genotype IDs missing phenotype: ", paste(head(missing_from_pheno, 10), collapse = ", "))
}
if (length(missing_from_geno) > 0) {
  message("Number of phenotype IDs missing genotype: ", length(missing_from_geno))
  message("Examples (up to 10) phenotype IDs missing genotype: ", paste(head(missing_from_geno, 10), collapse = ", "))
}

# subset & reorder phenotype to genotype ordering (keeps only genotype-present samples)
geno_order_present <- geno_ids[geno_ids %in% names(pheno1)]
pheno1_sub <- pheno1[geno_order_present]
message("Using ", length(pheno1_sub), " samples (phenotype ordered to genotype sample order)")

# Utilize rank-based normal transformation
message("Applying rank-based inverse normal transform to phenotype: ", PHENO_COL)
pheno1_sub <- norm_rank_transform(pheno1_sub, c = 0)
message("Phenotype transform complete.")

# Free memory from large temporary objects that are no longer needed,
# then run the garbage collector before the heavy scan1 call.
# Keep only objects needed for scan1: genoprobs, k, addcovar, pheno1_sub, and plot map/pmap if required.
# Remove large originals if present.
message("Removing things for memory efficiency...")
to_rm <- c("pheno.do", "covar", "pheno.obj", "grid", "gm_after_qc")
for (nm in to_rm) {
  if (exists(nm, inherits = FALSE)) {
    rm(list = nm)
  }
}
# call gc to release memory back to the OS if possible
invisible(gc())
message("Removals completed.")

# run scan1 (single phenotype) with diagnostics
message("Running scan1 for phenotype: ", PHENO_COL)
t0 <- Sys.time()
scan_out <- tryCatch({
  scan1(genoprobs = genoprobs, pheno = pheno1_sub, kinship = k, addcovar = addcovar, cores = ncores)
}, error = function(e) {
  stop("scan1 failed: ", conditionMessage(e))
})
t1 <- Sys.time()
message("scan1 completed in ", round(as.numeric(difftime(t1, t0, units = "secs")), 1), " sec")

# print info about scan markers & map markers
scan_markers <- rownames(scan_out)
message("Number of markers in scan1 output: ", length(scan_markers))
# markers from pmap and grid (show first 10)
pmap_markers <- unlist(lapply(pmap, names))
grid_markers <- unlist(lapply(grid.map.list, names))
message("Number of markers in pmap: ", length(pmap_markers), "; in grid map: ", length(grid_markers))
# intersection counts
n_common_pmap  <- length(intersect(scan_markers, pmap_markers))
n_common_grid  <- length(intersect(scan_markers, grid_markers))
message("Markers in common: scan vs pmap = ", n_common_pmap, "; scan vs grid = ", n_common_grid)
if (n_common_pmap == 0 && n_common_grid == 0) {
  message("Example scan markers (first 10): ", paste(head(scan_markers, 10), collapse = ", "))
  message("Example pmap markers (first 10): ", paste(head(pmap_markers, 10), collapse = ", "))
  message("Example grid markers (first 10): ", paste(head(grid_markers, 10), collapse = ", "))
  stop("No marker name overlap between scan1 output and either pmap or grid map. Use matching genoprobs/map.")
}

# choose map for plotting (prefer pmap if it overlaps)
plot_map <- if (n_common_pmap > 0) {
  message("Using pmap for plotting (has overlap).")
  pmap
} else {
  message("Using grid.map.list for plotting (pmap had zero overlap).")
  grid.map.list
}

# plotting and saving
outdir_ph <- file.path(OUT_DIR, paste0("HB_", safe_name(PHENO_COL)))
dir.create(outdir_ph, recursive = TRUE, showWarnings = FALSE)
plot_file <- file.path(outdir_ph, paste0(safe_name(PHENO_COL), "_scan1.png"))

message("Plotting scan1 to: ", plot_file)
png(plot_file, width = 1400, height = 900)
tryCatch({
  plot_scan1(scan_out, map = plot_map, main = paste0("scan1: ", PHENO_COL))
}, error = function(e) {
  dev.off()
  stop("plot_scan1 failed: ", conditionMessage(e))
})
dev.off()

# save scan object
saveRDS(scan_out, file = file.path(outdir_ph, paste0(safe_name(PHENO_COL), "_scan1_out.rds")))
message("Saved scan1 object to: ", file.path(outdir_ph, paste0(safe_name(PHENO_COL), "_scan1_out.rds")))

message("Single-phenotype pipeline complete.")

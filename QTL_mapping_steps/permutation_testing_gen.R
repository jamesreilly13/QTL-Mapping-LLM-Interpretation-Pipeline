#!/usr/bin/env Rscript
# permutation_testing_gen.R
#PERMUTATION TESTING SCRIPT

library(qtl2)

pheno <- "total.time.in.corner"

DATA_DIR <- "/projects/csna/colby2025/QTLdata"
OUT_DIR  <- "/projects/csna/colby2025/results"

run_permutations_for_pheno <- function(pheno_vec, genoprobs, kinship, addcovar, n_perm = 100, cores = 1) {
  message("  [perm] Starting permutation test (n_perm=", n_perm, ") ...")
  t0 <- Sys.time()
  operm <- scan1perm(genoprobs = genoprobs, pheno = pheno_vec, kinship = kinship, addcovar = addcovar,
                     n_perm = n_perm, cores = cores, verbose = TRUE)
  t1 <- Sys.time()
  message("  [perm] Completed in ", round(as.numeric(difftime(t1, t0, units = "secs")), 1), " sec")
  return(operm)
}

# normal rank transformation function
norm_rank_transform <- function(x, c = 0) {
  stopifnot(is.numeric(x) & is.vector(x))
  x_noNA = which(!is.na(x))
  N = length(x_noNA)
  x[x_noNA] = qnorm((rank(x[x_noNA], ties.method = "average") - c) / (N - (2 * c) + 1))
  return(x)
}

ncores_env <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA)
if (!is.na(ncores_env) && nzchar(ncores_env)) {
  ncores <- suppressWarnings(as.integer(ncores_env))
} else {
  ncores <- parallel::detectCores(logical = TRUE)
}
if (is.na(ncores) || ncores < 1) ncores <- 1L
message("Using ", ncores, " cores")
options(mc.cores = ncores)

rdata_gm    <- file.path(DATA_DIR, "gm_13batches_newid_qc.RData")
rds_apr     <- file.path(DATA_DIR, "allele_probs.rds")
rdata_apr   <- file.path(DATA_DIR, "OFA_residuals_07012020.csv")
grid_file   <- file.path(DATA_DIR, "marker_grid_0.02cM_plus.txt")
ovr_kin     <- file.path(DATA_DIR, "kinship_overall.rds")
loco_kin    <- file.path(DATA_DIR, "kinship_loco.rds")
sqlite_vari <- file.path(DATA_DIR, "cc_variants.sqlite")
genes_db    <- file.path(DATA_DIR, "mouse_genes.sqlite")

message("Loading RData: ", rdata_gm)
load(rdata_gm)
if (!exists("gm_after_qc")) stop("Expected object 'gm_after_qc' not found in ", rdata_gm)

message("Loading genotype probs from: ", rds_apr)
genoprobs <- readRDS(rds_apr)

message("Loading kinships...")
k <- readRDS(loco_kin)

message("Reading phenotype CSV: ", rdata_apr)
pheno.do <- read.csv(rdata_apr, stringsAsFactors = FALSE, check.names = FALSE)

message("Merging covariates (gm_after_qc$covar) with phenotype CSV by covar$name <-> pheno.Mouse.ID")
covar <- gm_after_qc$covar
if (!("name" %in% colnames(covar))) stop("covar missing 'name' column")
if (!("Mouse.ID" %in% colnames(pheno.do))) stop("pheno CSV missing 'Mouse.ID' column")
pheno.obj <- merge.data.frame(covar, pheno.do, by.x = "name", by.y = "Mouse.ID", all.x = FALSE, all.y = FALSE)
message("Merged data has ", nrow(pheno.obj), " rows and ", ncol(pheno.obj), " columns")

message("Creating addcovar object...")
sex <- (pheno.obj$sex == "M") * 1; names(sex) <- pheno.obj$name
gen <- pheno.obj$ngen; names(gen) <- pheno.obj$name
addcovar <- model.matrix(~sex + gen, data = pheno.obj)[, -1, drop = FALSE]
row.names(addcovar) <- pheno.obj$name

if (! (pheno %in% colnames(pheno.obj)) ) {
  stop("Phenotype '", pheno, "' not found in pheno.obj columns. Available columns: ",
       paste(colnames(pheno.obj), collapse = ", "))
}

pheno_raw <- pheno.obj[[pheno]]
pheno_norm <- norm_rank_transform(pheno_raw) 

names(pheno_norm) <- pheno.obj$name

n_perm_to_run <- 1000
operm <- run_permutations_for_pheno(pheno_vec = pheno_norm,
                                    genoprobs = genoprobs,
                                    kinship = k,
                                    addcovar = addcovar,
                                    n_perm = n_perm_to_run,
                                    cores = ncores)

saveRDS(operm, file = file.path(DATA_DIR, paste0("permutation_", pheno, "_n", n_perm_to_run, ".rds")))
message("Permutation results saved.")

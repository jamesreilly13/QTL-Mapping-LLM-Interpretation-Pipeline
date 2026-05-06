This subdirectory contains the process of QTL mapping given our input data. Input data required to run these analyses are as follows:

1. apr.revlearn.69K.02182020 -> genotype probabilities
2. cc_variants.sqlite -> variant map utilized to create variant query function
3. histoUCF_Grcm38_map -> gene name to ENSEMBL ID map
4. marker_grid_0.02cM_plus -> marker map for DO mice
5. mouse_genes.sqlite -> mouse gene annotations
6. gm_DO2816_qc -> cross2 QTL data object

As well as the residual data files for the Holeboard (HB), Light dark (LD), Novelty Place Preference (NPP), and Open Field Activity (OFA) phenotype 
groups (i.e. HB_residuals_07012020.csv).

The analysis performed by these scripts must be performed in the correct order as they build on eachother:

1. Convert genotype probabilities to allele probabilities (geno_to_allele.R)
2. Calculate kinship matrix (kinship_calc.R)
3. Permutation testing to get a threshold (permutation_testing_gen.R)
4. scan1 genome scan (for each group (i.e. NPP, OFA, etc.) it is the "scan1_*.R" file))
5. Determine the coefficient effects, SNPs, variants, genes within QTL region with "effects_genes_filter*.R" files for each group
6. Extract LOD intervals (extract_lod_intervals.sh)
7. Sync HPC outputs from the above analyses to local computer
8. Quality control to remove coefficient effects on invalid marker regions (Coef_LOD_trimmer.R)
   and remove scan1 outputs less than 6 LOD (filter_scan1_out.R)
9. Merge all phenotypes by group (merge folder)
10. Create qtl_data_trimmed object in R (build_qtl_data_trimmed.R)
11. Perform gene enrichment analysis from above object (gene_enrichment_v2_fromTrimmedObj.R)
12. Rename ENSEMBL IDs within qtl_trim object (mapping_ensembles.R)
13. Export object to Python for LLM process!

Important note regarding this subdirectory:

The following R-scripts included in this directory were ran on JAX's HPC and customized to specific phenotypes,
so they therefore will require a SUBMISSION script to be ran and a good amount of finagling with file names for successful running:

gene_to_allele.R, kinship_calc.R, permutation_testing_gen.R, all versions of effects_genes_filter*.R, and all versions of scan1_*.R.

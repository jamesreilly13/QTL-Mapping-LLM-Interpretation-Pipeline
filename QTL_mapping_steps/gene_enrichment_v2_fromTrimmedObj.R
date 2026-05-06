# James Reilly

# ----------------------------
# Re-run GO enrichment on cleaned gene sets
# ----------------------------

library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(stringr)
library(org.Mm.eg.db)
library(tibble)

# ----------------------------
# Helper: extract gene IDs
# ----------------------------
extract_gene_ids <- function(df) {
  tibble(
    ENSEMBL = str_extract(df$Dbxref, "ENSMUSG[0-9]+"),
    ENTREZ  = str_extract(df$Dbxref, "NCBI_Gene:[0-9]+") %>%
      str_remove("NCBI_Gene:")
  )
}

# ----------------------------
# Run GO enrichment for one phenotype
# ----------------------------
run_go_for_pheno <- function(gene_df, phenotype_name) {
  
  ids <- extract_gene_ids(gene_df)
  
  # Map ENSEMBL → ENTREZ
  mapped <- bitr(
    ids$ENSEMBL[!is.na(ids$ENSEMBL)],
    fromType = "ENSEMBL",
    toType   = "ENTREZID",
    OrgDb    = org.Mm.eg.db
  )
  
  gene_list <- unique(na.omit(c(
    mapped$ENTREZID,
    ids$ENTREZ
  )))
  
  message("  ", phenotype_name, ": ", length(gene_list), " mapped genes")
  
  if (length(gene_list) == 0) return(NULL)
  
  run_one <- function(ont) {
    res <- enrichGO(
      gene          = gene_list,
      OrgDb         = org.Mm.eg.db,
      keyType       = "ENTREZID",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.2,
      readable      = TRUE
    )
    
    if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
    
    df <- as.data.frame(res)
    df$ontology  <- ont
    df$phenotype <- phenotype_name
    
    return(df)
  }
  
  bind_rows(
    run_one("BP"),
    run_one("MF"),
    run_one("CC")
  )
}

# ----------------------------
# MAIN: loop over groups
# ----------------------------

for (g in names(qtl_data_trimmed)) {
  
  message("Re-running GO enrichment for group: ", g)
  
  genes_df <- qtl_data_trimmed[[g]]$genes
  
  if (!"phenotype" %in% colnames(genes_df)) {
    stop("Missing 'phenotype' column in genes table for group: ", g)
  }
  
  phenos <- unique(genes_df$phenotype)
  
  go_results_all <- list()
  
  for (p in phenos) {
    
    message(" Processing phenotype: ", p)
    
    sub_genes <- genes_df %>%
      filter(phenotype == p)
    
    res <- run_go_for_pheno(sub_genes, p)
    
    if (!is.null(res)) {
      go_results_all[[p]] <- res
    }
  }
  
  combined_go <- bind_rows(go_results_all)
  
  message(" Total GO rows for ", g, ": ", nrow(combined_go))
  
  # ----------------------------
  # Replace GO slot
  # ----------------------------
  qtl_data_trimmed[[g]]$go <- combined_go
}

message("✅ GO enrichment successfully updated for all groups.")
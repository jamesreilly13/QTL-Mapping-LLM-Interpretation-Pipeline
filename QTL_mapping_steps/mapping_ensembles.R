library(dplyr)
library(stringr)
library(tidyverse)

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



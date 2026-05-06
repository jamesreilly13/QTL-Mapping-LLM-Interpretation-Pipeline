# James Reilly
# scan1snps_merger.R

library(dplyr)
library(purrr)
library(stringr)
library(tidyr)

base_dir <- "/Users/reillj/Documents/QTL/results"

# Get all group directories
group_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
group_dirs <- group_dirs[grepl("^HPC_filter_output_", basename(group_dirs))]

for (group_path in group_dirs) {
  
  group <- str_remove(basename(group_path), "^HPC_filter_output_")
  message("Processing scan1snps for group: ", group)
  
  phenotype_dirs <- list.dirs(group_path, recursive = FALSE, full.names = TRUE)
  phenotype_dirs <- phenotype_dirs[!grepl("_filter_combined$", basename(phenotype_dirs))]
  
  # Combine all scan1snps into one dataframe
  scan1snps_combined <- map_dfr(phenotype_dirs, function(pheno_dir) {
    
    phenotype <- basename(pheno_dir)
    
    files <- list.files(
      pheno_dir,
      pattern = "scan1snps_filter_chr.*\\.rds$",
      full.names = TRUE
    )
    
    if (length(files) == 0) return(NULL)
    
    map_dfr(files, function(f) {
      
      snp_obj <- tryCatch(readRDS(f), error = function(e) NULL)
      if (is.null(snp_obj)) return(NULL)
      
      # Convert snpinfo to dataframe
      snp_df <- snp_obj$snpinfo %>%
        as.data.frame() 
      
      # Convert LOD matrix to long dataframe
      lod_df <- snp_obj$lod %>%
        as.data.frame() %>%
        mutate(snp = rownames(.)) %>%
        pivot_longer(
          cols = -snp,
          names_to = "phenotype_column",
          values_to = "lod"
        )
      
      # Join LODs to snpinfo
      df <- snp_df %>%
        left_join(lod_df, by = "snp") %>%
        mutate(
          group = group,
          phenotype = phenotype,
          source_file = basename(f)
        )
      
      return(df)
    })
  })
  
  if (nrow(scan1snps_combined) == 0) {
    message("  No scan1snps files found for ", group)
    next
  }
  
  # Save to combined directory
  combined_dir <- file.path(group_path, paste0(group, "_filter_combined"))
  if (!dir.exists(combined_dir)) dir.create(combined_dir)
  
  out_file <- file.path(
    combined_dir,
    paste0(group, "_combined_scan1snps.rds")
  )
  
  saveRDS(scan1snps_combined, out_file)
  message("  Saved: ", out_file)
}

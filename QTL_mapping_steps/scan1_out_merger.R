# James Reilly
# scan1_out_merger.R

library(dplyr)
library(purrr)
library(stringr)
library(tidyr)

base_dir <- "/Users/reillj/Documents/QTL/results"

# Find group directories
group_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
group_dirs <- group_dirs[grepl("^HPC_output_", basename(group_dirs))]

group_dirs <- "/Users/reillj/Documents/QTL/results/HPC_output_NPP"

for (group_path in group_dirs) {
  
  group <- str_remove(basename(group_path), "^HPC_output_")
  message("Processing scan1_out for group: ", group)
  
  phenotype_dirs <- list.dirs(
    group_path,
    recursive = FALSE,
    full.names = TRUE
  )
  
  # Skip combined folders
  phenotype_dirs <- phenotype_dirs[
    !grepl("_combined$", basename(phenotype_dirs))
  ]
  
  scan1_combined <- map_dfr(phenotype_dirs, function(pheno_dir) {
    
    phenotype <- basename(pheno_dir)
    
    files <- list.files(
      pheno_dir,
      pattern = "_scan1_out\\.rds$",
      full.names = TRUE
    )
    
    if (length(files) == 0) return(NULL)
    
    map_dfr(files, function(f) {
      
      scan1_obj <- tryCatch(readRDS(f), error = function(e) NULL)
      if (is.null(scan1_obj)) return(NULL)
      
      # Convert to dataframe
      df <- as.data.frame(scan1_obj) %>%
        mutate(marker = rownames(scan1_obj)) %>%
        relocate(marker)
      
      # Handle multiple phenotype columns safely
      df_long <- df %>%
        pivot_longer(
          cols = -marker,
          names_to = "pheno_column",
          values_to = "lod"
        ) %>%
        mutate(
          group = group,
          phenotype = phenotype,
          source_file = basename(f)
        )
      
      return(df_long)
    })
  })
  
  if (nrow(scan1_combined) == 0) {
    message("  No scan1_out files found for ", group)
    next
  }
  
  # Save output
  combined_dir <- file.path(group_path, paste0(group, "_combined"))
  if (!dir.exists(combined_dir)) dir.create(combined_dir)
  
  out_file <- file.path(
    combined_dir,
    paste0(group, "_combined_scan1_out.rds")
  )
  
  saveRDS(scan1_combined, out_file)
  message("  Saved: ", out_file)
}

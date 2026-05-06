# James Reilly
# coef_merger.R

library(dplyr)
library(purrr)
library(stringr)

base_dir <- "/Users/reillj/Documents/QTL/results"

# Get group directories
group_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
group_dirs <- group_dirs[grepl("^HPC_filter_output_", basename(group_dirs))]

for (group_path in group_dirs) {
  
  group <- str_remove(basename(group_path), "^HPC_filter_output_")
  message("Processing coef for group: ", group)
  
  phenotype_dirs <- list.dirs(
    group_path,
    recursive = FALSE,
    full.names = TRUE
  )
  
  phenotype_dirs <- phenotype_dirs[
    !grepl("_filter_combined$", basename(phenotype_dirs))
  ]
  
  # Combine everything into one dataframe
  coef_combined <- map_dfr(phenotype_dirs, function(pheno_dir) {
    
    phenotype <- basename(pheno_dir)
    
    files <- list.files(
      pheno_dir,
      pattern = "coef_filter_chr.*\\.rds$",
      full.names = TRUE
    )
    
    if (length(files) == 0) return(NULL)
    
    map_dfr(files, function(f) {
      
      coef_obj <- tryCatch(readRDS(f), error = function(e) NULL)
      if (is.null(coef_obj)) return(NULL)
      
      # Convert matrix to dataframe
      df <- as.data.frame(coef_obj) %>%
        mutate(sample = rownames(coef_obj)) %>%
        relocate(sample) %>%
        mutate(
          group = group,
          phenotype = phenotype,
          source_file = basename(f)
        )
      
      return(df)
    })
  })
  
  if (nrow(coef_combined) == 0) {
    message("  No coef files found for ", group)
    next
  }
  
  # Save to combined directory
  combined_dir <- file.path(group_path, paste0(group, "_filter_combined"))
  if (!dir.exists(combined_dir)) dir.create(combined_dir)
  
  out_file <- file.path(
    combined_dir,
    paste0(group, "_combined_coef.rds")
  )
  
  saveRDS(coef_combined, out_file)
  message("  Saved: ", out_file)
}

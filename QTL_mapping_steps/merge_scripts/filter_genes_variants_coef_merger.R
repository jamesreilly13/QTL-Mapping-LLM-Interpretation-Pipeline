# James Reilly
# genes_variants_merger.R

library(dplyr)
library(purrr)
library(stringr)

base_dir <- "/Users/reillj/Documents/QTL/results"

# All group folders
group_dirs <- list.dirs(
  base_dir,
  recursive = FALSE,
  full.names = TRUE
)

# Keep only directories named like "HPC_filter_output_*"
group_dirs <- group_dirs[
  grepl("^HPC_filter_output_", basename(group_dirs))
]

# Patterns of files to combine
file_patterns <- c(
  genes       = "genes_filter_chr.*\\.rds$",
  variants    = "variants_filter_chr.*\\.rds$"
)

# Helper function to safely read RDS
safe_read_rds <- function(file) {
  tryCatch(
    readRDS(file),
    error = function(e) NULL
  )
}

for (group_path in group_dirs) {
  
  group_name <- str_remove(basename(group_path), "^HPC_filter_output_")
  message("Processing group: ", group_name)
  
  phenotype_dirs <- list.dirs(group_path, recursive = FALSE, full.names = TRUE)
  
  phenotype_dirs <- phenotype_dirs[
    !grepl("_combined$", basename(phenotype_dirs))
  ]
  
  for (type_name in names(file_patterns)) {
    
    message("  Combining: ", type_name)
    
    combined <- map_dfr(phenotype_dirs, function(pheno_dir) {
      
      phenotype_name <- basename(pheno_dir)
      
      files <- list.files(
        pheno_dir,
        pattern = file_patterns[[type_name]],
        full.names = TRUE
      )
      
      if (length(files) == 0) return(NULL)
      
      map_dfr(files, function(f) {
        obj <- safe_read_rds(f)
        if (is.null(obj)) return(NULL)
        
        obj %>%
          mutate(
            group = group_name,
            phenotype = phenotype_name,
            source_file = basename(f)
          )
      })
    })
    
    if (nrow(combined) == 0) {
      message("    No files found for ", type_name)
      next
    }
    
    # Save to combined directory
    combined_dir <- file.path(group_path, paste0(group_name, "_filter_combined"))
    
    out_file <- file.path(
      combined_dir,
      paste0(group_name, "_combined_", type_name, ".rds")
    )
    
    saveRDS(combined, out_file)
    message("    Saved: ", out_file)
  }
}

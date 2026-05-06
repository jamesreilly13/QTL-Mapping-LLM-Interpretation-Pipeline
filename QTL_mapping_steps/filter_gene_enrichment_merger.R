# James Reilly
# gene_enrichment_merger.R

library(clusterProfiler)
library(dplyr)
library(purrr)
library(stringr)

base_dir <- "/Users/reillj/Documents/QTL/results"

group_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
group_dirs <- group_dirs[grepl("^HPC_filter_output_", basename(group_dirs))]

go_patterns <- c(
  MF = "_GO_MF\\.rds$",
  BP = "_GO_BP\\.rds$",
  CC = "_GO_CC\\.rds$"
)

for (group_path in group_dirs) {
  
  group <- str_remove(basename(group_path), "^HPC_filter_output_")
  message("Processing GO for group: ", group)
  
  phenotype_dirs <- list.dirs(
    group_path,
    recursive = FALSE,
    full.names = TRUE
  )
  
  phenotype_dirs <- phenotype_dirs[
    !grepl("_filter_combined$", basename(phenotype_dirs))
  ]
  
  all_go <- list()
  
  for (pheno in phenotype_dirs) {
    
    phenotype <- basename(pheno)
    
    for (ont in names(go_patterns)) {
      
      files <- list.files(
        pheno,
        pattern = go_patterns[[ont]],
        full.names = TRUE
      )
      
      if (length(files) == 0) next
      
      for (f in files) {
        
        er <- tryCatch(readRDS(f), error = function(e) NULL)
        
        if (is.null(er)) {
          message("    Skipping NULL GO object: ", f)
          next
        }
        
        if (!isS4(er) || !"result" %in% slotNames(er)) {
          message("    Skipping invalid GO object: ", f)
          next
        }
        
        if (nrow(er@result) == 0) {
          message("    Empty GO result: ", f)
          next
        }
        
        df <- er@result %>%
          mutate(
            group = group,
            phenotype = phenotype,
            ontology = ont,
            source_file = basename(f)
          )
        
        all_go[[length(all_go) + 1]] <- df
      }
    }
  }
  
  if (length(all_go) == 0) {
    message("  No GO files found for ", group)
    next
  }
  
  go_combined <- bind_rows(all_go)
  
  combined_dir <- file.path(group_path, paste0(group, "_filter_combined"))
  if (!dir.exists(combined_dir)) dir.create(combined_dir)
  
  out_file <- file.path(
    combined_dir,
    paste0(group, "_combined_GO.rds")
  )
  
  saveRDS(go_combined, out_file)
  message("  Saved: ", out_file)
}

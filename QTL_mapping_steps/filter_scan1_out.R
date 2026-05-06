# James Reilly - takes all scan1_outs and gets rid of LODs < 6

library()

OFA_scan1 <- readRDS("/Users/reillj/Documents/QTL/results/HPC_filter_output_OFA/OFA_filter_combined/OFA_combined_scan1_out.rds")
HB_scan1 <- readRDS("/Users/reillj/Documents/QTL/results/HPC_filter_output_HB/HB_filter_combined/HB_combined_scan1_out.rds")
LD_scan1 <- readRDS("/Users/reillj/Documents/QTL/results/HPC_filter_output_LD/LD_filter_combined/LD_combined_scan1_out.rds")
NPP_scan1 <- readRDS("/Users/reillj/Documents/QTL/results/HPC_output_NPP/NPP_combined/NPP_combined_scan1_out.rds")

filtered_OFA <- OFA_scan1[OFA_scan1$lod >= 6, ]
filtered_HB <- HB_scan1[HB_scan1$lod >= 6, ]
filtered_LD <- LD_scan1[LD_scan1$lod >= 6, ]
filtered_NPP <- NPP_scan1[NPP_scan1$lod >= 6, ]

saveRDS(filtered_OFA, "/Users/reillj/Documents/QTL/results/HPC_filter_output_OFA/OFA_filter_combined/OFA_filter_combined_scan1_out.rds")
saveRDS(filtered_HB, "/Users/reillj/Documents/QTL/results/HPC_filter_output_HB/HB_filter_combined/HB_filter_combined_scan1_out.rds")
saveRDS(filtered_OFA, "/Users/reillj/Documents/QTL/results/HPC_filter_output_LD/LD_filter_combined/LD_filter_combined_scan1_out.rds")
saveRDS(filtered_NPP, "/Users/reillj/Documents/QTL/results/HPC_filter_output_NPP/NPP_filter_combined/NPP_filter_combined_scan1_out.rds")

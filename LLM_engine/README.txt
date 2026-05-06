This subdirectory contains the two files needed to run the LLM engine for Biological interpretation.

This LLM bases its entire interpretation on the output from the QTL Mapping process so it is imperative 
to include the file "qtl_data.pkl" (which should be exported from R following the analysis) to get 
successful results. 

The engine functions in the following way:

1. Load all precomputed biological data
2. User selects group + phenotype in Streamlit
3. Build the gene table (core dataset)
4. Build gene context for the LLM
5. Identify enriched biological processes (GO terms)
6. Analyze genetic variants (SNPs)
7. Analyze founder strain effects
8. Retrieve relevant literature from PubMed
9. Aggregate literature into gene-level evidence
10. Generate final biological interpretation with LLM


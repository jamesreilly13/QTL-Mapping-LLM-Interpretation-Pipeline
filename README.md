# QTL-Mapping-LLM-Interpretation-Pipeline

Substance use disorders are driven by complex interactions between genetic, environmental, and behavioral factors, but the underlying gene networks remain poorly understood. To address this, we use Quantitative Trait Loci (QTL) mapping in Diversity Outbred (DO) mice, a highly recombinant population derived from eight founder strains that enables high-resolution analysis of polygenic traits.

We analyze 73 behavioral phenotypes across ~3,400 mice, including assays such as open field activity, light–dark box, holeboard exploration, and novelty place preference. Using an R/qtl2-based pipeline, we perform genome-wide association scans with kinship correction and permutation-based significance thresholds to identify trait-associated loci.

To move beyond locus discovery, we integrate:

- gene annotation
- variant-level analysis
- gene ontology (GO) enrichment

Finally, we introduce a retrieval-augmented LLM framework that combines QTL results with automated literature mining to generate biologically grounded interpretations. This system enables interactive exploration of genotype–phenotype relationships and supports hypothesis generation in addiction genetics.

The files for this process are contained in this repository in two subdirectories (LLM_engine and QTL_mapping_steps). QTL_mapping_steps is the first part of this process and LLM engine utlizes its outputs for analysis. Please utilize the README.txt files found within each directory for instructions.

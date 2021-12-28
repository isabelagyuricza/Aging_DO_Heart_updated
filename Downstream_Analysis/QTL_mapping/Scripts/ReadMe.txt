# This folder contains the scripts used to analyze the DO aging heart QTL mapping. QTL mapping was run using our tool at https://qtlviewer.jax.org/ and, for this reason, mapping scripts are not reported here. 

1. Residual_perm_on_cluster.R: First step of the analysis. Script runs age-interactive residual permutations on our HPC cluster to find significance thresholds for the age-QTLs.

2. Residual_permutation_analysis_mrna.R: Second step of the analysis. Script gathers the permuted lod scores for age-eQTLs to find the 95% quantile to use as significancy threshold.

3. Residual_permutation_analysis_protein.R: Third step of the analysis. Script gathers the permuted lod scores for age-pQTLs to find the 95% quantile to use as significancy threshold.

4. plot_transcriptome_age_int.R: Script plots the transcriptome map and the density plots for both distal age-eQTLs and age-pQTLs. 

5. Hotspot_chr3.R: Script runs the analysis on the age-pQTL hotspot on chromosome 3, including filtering and enrichment.

6. Hotspot_chr12.R: Script runs the analysis on the age-pQTL hotspot on chromosome 12, including filtering and enrichment.

7. Generating_supplemental_table.R: Script parses age-QTLs results and hotspot analysis to build a supplemental table with summaries for the QTL mapping. 

8. Kidney_vs_Heart.R: Script gathers QTL mapping summaries for both kidney and heart in order to compare the results.


# Isabela Gerdes Gyuricza - Churchill Lab - JAX # 
# Date: 12-28-2021 #
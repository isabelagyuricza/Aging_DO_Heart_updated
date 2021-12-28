# This folder contains the scripts used to analyze the DO aging heart. 

1. DEA_and_enrichment.R: First step of the analysis. Script runs the differential expression and the enrichment analysis on both transcripts and proteins. 

2. gather_enrichments.R: Second step of the analysis. Script combines the enrichment results of transcripts and proteins based on gene IDs. 

# The scripts below can be run independently #

3. proteincomplex_corplot.R: Script generates age and sex-specific correlation plots for multiple protein complexes to evaluate how the correlation changes in the heart.

4. proteincomplex_overall_regression.R: Scripts computes the overall age effect for each protein complex to check how age affects changes in correlation in the heart.

5. proteincomplex_regression.R: Scripts computes the age effect on the correlation of each gene-pair member of protein complexes in the heart.

6. KIDNEY_proteincomplex_corplot.R: Script generates age and sex-specific correlation plots for multiple protein complexes to evaluate how the correlation changes in the kidney.

7. KIDNEY_proteincomplex_overall_regression.R: Scripts computes the overall age effect for each protein complex to check how age affects changes in correlation in the kidney.

8. KIDNEY_proteincomplex_regression.R: Scripts computes the age effect on the correlation of each gene-pair member of protein complexes in the kidney.

9. Transcript_Protein_concordance: Script plots the concordance between the age/sex effects of transcripts and their corresponding proteins. 

10. Plotting_selected_genes.R: Script generates box plots and scatter plots for the expression  of selected transcripts and proteins.

11. Compare_Heart_Kidney.R: Script compares the transcript/protein age effects between heart and kidney and run enrichment analysis on the concordant and discordant groups.

12. Enrichment_functions.R: Script defines functions to summarize results from the enrichment analysis by ClusterProfiler.

13. Variance_RNA_Protein.R: Script computes age-specific variances for the expression/abundance of each transcript and protein in our data.


# Isabela Gerdes Gyuricza - Churchill Lab - JAX # 
# Date: 12-28-2021 #
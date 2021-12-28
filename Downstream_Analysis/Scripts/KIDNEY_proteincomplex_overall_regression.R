##################### Protein complexes analysis overall #######################

# KIDNEY

# This script uses the information on Greg's table 
# to gather the genes and proteins of each complex
# and check how their correlation is affected by age.

# Gathering the correlation values into protein complexes, applying the regression
# and permutting..

# Author: Isabela Gerdes Gyuricza
# Date: 07_28_2021

################################################################################
############ clear workspace

rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

################################################################################
############ libraries

library(tidyverse)
library(DGCA)

################################################################################
############ load data

#setwd("~/Box/JAC_Heart_Data/Heart_Data_June2021")

# Load QTLviewer data 

load("Downstream_Analysis/Data/JAC_DO_kidney_v5_03_21_2020.RData")

# Cleaning up what I don't need

rm(genoprobs, K, map, markers, ensembl.version,dataset.protein.zam)

complex_info <- readRDS(
  "Downstream_Analysis/Data/ori_protein_complex_long_dat.RDS")

complex_info <- complex_info %>% 
  select(complex_name,mouse_ensembl_id,symbol) %>% 
  unique()

complex_info <- complex_info %>% 
  rename(gene.id = mouse_ensembl_id) %>% 
  inner_join(dataset.protein.zaz$annot.protein %>% 
               select(gene.id, protein.id), 
             by = "gene.id") %>% 
  select(complex_name,gene.id,protein.id,symbol) %>% 
  mutate(complex_name = gsub(" ","_",complex_name))

# Keeping only complexes with 4 members or more 

member_count <- complex_info %>% 
  group_by(complex_name) %>% 
  tally %>% 
  filter(n >= 4)

# Filtering the dataset

complex_info <- complex_info %>% 
  filter(complex_name %in% member_count$complex_name)

# We have 187 total protein complexes that match the criteria

rm(member_count)

# Adjusting the mitochondrial complexes names 

complex_info$complex_name[
  which(
    complex_info$complex_name == "F0/F1_ATP_synthase_(complex_V)")] <- "Mitochondrial_complex_V"

complex_info$complex_name[
  which(
    complex_info$complex_name == "Cytochrome_bc1_complex_(Ubiquinol-cytochrome_c_reductase_complex,_complex_III)")] <- "Mitochondrial_complex_III"

complex_info$complex_name[
  which(
    complex_info$complex_name == "Cytochrome_c_oxidase_(complex_IV)")] <- "Mitochondrial_complex_IV"

complex_info$complex_name[
  which(
    complex_info$complex_name == "mitochondrial_inner_membrane_presequence_translocase_complex")] <- "Mitochondrial_inner_membrane_translocase"

complex_info$complex_name[
  which(
    complex_info$complex_name == "mitochondrial_outer_membrane_translocase_complex")] <- "Mitochondrial_outer_membrane_translocase"

complex_info$complex_name[
  which(
    complex_info$complex_name == "Respiratory_chain_complex_I_(early_intermediate_NDUFAF1_assembly),__mitochondrial")] <- "Mitochondrial_complex_I"

# I could not find Mitochondrial Complex II, creating my own..

genes_info <- ensimplR::batchGenes(c("ENSMUSG00000021577","ENSMUSG00000009863",
                                     "ENSMUSG00000058076","ENSMUSG00000000171"), 
                                             species = "Mm", release = 94)


more_info <- data_frame(complex_name = "Mitochondrial_complex_II",
                        gene.id = genes_info$id)

more_info <- more_info %>% 
  left_join(dataset.protein.zaz$annot.protein, by = "gene.id") %>% 
  select(complex_name,gene.id,protein.id,symbol)

complex_info <- rbind(complex_info,more_info)

protein_expression <- data.frame(dataset.protein.zaz$data$norm[,
                                                               which(colnames(dataset.protein.zaz$data$norm)
                                                                     %in% unique(complex_info$protein.id))]) %>% 
  rownames_to_column("Mouse.ID") %>% 
  gather("protein.id","protein_expression",-Mouse.ID)

gene_expression <- data.frame(dataset.mrna$data$norm[,
                                                     which(colnames(dataset.mrna$data$norm)
                                                           %in% unique(complex_info$gene.id))]) %>% 
  rownames_to_column("Mouse.ID") %>% 
  gather("gene.id","gene_expression",-Mouse.ID)

df <- complex_info %>% 
  left_join(protein_expression, by = "protein.id") %>% 
  left_join(gene_expression, by = c("gene.id","Mouse.ID")) %>% 
  left_join(dataset.mrna$annot.sample %>% 
              select(Mouse.ID, Age, Sex), by = "Mouse.ID") %>% 
  mutate(Age = (Age/12)-1) #Rescaling age to the same as DEA

# Correcting symbols, so that aren't any repeated ones (Adding .1 to duplicates)

symbols <- df %>% 
  select(protein.id, symbol) %>% 
  unique() %>% 
  mutate(gene_symbol = make.unique(symbol)) %>% 
  select(protein.id, gene_symbol)

df <- df %>% 
  left_join(symbols, by = "protein.id") %>% 
  select(complex_name, gene.id, protein.id, gene_symbol, Mouse.ID, Age, Sex,
         protein_expression, gene_expression) %>% 
  rename(mouse.id = Mouse.ID)

rm(more_info, gene_expression, protein_expression, complex_info, 
   genes_info, dataset.mrna, dataset.protein.zaz, symbols)

################## Fitting regression for each complex #########################

# 1) Proteins

age_trend_real_protein_list <- list()

for (i in unique(df$complex_name)){
  
  df_test_6 <- df %>% 
    filter(complex_name == i,
           Age == -0.5) %>% 
    select(mouse.id,gene_symbol,protein_expression) %>% 
    spread(gene_symbol,protein_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_12 <- df %>% 
    filter(complex_name == i,
           Age == 0) %>% 
    select(mouse.id, gene_symbol,protein_expression) %>% 
    spread(gene_symbol,protein_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_18 <- df %>% 
    filter(complex_name == i,
           Age == 0.5) %>% 
    select(mouse.id, gene_symbol,protein_expression) %>% 
    spread(gene_symbol,protein_expression) %>% 
    column_to_rownames("mouse.id")
  
  COR_06 <- cor(df_test_6, use = "pairwise.complete.obs")
  COR_12 <- cor(df_test_12, use = "pairwise.complete.obs")
  COR_18 <- cor(df_test_18, use = "pairwise.complete.obs")
  
  names <- which( upper.tri(COR_06,diag=F) , arr.ind = TRUE )
  
  df_cor <- tibble(cor_6 = COR_06[upper.tri(COR_06)],
                   cor_12 = COR_12[upper.tri(COR_12)],
                   cor_18 = COR_18[upper.tri(COR_18)],
                   comparing_this = dimnames(COR_06)[[1]][names[,1]],
                   to_this = dimnames(COR_06)[[2]][names[,2]])
  
  df_cor <- df_cor %>% 
    unite(col = "pair_comparison", c("comparing_this","to_this"),sep = "_") %>% 
    gather("Age","cor",-pair_comparison) %>% 
    mutate(Age = gsub("cor_","",Age)) %>% 
    mutate(Age = as.numeric(Age)) %>% 
    mutate(Age = (Age/12)-1)
  
  
  MODEL <- suppressMessages(lmerTest::lmer(
    cor ~ (1 + Age | pair_comparison) + Age,
    data = df_cor))
  
  print(lmerTest::ranova(MODEL))
  
  FIXEF <- lme4::fixef(MODEL) 
  
  coefs <- summary(MODEL)$coefficients
  
  coefs <- coefs["Age",which(colnames(coefs) %in% c("Std. Error","Pr(>|t|)"))]
  
  names(coefs) <- c("SE","pvalue_model")
  
  AGE_EFFECT <- FIXEF[2]
  
  age_trend_real_protein_list[[i]] <- c(AGE_EFFECT, coefs)
  
}

rm(i,df_test_12,df_test_18,df_test_6,COR_06,COR_12,COR_18,df_cor,
   names, AGE_EFFECT, coefs, FIXEF,MODEL)

# Creating a function to do that permuting mice 1000 times, disrupting only 
#the mouse - age association.

permut_age_trend <- function(df, times, seed = 01) {
  
  set.seed(seed = seed)
  age_trend_permut_protein_list <- list()
  
  for (i in unique(df$complex_name)){
    
    age_trend_perm <- vector(mode = 'numeric', length = times)
    
    for (j in 1:times) {
      
      df_test_6 <- df %>% 
        filter(complex_name == i,
               mouse.id %in% sample(unique(df$mouse.id),62)) %>% 
        select(mouse.id,gene_symbol,protein_expression) %>% 
        spread(gene_symbol,protein_expression) %>% 
        column_to_rownames("mouse.id")
      
      df_test_12 <- df %>% 
        filter(complex_name == i,
               mouse.id %in% sample(unique(df$mouse.id[
                 !df$mouse.id %in% rownames(df_test_6)
               ]),62)) %>% 
        select(mouse.id, gene_symbol,protein_expression) %>% 
        spread(gene_symbol,protein_expression) %>% 
        column_to_rownames("mouse.id")
      
      df_test_18 <- df %>% 
        filter(complex_name == i,
               !mouse.id %in% rownames(df_test_6) &
                 !mouse.id %in% rownames(df_test_12)) %>% 
        select(mouse.id, gene_symbol,protein_expression) %>% 
        spread(gene_symbol,protein_expression) %>% 
        column_to_rownames("mouse.id")
      
      COR_06 <- cor(df_test_6, use = "pairwise.complete.obs")
      COR_12 <- cor(df_test_12, use = "pairwise.complete.obs")
      COR_18 <- cor(df_test_18, use = "pairwise.complete.obs")
      
      names <- which( upper.tri(COR_06,diag=F) , arr.ind = TRUE )
      
      df_cor <- tibble(cor_6 = COR_06[upper.tri(COR_06)],
                       cor_12 = COR_12[upper.tri(COR_12)],
                       cor_18 = COR_18[upper.tri(COR_18)],
                       comparing_this = dimnames(COR_06)[[1]][names[,1]],
                       to_this = dimnames(COR_06)[[2]][names[,2]])
      
      df_cor <- df_cor %>% 
        unite(col = "pair_comparison", c("comparing_this","to_this"),sep = "_") %>% 
        gather("Age","cor",-pair_comparison) %>% 
        mutate(Age = gsub("cor_","",Age))  %>% 
        mutate(Age = as.numeric(Age)) %>% 
        mutate(Age = (Age/12)-1)
      
      MODEL <- suppressMessages(lme4::lmer(
        cor ~ (1 + Age | pair_comparison) + Age,
        data = df_cor))
      
      FIXEF <- lme4::fixef(MODEL)
      
      AGE_EFFECT <- FIXEF[2]
      
      age_trend_perm[j] <- as.numeric(AGE_EFFECT)
      
    }
    
    age_trend_permut_protein_list[[i]] <- age_trend_perm
    
    
  }
  
  return(age_trend_permut_protein_list)
}

age_trend_permut_protein_results <- permut_age_trend(df,1000)

# save(age_trend_permut_protein_results, 
#      file = "Downstream_Analysis/Results/KIDNEY_proteincomplex_overall_perms_protein.RData")


# Using bigEmpPVals function from DGCA to compute the pvals from the permutation

empirical_pvals <- list() 

for (i in names(age_trend_permut_protein_results)){
  
  complex_df_perm <- age_trend_permut_protein_results[[i]]
  
  complex_df_real <- age_trend_real_protein_list[[i]][[1]]
  
  
  pval <- bigEmpPVals(stat = abs(complex_df_real), 
                      stat0 = abs(complex_df_perm))
  
  empirical_pvals[[i]] <- pval
  
}

rm(complex_df_perm, complex_df_real,i,pval)

empirical_pvals <- bind_rows(empirical_pvals) %>% 
  gather("Complex_name","p_Protein")

final_df <- bind_rows(age_trend_real_protein_list, .id = "Complex_name") %>% 
  rename(Age_effect_Protein = Age,
         pvalue_model_Protein = pvalue_model,
         Age_effect_SE_Protein = SE) %>% 
  left_join(empirical_pvals, by = "Complex_name") %>% 
  mutate(adjusted_p_Protein = adjustPVals(p_Protein, adjust = "BH"))

rm(empirical_pvals)

final_df %>% 
  filter(adjusted_p_Protein < 0.1) %>% 
  arrange(desc(abs(Age_effect_Protein))) %>% 
  View()

# There are 43 protein complexes that change their overall correlation with age, 
# and they mostly increase with age.

################################################################################


# 1) Genes

age_trend_real_gene_list <- list()

for (i in unique(df$complex_name)){
  
  df_test_6 <- df %>% 
    filter(complex_name == i,
           Age == -0.5) %>% 
    select(mouse.id,gene_symbol,gene_expression) %>% 
    unique() %>% 
    spread(gene_symbol,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_12 <- df %>% 
    filter(complex_name == i,
           Age == 0) %>% 
    select(mouse.id, gene_symbol,gene_expression) %>% 
    unique() %>% 
    spread(gene_symbol,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_18 <- df %>% 
    filter(complex_name == i,
           Age == 0.5) %>% 
    select(mouse.id, gene_symbol,gene_expression) %>% 
    unique() %>% 
    spread(gene_symbol,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  COR_06 <- cor(df_test_6, use = "pairwise.complete.obs")
  COR_12 <- cor(df_test_12, use = "pairwise.complete.obs")
  COR_18 <- cor(df_test_18, use = "pairwise.complete.obs")
  
  names <- which( upper.tri(COR_06,diag=F) , arr.ind = TRUE )
  
  df_cor <- tibble(cor_6 = COR_06[upper.tri(COR_06)],
                   cor_12 = COR_12[upper.tri(COR_12)],
                   cor_18 = COR_18[upper.tri(COR_18)],
                   comparing_this = dimnames(COR_06)[[1]][names[,1]],
                   to_this = dimnames(COR_06)[[2]][names[,2]])
  
  df_cor <- df_cor %>% 
    unite(col = "pair_comparison", c("comparing_this","to_this"),sep = "_") %>% 
    gather("Age","cor",-pair_comparison) %>% 
    mutate(Age = gsub("cor_","",Age)) %>% 
    mutate(Age = as.numeric(Age)) %>% 
    mutate(Age = (Age/12)-1)
  
  MODEL <- suppressMessages(lmerTest::lmer(                       
    cor ~ (1 + Age | pair_comparison) + Age,
    data = df_cor))
  
  print(lmerTest::ranova(MODEL))
  
  FIXEF <- lme4::fixef(MODEL)
  
  coefs <- summary(MODEL)$coefficients
  
  coefs <- coefs["Age",which(colnames(coefs) %in% c("Std. Error","Pr(>|t|)"))]
  
  names(coefs) <- c("SE","pvalue_model")
  
  
  AGE_EFFECT <- FIXEF[2]
  
  age_trend_real_gene_list[[i]] <- c(AGE_EFFECT, coefs) 
}

rm(i,df_test_12,df_test_18,df_test_6,COR_06,COR_12,COR_18,df_cor,
   names, MODEL, FIXEF, coefs, AGE_EFFECT)

# Creating a function to do that permuting mice 1000 times, disrupting only 
#the mouse - age association.

permut_age_trend <- function(df, times, seed = 10) {
  
  set.seed(seed = seed)
  age_trend_permut_gene_list <- list()
  
  for (i in unique(df$complex_name)){
    
    age_trend_perm <- vector(mode = 'numeric', length = times)
    
    for (j in 1:times) {
      
      df_test_6 <- df %>% 
        filter(complex_name == i,
               mouse.id %in% sample(unique(df$mouse.id),62)) %>% 
        select(mouse.id,gene_symbol,gene_expression) %>% 
        unique() %>% 
        spread(gene_symbol,gene_expression) %>% 
        column_to_rownames("mouse.id")
      
      df_test_12 <- df %>% 
        filter(complex_name == i,
               mouse.id %in% sample(unique(df$mouse.id[
                 !df$mouse.id %in% rownames(df_test_6)
               ]),62)) %>% 
        select(mouse.id, gene_symbol,gene_expression) %>% 
        unique() %>% 
        spread(gene_symbol,gene_expression) %>% 
        column_to_rownames("mouse.id")
      
      df_test_18 <- df %>% 
        filter(complex_name == i,
               !mouse.id %in% rownames(df_test_6) &
                 !mouse.id %in% rownames(df_test_12)) %>% 
        select(mouse.id, gene_symbol,gene_expression) %>% 
        unique() %>% 
        spread(gene_symbol,gene_expression) %>% 
        column_to_rownames("mouse.id")
      
      COR_06 <- cor(df_test_6, use = "pairwise.complete.obs")
      COR_12 <- cor(df_test_12, use = "pairwise.complete.obs")
      COR_18 <- cor(df_test_18, use = "pairwise.complete.obs")
      
      names <- which( upper.tri(COR_06,diag=F) , arr.ind = TRUE )
      
      df_cor <- tibble(cor_6 = COR_06[upper.tri(COR_06)],
                       cor_12 = COR_12[upper.tri(COR_12)],
                       cor_18 = COR_18[upper.tri(COR_18)],
                       comparing_this = dimnames(COR_06)[[1]][names[,1]],
                       to_this = dimnames(COR_06)[[2]][names[,2]])
      
      df_cor <- df_cor %>% 
        unite(col = "pair_comparison", c("comparing_this","to_this"),sep = "_") %>% 
        gather("Age","cor",-pair_comparison) %>% 
        mutate(Age = gsub("cor_","",Age)) %>% 
        mutate(Age = as.numeric(Age)) %>% 
        mutate(Age = (Age/12)-1)
      
      MODEL <- suppressMessages(lme4::lmer(
        cor ~ (1 + Age | pair_comparison) + Age,
        data = df_cor ))
      
      FIXEF <- lme4::fixef(MODEL)
      
      AGE_EFFECT <- FIXEF[2]
      
      age_trend_perm[j] <- as.numeric(AGE_EFFECT)
      
    }
    
    age_trend_permut_gene_list[[i]] <- age_trend_perm
    
    
  }
  
  return(age_trend_permut_gene_list)
}

age_trend_permut_gene_results <- permut_age_trend(df,1000)


# save(age_trend_permut_gene_results, 
#      file = "Downstream_Analysis/Results/proteincomplex_overall_perms_transcripts.RData")
# 

# Using bigEmpPVals function from DGCA to compute the pvals from the permutation

empirical_pvals <- list() 

for (i in names(age_trend_permut_gene_results)){
  
  complex_df_perm <- age_trend_permut_gene_results[[i]]
  
  complex_df_real <- age_trend_real_gene_list[[i]][[1]]
  
  
  pval <- bigEmpPVals(stat = abs(complex_df_real), 
                      stat0 = abs(complex_df_perm))
  
  empirical_pvals[[i]] <- pval
  
}

empirical_pvals <- bind_rows(empirical_pvals) %>% 
  gather("Complex_name","p_Transcript")


final_df_2 <- bind_rows(age_trend_real_gene_list, .id = "Complex_name") %>% 
  rename(Age_effect_Transcript = Age,
         pvalue_model_Transcript = pvalue_model,
         Age_effect_SE_Transcript = SE) %>% 
  left_join(empirical_pvals, by = "Complex_name") %>% 
  mutate(adjusted_p_Transcript = adjustPVals(p_Transcript, adjust = "BH"))

# There are no significant overall changes at the transcript level!

final_df <- final_df_2 %>% 
  left_join(final_df, by = "Complex_name") %>% 
  select(Complex_name, 
         Age_effect_Transcript,
         Age_effect_SE_Transcript,
         p_Transcript,
         adjusted_p_Transcript,
         Age_effect_Protein,
         Age_effect_SE_Protein,
         p_Protein,
         adjusted_p_Protein)

# Create dataframe for plotting and considering permutations adjusted p-values
# as threshold (FDR < 0.1)

final_df <- final_df %>% 
  mutate(STD_Age_effect_Transcript = Age_effect_Transcript/Age_effect_SE_Transcript,
         STD_Age_effect_Protein = Age_effect_Protein/Age_effect_SE_Protein,
         Signif_Protein = ifelse(adjusted_p_Protein < 0.1, TRUE,FALSE),
         Signif_Transcript = ifelse(adjusted_p_Transcript < 0.1, TRUE,FALSE))

# Now, trying to add the number of subunits for each complex

df_new <- df %>% 
  select(complex_name, protein.id) %>% 
  unique() %>% 
  group_by(complex_name) %>% 
  summarise(size = length(protein.id)) %>% 
  rename(Complex_name = complex_name) %>% 
  inner_join(final_df, by = "Complex_name")

# Save for supplemental material  

# write.csv(df_new %>% 
#             select(-STD_Age_effect_Transcript, -STD_Age_effect_Protein,
#                    -Signif_Transcript, -Signif_Protein) %>%
#             rename(Size = size),
#           file = "Downstream_Analysis/Results/KIDNEY_proteincomplex_overall_age_effects.csv",
#           row.names = FALSE)

pdf("Downstream_Analysis/Results/KIDNEY_proteincomplex_STD_overall_age_effects_2.pdf",
    width = 13,height = 10)
#quartz()
df_new %>%
  mutate(
    `Significance (FDR < 0.1)` = ifelse(
      Signif_Transcript == TRUE & Signif_Protein == FALSE,
      'Transcript Only',
      ifelse(
        Signif_Transcript == FALSE & Signif_Protein == TRUE, 
        'Protein Only',
        ifelse(
          Signif_Transcript == TRUE & Signif_Protein == TRUE,
          'Both',
          'Neither'
        )
      )
    ),
    `Significance (FDR < 0.1)` = factor(
      `Significance (FDR < 0.1)`,
      levels = c('Transcript Only', 'Protein Only', 'Both', 'Neither')
    )
  ) %>% 
  ggplot(aes(x = STD_Age_effect_Transcript, y = STD_Age_effect_Protein)) +
  geom_point(aes(colour = `Significance (FDR < 0.1)`, size = size), alpha = 0.7) + 
  scale_size(name = "Size", range = c(1,15)) +
  geom_text(data = df_new %>% 
              filter(abs(STD_Age_effect_Protein) > 2 |
                       abs(STD_Age_effect_Transcript) > 2),
            aes(label = Complex_name), hjust = 0, vjust = 0, size = 1) +
  scale_x_continuous(
    lim = c(-20,20),
    breaks = c(-20, -10,0,10,20)) +
  scale_y_continuous(
    lim = c(-20, 20),
    breaks = c(-20, -10, 0, 10, 20)) +
  scale_colour_manual(breaks = c('Transcript Only', 'Protein Only', 'Both', 'Neither'),
                      values = c("firebrick","royalblue","mediumorchid4","gray")) +
  geom_vline(xintercept = 0,linetype = "dashed") +
  geom_hline(yintercept = 0,linetype = "dashed") +
  theme_minimal() 

dev.off()

################################################################################

# Now, plotting the age effect against the size of the protein complex.

# rm(list = setdiff(ls(), "df"))
# 
# results <- read_csv("Downstream_Analysis/Results/KIDNEY_proteincomplex_overall_age_effects.csv")
# 
# df_new <- df %>% 
#   select(complex_name, protein.id) %>% 
#   unique() %>% 
#   group_by(complex_name) %>% 
#   summarise(size = length(protein.id)) %>% 
#   rename(Complex_name = complex_name) %>% 
#   inner_join(results) %>% 
#   mutate(STD_Age_effect_Transcript = Age_effect_Transcript/Age_effect_SE_Transcript,
#          STD_Age_effect_Protein = Age_effect_Protein/Age_effect_SE_Protein)
# 
# pdf("Downstream_Analysis/Results/KIDNEY_proteincomplex_size_vs_effect.pdf",
#     width = 14,height =4)
# 
# plot1 <- df_new %>% 
#   mutate(`Significant (FDR < 0.1)` = ifelse(adjusted_p_Transcript < 0.1,
#                                             TRUE,
#                                             FALSE)) %>% 
#   ggplot(aes(x = Age_effect_Transcript, y = size)) +
#   geom_point(aes(color = `Significant (FDR < 0.1)`), size = 5, alpha = 0.7) + 
#   geom_text(aes(label = ifelse(size > 15, Complex_name,
#                                "")), hjust = 0, vjust = 0, size = 3) +
#   scale_colour_manual(breaks = c('TRUE',"FALSE"),
#                       values = c("firebrick","gray")) +
#   theme_bw() +
#   ggtitle("Transcripts")
# 
# plot2 <- df_new %>%
#   mutate(`Significant (FDR < 0.1)` = ifelse(adjusted_p_Protein < 0.1,
#                                             TRUE,
#                                             FALSE)) %>% 
#   ggplot(aes(x = Age_effect_Protein, y = size)) +
#   geom_point(aes(color = `Significant (FDR < 0.1)`), size = 5, alpha = 0.7) + 
#   geom_text(aes(label = ifelse(size > 14, Complex_name,
#                                "")), hjust = 0, vjust = 0, size = 3) +
#   scale_colour_manual(breaks = c('TRUE',"FALSE"),
#                       values = c("royalblue","gray")) +
#   theme_bw() +
#   ggtitle("Proteins")
# 
# cowplot::plot_grid(plot1,plot2, ncol = 2)
# 
# dev.off()
# 

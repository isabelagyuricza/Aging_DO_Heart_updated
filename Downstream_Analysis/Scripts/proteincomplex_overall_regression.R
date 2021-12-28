##################### Protein complexes analysis overall #######################

# This script uses the information on Ori's table 
# to gather the genes and proteins of each complex
# and check how their correlation is affected by age.

# Gathering the correlation values into protein complexes, applying the regression
# and permutting..

# Author: Isabela Gerdes Gyuricza
# Date: 07_20_2021

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

load("QTLviewer_dataset/JAC_DO_heart_v9.gz")

# Cleaning up what I don't need

rm(genoprobs, K, map, markers, ensembl.version)

complex_info <- readRDS(
  "Downstream_Analysis/Data/ori_protein_complex_long_dat.RDS")

complex_info <- complex_info %>% 
  select(complex_name,mouse_ensembl_id,symbol) %>% 
  unique()

complex_info <- complex_info %>% 
  rename(gene.id = mouse_ensembl_id) %>% 
  inner_join(dataset.protein$annot.protein %>% 
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

# We have 120 total protein complexes that match the criteria

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

genes_info <- dataset.protein$annot.protein %>% 
  filter(symbol %in% c("Sdha","Sdhb","Sdhc","Sdhd"))

more_info <- data_frame(complex_name = "Mitochondrial_complex_II",
                        gene.id = genes_info$gene.id)

# Also, I noticed that there are some genes missing from the proteasome, making 
# sure we have all of them

match <- c("Psma", "Psmb","Psmc","Psmd","Psme","Psmf","Psmg")

genes_info <-  dataset.protein$annot.protein %>% 
  filter(grepl(paste(match,collapse = "|"),symbol))

more_info <- rbind(more_info, 
                   data_frame(complex_name = "26S_Proteasome",
                              gene.id = genes_info$gene.id))

more_info <- more_info %>% 
  left_join(dataset.protein$annot.protein, by = "gene.id") %>% 
  select(complex_name,gene.id,protein.id,symbol)

complex_info <- complex_info %>% 
  full_join(more_info, by = c("complex_name","gene.id","protein.id","symbol")) %>% 
  unique()

protein_expression <- data.frame(dataset.protein$data$norm[,
                                                           which(colnames(dataset.protein$data$norm)
                                                                 %in% unique(complex_info$protein.id))]) %>% 
  rownames_to_column("mouse.id") %>% 
  gather("protein.id","protein_expression",-mouse.id)

gene_expression <- data.frame(dataset.mrna$data$norm[,
                                                     which(colnames(dataset.mrna$data$norm)
                                                           %in% unique(complex_info$gene.id))]) %>% 
  rownames_to_column("mouse.id") %>% 
  gather("gene.id","gene_expression",-mouse.id)

df <- complex_info %>% 
  left_join(protein_expression, by = "protein.id") %>% 
  left_join(gene_expression, by = c("gene.id","mouse.id")) %>% 
  left_join(dataset.mrna$annot.sample %>% 
              select(mouse.id, Age, Sex), by = "mouse.id") %>% 
  mutate(Age = (Age/12)-1) #Rescaling age to the same as DEA

# Correcting symbols, so that aren't any repeated ones (Adding .1 to duplicates)

symbols <- df %>% 
  select(protein.id, symbol) %>% 
  unique() %>% 
  mutate(gene_symbol = make.unique(symbol)) %>% 
  select(protein.id, gene_symbol)

df <- df %>% 
  left_join(symbols, by = "protein.id") %>% 
  select(complex_name, gene.id, protein.id, gene_symbol, mouse.id, Age, Sex,
         protein_expression, gene_expression) %>% 
  unique()

rm(more_info, gene_expression, protein_expression, complex_info, 
   genes_info, dataset.mrna, dataset.protein, symbols,match)

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

permut_age_trend <- function(df, times, seed = 19940418) {
  
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
#      file = "Downstream_Analysis/Results/proteincomplex_overall_perms_protein.RData")


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

# There are 38 protein complexes that change their overall correlation with age, 
# al they all decrease!

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

permut_age_trend <- function(df, times, seed = 20200701) {
  
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
#   file = "Downstream_Analysis/Results/proteincomplex_overall_age_effects.csv",
#   row.names = FALSE)

pdf("Downstream_Analysis/Results/proteincomplex_STD_overall_age_effects_2.pdf",
    width = 8,height = 7)
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
              filter(STD_Age_effect_Protein < -8),
            aes(label = Complex_name), hjust = 0, vjust = 0, size = 2) +
  scale_x_continuous(
    lim = c(-35,35),
    breaks = seq(-200, 200, 20),
    minor_breaks =  seq(-200, 200, 10)
  ) +
  scale_y_continuous(
    lim = c(-35, 35),
    breaks = seq(-200, 200, 20),
    minor_breaks = seq(-200, 200, 10)
  ) +
  scale_colour_manual(breaks = c('Transcript Only', 'Protein Only', 'Both', 'Neither'),
                      values = c("firebrick","royalblue","mediumorchid4","gray")) +
  geom_vline(xintercept = 0,linetype = "dashed") +
  geom_hline(yintercept = 0,linetype = "dashed") +
  theme_minimal() +
  theme(legend.position = "none")

dev.off()

################################################################################

# Checking if global expression changes with age, probably not. Gathering the 
# LFC and slope for each subunit of the protein complexes that show the most
# pronounced overall changes (26S_Proteasome, Large Drosha complex, 
# Chaperon containing T complex, nuclear pore complex and
# cytoplasmic ribosomal large subunit.)

rm(list=setdiff(ls(), "df"))

DEA_transcripts <- read_csv("Downstream_Analysis/Results/DEA_age_transcripts.csv")

DEA_proteins <- read_csv("Downstream_Analysis/Results/DEA_age_proteins.csv")


df_new <- df %>%
  left_join(DEA_transcripts %>% 
              rename(gene.id = Gene_ID) %>% 
              select(gene.id, Effect, Effect_SE, adjusted_p)) %>% 
  rename(Effect_transcript = Effect, 
         SE_transcript = Effect_SE,
         pvalue_transcript = adjusted_p) %>% 
  left_join(DEA_proteins %>% 
              rename(protein.id = Protein_ID) %>% 
              select(protein.id, Effect, Effect_SE, adjusted_p)) %>% 
  rename(Effect_protein = Effect, 
         SE_protein = Effect_SE,
         pvalue_protein = adjusted_p) %>% 
  select(complex_name, gene_symbol, Effect_transcript, Effect_protein, 
         SE_transcript, SE_protein, pvalue_transcript, pvalue_protein) %>% 
  unique() %>% 
  mutate(STD_effect_transcript = Effect_transcript/SE_transcript,
         STD_effect_protein = Effect_protein/SE_protein) %>% 
  gather("Type_effect","Age_effect", Effect_transcript, Effect_protein) %>% 
  gather("Type_STD_effect","STD_age_effect",STD_effect_transcript, 
         STD_effect_protein) %>% 
  gather("Type_SE","SE", SE_transcript, SE_protein) %>% 
  gather("Type_pvalue","pvalue", pvalue_transcript, pvalue_protein) %>% 
  mutate(Type_effect = gsub("Effect_","", Type_effect),
         Type_STD_effect = gsub("STD_effect_","", Type_STD_effect),
         Type_SE = gsub("SE_","", Type_SE),
         Type_pvalue = gsub("pvalue_","", Type_pvalue)) %>% 
  filter(Type_effect == Type_STD_effect & 
           Type_STD_effect == Type_SE &
           Type_SE == Type_pvalue) %>% 
  rename(Type = Type_effect) %>% 
  select(complex_name, gene_symbol, Type, Age_effect, SE, STD_age_effect, pvalue) %>% 
  mutate(log10pvalue = -log10(pvalue)) %>% 
  filter(complex_name %in% c("26S_Proteasome",
                             "Nuclear_pore_complex_(NPC)",
                             "cytoplasmic_ribosomal_large_subunit",
                             "chaperonin-containing_T-complex",
                             "Large_Drosha_complex"))
 

pdf("Downstream_Analysis/Results/proteincomplex_age_effects.pdf", 
    width = 9, 
    height = 5)

df_new %>% 
  filter(Type == "transcript") %>% 
  ggplot(aes(x = Age_effect, y = log10pvalue)) +
  geom_point(aes(color = complex_name), size = 3, alpha = 0.7) +
  geom_vline(xintercept = 0, col="gray", linetype = "dashed", size = 1) +
  geom_hline(yintercept = -log10(0.01), col="red", linetype = "dashed", size = 1) +
  geom_text(aes(label = ifelse(complex_name == "26S_Proteasome" &
                                 gene_symbol %in% c("Psmb8","Psmb9",
                                                    "Psmb3","Psmd7"),
                               gene_symbol, "")), size = 3) +
  xlim(c(-0.55, 0.55)) +
  labs(title = "Transcript level",
       x = "Age_effect", 
       y = "-log10(p-value)") +
  scale_color_manual(values = c(
    '#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00')) +
  theme_bw()
  
df_new %>% 
  filter(Type == "protein") %>% 
  ggplot(aes(x = Age_effect, y = log10pvalue)) +
  geom_point(aes(color = complex_name), size = 3, alpha = 0.7) +
  geom_vline(xintercept = 0, col="gray", linetype = "dashed", size = 1) +
  geom_hline(yintercept = -log10(0.01), col="red", linetype = "dashed", size = 1) +
  geom_text(aes(label = ifelse(complex_name == "26S_Proteasome" &
                                 log10pvalue > -log10(0.01),
                               gene_symbol, "")), size = 3) +
  xlim(c(-1, 1)) +
  labs(title = "Protein level",
       x = "Age_effect", 
       y = "-log10(p-value)") +
  scale_color_manual(values = c(
    '#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','black')) +
  theme_bw()

dev.off()

################################################################################

# Now, plotting the age effect against the size of the protein complex.

# rm(list = setdiff(ls(), "df"))
# 
# results <- read_csv("Downstream_Analysis/Results/proteincomplex_overall_age_effects.csv")
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
# pdf("Downstream_Analysis/Results/proteincomplex_size_vs_effect.pdf",
#     width = 12,height =4)
# 
# plot1 <- df_new %>%
#   ggplot(aes(x = Age_effect_Transcript, y = size)) +
#   geom_point(color = "gray", size = 5, alpha = 0.7) + 
#   geom_text(aes(label = ifelse(size > 15, Complex_name,
#                                "")), hjust = 0, vjust = 0, size = 3) +
#   theme_bw() +
#   ggtitle("Transcripts")
# 
# plot2 <- df_new %>%
#   mutate(`Significant (FDR < 0.1)` = ifelse(adjusted_p_Protein < 0.1,
#                                             TRUE,
#                                             FALSE)) %>% 
#   ggplot(aes(x = Age_effect_Protein, y = size)) +
#   geom_point(aes(color = `Significant (FDR < 0.1)`), size = 5, alpha = 0.7) + 
#   geom_text(aes(label = ifelse(size > 11, Complex_name,
#                                "")), hjust = 0, vjust = 0, size = 3) +
#   scale_colour_manual(breaks = c('TRUE',"FALSE"),
#                       values = c("royalblue","gray")) +
#   theme_bw() +
#   ggtitle("Proteins")
# 
# cowplot::plot_grid(plot1,plot2, ncol = 2,
#                    rel_widths = c(1.5,2))
# 
# dev.off()
# 

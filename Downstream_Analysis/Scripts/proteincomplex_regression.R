################### Protein complexes analysis - regression ####################

# Running permutations on changes in correlation with age to find significant 
# differences

# Author: Isabela Gerdes Gyuricza
# Date: 06_07_2021

################################################################################
############ clear workspace

rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

################################################################################
############ libraries

library(DGCA)
library(tidyverse)

################################################################################
############ load data

setwd("~/Box/JAC_Heart_Data/Heart_Data_June2021")

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

##### Overall regression for each pair_comparison within a protein complex #####
################################################################################

est_age_trend <- function(x, age){
  MODEL <- lm(x ~ age)
  MODEL <- summary(MODEL)
  AGE_TREND <- MODEL$coefficients[2,1]
  return(AGE_TREND)
}


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
    mutate(Age = gsub("cor_","",Age)) 
  
  age_trend_real <- df_cor %>% 
    mutate(Age = as.numeric(Age)) %>% 
    mutate(Age = (Age/12)-1) %>% 
    group_by(pair_comparison) %>% 
    summarise(
      CorrTrend = est_age_trend(cor, Age)) %>% 
    ungroup()
  
  age_trend_real_protein_list[[i]] <- age_trend_real
  
}

rm(i,df_test_12,df_test_18,df_test_6,age_trend_real,COR_06,COR_12,COR_18,df_cor,
   names)

# Creating a function to do that permuting mice 1000 times, disrupting only 
# the mouse - age association.

# First, checking how many mice per group

df %>% select(mouse.id, Age, Sex) %>% unique() %>% filter(Age == -0.5) %>% dim()
# [1] 62 3

df %>% select(mouse.id, Age, Sex) %>% unique() %>% filter(Age == 0) %>% dim()
# [1] 62 3

set.seed(2021)

permut_age_trend <- function(df, times) {
  
  age_trend_permut_protein_list <- list()
  
  for (i in unique(df$complex_name)){
    
    age_trend_perm <- list()
    
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
        mutate(Age = gsub("cor_","",Age))
      
      age_trend <- df_cor %>% 
        mutate(Age = as.numeric(Age)) %>% 
        mutate(Age = (Age/12)-1) %>% 
        group_by(pair_comparison) %>% 
        summarise(
          CorrTrend = est_age_trend(cor, Age)) %>% 
        ungroup()
      
      age_trend_perm[[j]] <- age_trend
      
    }
    
    age_trend_permut_protein_list[[i]] <- age_trend_perm
    
    
  }
  
  return(age_trend_permut_protein_list)
}

age_trend_permut_protein_results <- permut_age_trend(df,1000)

# For each complex and each pair, compute the 95% quantile confidence

age_trend_permut_protein_results <- sapply(
  age_trend_permut_protein_results,bind_rows, 
                USE.NAMES = TRUE, simplify = FALSE)

save(age_trend_permut_protein_results,
     file = "Downstream_Analysis/Results/proteincomplex_perms_protein.RData")


# Using bigEmpPVals function from DGCA to compute the pvals from the permutation

empirical_pvals <- list() 

for (i in names(age_trend_permut_protein_results)){
  
  complex_df_perm <- age_trend_permut_protein_results[[i]]
  
  complex_df_real <- age_trend_real_protein_list[[i]]
  
  complex_pair_comparison <- tibble(
    pair_comparison = complex_df_real$pair_comparison,
    CorrTrend = complex_df_real$CorrTrend,
    pval = 0
  )
  
  for (j in unique(complex_df_real$pair_comparison)) {
    
    complex_df_perm_subset <- complex_df_perm %>% 
      filter(pair_comparison == j) %>% 
      select(CorrTrend) %>% 
      unlist()
    
    complex_df_real_subset <- complex_df_real %>% 
      filter(pair_comparison == j) %>% 
      select(CorrTrend) %>% 
      unlist()
    
    pval <- bigEmpPVals(stat = abs(complex_df_real_subset), 
                        stat0 = abs(complex_df_perm_subset))
    
    complex_pair_comparison$pval[
      complex_pair_comparison$pair_comparison == j
      ] <- pval
    
  }
  empirical_pvals[[i]] <- complex_pair_comparison
  
}

rm(i, j, complex_df_perm, complex_df_real, complex_pair_comparison,
   complex_df_perm_subset, complex_df_real_subset, pval)

adjust <- function (x){
  
  padj <- adjustPVals(x$pval, "BH")
  
  x$pval.adj <- padj
  
  return(x)
  
}


complex_protein_age_trend_final <- sapply(empirical_pvals, adjust, 
                                          simplify = FALSE)

complex_protein_age_trend_final <- bind_rows(complex_protein_age_trend_final,
                                             .id = "complex_name")

complex_protein_age_trend_final %>% 
  filter(pval.adj < 0.1) %>% 
  View()

complex_protein_age_trend_final_significant <- complex_protein_age_trend_final %>%
  filter(pval.adj < 0.1)

# From 120 protein complexes we have analyzed, 52 showed significant pairs
# changing their correlation with age (total of 354 correlations changing).
# From those, 325 decrease with age!!

# Saving data for supplemental material

save <- complex_protein_age_trend_final %>% 
  separate(pair_comparison, into = c("gene_symbol","symbol_2"),sep = "_") %>% 
  left_join(df %>% select(protein.id,gene_symbol,complex_name) %>% unique(), 
            by =c("gene_symbol","complex_name")) %>% 
  select(complex_name,gene_symbol,symbol_2,CorrTrend, pval, pval.adj, protein.id) %>% 
  rename(symbol_1 = gene_symbol,
         gene_symbol = symbol_2,
         protein.id_1 = protein.id) %>% 
  left_join(df %>% select(protein.id,gene_symbol,complex_name) %>% unique()
            , by = c("gene_symbol","complex_name")) %>% 
  select(complex_name, symbol_1, protein.id_1, gene_symbol, protein.id, 
         CorrTrend, pval, pval.adj) %>% 
  rename(symbol_2 = gene_symbol) %>% 
  unite("Pair_comparison_ID", protein.id_1, protein.id, sep = "_") %>% 
  unite("Pair_comparison_symbol", symbol_1, symbol_2, sep = "_") %>% 
  rename(Complex_name = complex_name, 
         Age_effect = CorrTrend, 
         p = pval, 
         adjusted_p = pval.adj) %>% 
  select(Complex_name, Pair_comparison_ID, Pair_comparison_symbol, Age_effect, 
         p, adjusted_p)

# Check out which complexes have the largest proportion of pairs changing

save %>% 
  group_by(Complex_name) %>% 
  mutate(Total = length(Pair_comparison_ID)) %>% 
  ungroup() %>% 
  filter(adjusted_p < 0.1) %>% 
  group_by(Complex_name) %>% 
  mutate(Significant = length(Pair_comparison_ID)) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(Prop = (Significant/Total)*100) %>% 
  ungroup() %>% 
  select(Complex_name, Prop) %>% 
  unique() %>% 
  arrange(desc(Prop)) %>% 
  View()

# write.csv(save,
#   file = "Downstream_Analysis/Results/proteincomplex_protein_pairs.csv",
#       row.names = FALSE)

rm(save)

# Plot the heatmaps for the age trend for each pair

age_trend_real_protein <- bind_rows(age_trend_real_protein_list,
                                    .id = "complex_name")

col <- colorRampPalette(c("#2166ac", "#67a9cf", "#d1e5f0", "#f7f7f7", "#fddbc7",
                          "#ef8a62", "#b2182b"))

order_rule <- vector('list', length = length(unique(age_trend_real_protein$complex_name)))
names(order_rule) <- unique(age_trend_real_protein$complex_name)

pdf("Downstream_Analysis/Results/proteincomplex_age_trend_protein.pdf", width = 8)

for (i in unique(age_trend_real_protein$complex_name)){
  
  d <- age_trend_real_protein %>% 
    filter(complex_name == i) %>% 
    mutate(signif = ifelse(pair_comparison %in% 
                             complex_protein_age_trend_final_significant$pair_comparison[
                               complex_protein_age_trend_final_significant$complex_name == i], 
                           "TRUE", "FALSE")) %>% 
    mutate(signif = as.logical(signif)) %>% 
    select(pair_comparison,CorrTrend, signif) %>% 
    separate(pair_comparison, into = c("protein_1","protein_2"),sep = "_")
  
  # Define order
  IDs <- df %>% 
    filter(complex_name == i) %>% 
    group_by(complex_name, gene.id, protein.id, gene_symbol) %>% 
    summarise() %>% 
    ungroup()
  
  TREND_MAT <- matrix(
    data = 0,
    nrow = length(unique(c(d$protein_1, d$protein_2))),
    ncol = length(unique(c(d$protein_1, d$protein_2))),
    dimnames = list(unique(c(d$protein_1, d$protein_2)), 
                    unique(c(d$protein_1, d$protein_2)))
  )
  
  for(r in 1:nrow(d)){
    TREND_MAT[d$protein_1[r], d$protein_2[r]] <- d$CorrTrend[r]
    TREND_MAT[d$protein_2[r], d$protein_1[r]] <- d$CorrTrend[r]
  }
  
  CLUSTERING <- hclust(dist(TREND_MAT), method = 'ward.D2')
  
  ORDER <- colnames(TREND_MAT)[CLUSTERING$order]
  
  order_rule[[i]] <- IDs %>% 
    mutate(
      gene_symbol = factor(gene_symbol, levels = ORDER)
    ) %>% 
    arrange(gene_symbol) %>% 
    mutate(
      gene.id = factor(gene.id, levels = unique(gene.id))
    ) %>% 
    arrange(gene.id)
  
  d <- d %>% 
    mutate(
      protein_1_old = factor(protein_1, 
                             levels = levels(order_rule[[i]]$gene_symbol)),
      protein_2_old = factor(protein_2, 
                             levels = levels(order_rule[[i]]$gene_symbol))
    ) %>% 
    mutate(
      protein_1 = ifelse(
        as.numeric(protein_1_old) < as.numeric(protein_2_old),
        as.character(protein_1_old), as.character(protein_2_old)
      ),
      protein_2 = ifelse(
        as.numeric(protein_1_old) < as.numeric(protein_2_old),
        as.character(protein_2_old), as.character(protein_1_old)
      )
    ) %>% 
    mutate(
      protein_1 = factor(protein_1, 
                         levels = levels(order_rule[[i]]$gene_symbol)),
      protein_2 = factor(protein_2, 
                         levels = levels(order_rule[[i]]$gene_symbol))
    )
  
  SIZE <- 32/length(unique(d$protein_1))
  
  #quartz()
  plot <- d %>% 
    ggplot(aes(protein_1, protein_2)) +
    geom_hline(yintercept = seq(2, length(levels(d$protein_2))) - .5, 
               color="#eeeeee", size = .25) + 
    geom_vline(xintercept = seq(2, length(levels(d$protein_1))) - .5, 
               color="#eeeeee", size = .25) +
    ggtitle(paste0('Age trend estimates - ',i)) +
    theme_bw() +
    xlab('Gene symbol') +
    ylab('Gene symbol') +
    geom_tile(aes(fill = CorrTrend), color='#eeeeee') +
    geom_point(data = d %>% filter(signif), color='black', shape = 8, 
               size = SIZE, stroke = SIZE/2) +
    scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b" , 
                         name = "Age estimate") +
    scale_x_discrete(
      drop = FALSE, 
      limits = levels(d$protein_1)[-length(levels(d$protein_1))]
    ) +
    scale_y_discrete(
      drop = FALSE, 
      limits = levels(d$protein_2)[-1]
    ) +
    theme(
      axis.text.x=element_text(angle=90, size = 10, colour = "black"),
      axis.text.y=element_text(size = 10, colour = 'black'),
      axis.ticks=element_blank(),
      axis.line=element_blank(),
      panel.border=element_blank(),
      panel.grid.major = element_blank()
    )
  
  print(plot)
  
}

dev.off()

rm(age_trend_real_protein_list, d, empirical_pvals, IDs, CLUSTERING, plot, 
   TREND_MAT, i, ORDER, r,SIZE)

################################################################################

# 1) Transcripts

age_trend_real_gene_list <- list()

for (i in unique(df$complex_name)){
  
  df_test_6 <- df %>% 
    filter(complex_name == i,
           Age == -0.5) %>% 
    select(mouse.id,gene_symbol,gene_expression) %>% 
    unique() %>% 
    spread(gene_symbol,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_6 <- df_test_6[,!is.na(colSums(df_test_6))]

  
  df_test_12 <- df %>% 
    filter(complex_name == i,
           Age == 0) %>% 
    select(mouse.id, gene_symbol,gene_expression) %>% 
    unique() %>% 
    spread(gene_symbol,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_12 <- df_test_12[,!is.na(colSums(df_test_12))]
  
  df_test_18 <- df %>% 
    filter(complex_name == i,
           Age == 0.5) %>% 
    select(mouse.id, gene_symbol,gene_expression) %>% 
    unique() %>% 
    spread(gene_symbol,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_18 <- df_test_18[,!is.na(colSums(df_test_18))]
  
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
    mutate(Age = gsub("cor_","",Age)) 
  
  age_trend_real <- df_cor %>% 
    mutate(Age = as.numeric(Age)) %>% 
    mutate(Age = (Age/12)-1) %>% 
    group_by(pair_comparison) %>% 
    summarise(
      CorrTrend = est_age_trend(cor, Age)) %>% 
    ungroup()
  
  age_trend_real_gene_list[[i]] <- age_trend_real
  
}

rm(i,df_test_12,df_test_18,df_test_6,age_trend_real,COR_06,COR_12,COR_18,df_cor,
   names)


# Creating a function to do that permuting mice 1000 times, disrupting only 
#the mouse - age association.

set.seed(2022)

permut_age_trend <- function(df, times) {
  
  age_trend_permut_gene_list <- list()
  
  for (i in unique(df$complex_name)){
    
    age_trend_perm <- list()
    
    for (j in 1:times) {
      
      df_test_6 <- df %>% 
        filter(complex_name == i,
               mouse.id %in% sample(unique(df$mouse.id),62)) %>% 
        select(mouse.id,gene_symbol,gene_expression) %>% 
        unique() %>% 
        spread(gene_symbol,gene_expression) %>% 
        column_to_rownames("mouse.id")
      
      df_test_6 <- df_test_6[,!is.na(colSums(df_test_6))]
      
      df_test_12 <- df %>% 
        filter(complex_name == i,
               mouse.id %in% sample(unique(df$mouse.id[
                 !df$mouse.id %in% rownames(df_test_6)
                 ]),62)) %>% 
        select(mouse.id, gene_symbol,gene_expression) %>% 
        unique() %>% 
        spread(gene_symbol,gene_expression) %>% 
        column_to_rownames("mouse.id")
      
      df_test_12 <- df_test_12[,!is.na(colSums(df_test_12))]
      
      df_test_18 <- df %>% 
        filter(complex_name == i,
               !mouse.id %in% rownames(df_test_6) &
                 !mouse.id %in% rownames(df_test_12)) %>% 
        select(mouse.id, gene_symbol,gene_expression) %>% 
        unique() %>% 
        spread(gene_symbol,gene_expression) %>% 
        column_to_rownames("mouse.id")
      
      df_test_18 <- df_test_18[,!is.na(colSums(df_test_18))]
      
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
        mutate(Age = gsub("cor_","",Age)) 
      
      age_trend <- df_cor %>% 
        mutate(Age = as.numeric(Age)) %>% 
        mutate(Age = (Age/12)-1) %>% 
        group_by(pair_comparison) %>% 
        summarise(
          CorrTrend = est_age_trend(cor, Age)) %>% 
        ungroup()
      
      age_trend_perm[[j]] <- age_trend
      
    }
    
    age_trend_permut_gene_list[[i]] <- age_trend_perm
    
    
  }
  
  return(age_trend_permut_gene_list)
}

age_trend_permut_gene_results <- permut_age_trend(df,1000)

age_trend_permut_gene_results <- sapply(age_trend_permut_gene_results,bind_rows, 
                                        USE.NAMES = TRUE, simplify = FALSE)

save(age_trend_permut_gene_results, 
     file = "Downstream_Analysis/Results/proteincomplex_perms_transcript.RData")


# Using bigEmpPVals function from DGCA to compute the pvals from the permutation

empirical_pvals <- list() 

for (i in names(age_trend_permut_gene_results)){
  
  complex_df_perm <- age_trend_permut_gene_results[[i]]
  
  complex_df_real <- age_trend_real_gene_list[[i]]
  
  complex_pair_comparison <- tibble(
    pair_comparison = complex_df_real$pair_comparison,
    CorrTrend = complex_df_real$CorrTrend,
    pval = 0
  )
  
  for (j in unique(complex_df_real$pair_comparison)) {
    
    complex_df_perm_subset <- complex_df_perm %>% 
      filter(pair_comparison == j) %>% 
      select(CorrTrend) %>% 
      unlist()
    
    complex_df_real_subset <- complex_df_real %>% 
      filter(pair_comparison == j) %>% 
      select(CorrTrend) %>% 
      unlist()
    
    pval <- bigEmpPVals(stat = abs(complex_df_real_subset), 
                        stat0 = abs(complex_df_perm_subset))
    
    complex_pair_comparison$pval[
      complex_pair_comparison$pair_comparison == j
      ] <- pval
    
  }
  empirical_pvals[[i]] <- complex_pair_comparison
  
}

rm(i, j, complex_df_perm, complex_df_real, complex_pair_comparison,
   complex_df_perm_subset, complex_df_real_subset, pval)


complex_gene_age_trend_final <- sapply(empirical_pvals, adjust, simplify = FALSE)

complex_gene_age_trend_final <- bind_rows(complex_gene_age_trend_final,
                                          .id = "complex_name")

complex_gene_age_trend_final %>% 
  filter(pval.adj < 0.1) %>% 
  View()

# From 120 protein complexes we have analyzed, 20 showed significant pairs
# changing their correlation with age (total of 52 correlations changing).
# From those, 46 increase with age!!

complex_gene_age_trend_final_significant <- complex_gene_age_trend_final %>% 
  filter(pval.adj < 0.1) 

# Saving data for suplemental material

save <- complex_gene_age_trend_final %>% 
  separate(pair_comparison, into = c("gene_symbol","symbol_2"),sep = "_") %>% 
  left_join(df %>% select(gene.id,gene_symbol,complex_name) %>% unique(), 
            by =c("gene_symbol","complex_name")) %>% 
  select(complex_name,gene_symbol,symbol_2,CorrTrend, pval, pval.adj, gene.id) %>% 
  rename(symbol_1 = gene_symbol,
         gene_symbol = symbol_2,
         gene.id_1 = gene.id) %>% 
  left_join(df %>% select(gene.id,gene_symbol,complex_name) %>% unique()
            , by = c("gene_symbol","complex_name")) %>% 
  select(complex_name, symbol_1, gene.id_1, gene_symbol, gene.id, 
         CorrTrend, pval, pval.adj) %>% 
  rename(symbol_2 = gene_symbol) %>% 
  unite("Pair_comparison_ID", gene.id_1, gene.id, sep = "_") %>% 
  unite("Pair_comparison_symbol", symbol_1, symbol_2, sep = "_") %>% 
  rename(Complex_name = complex_name, 
         Age_effect = CorrTrend, 
         p = pval, 
         adjusted_p = pval.adj) %>% 
  select(Complex_name, Pair_comparison_ID, Pair_comparison_symbol, Age_effect, 
         p, adjusted_p)

# write.csv(save, 
# file = "Downstream_Analysis/Results/proteincomplex_transcripts_pairs.csv",
# row.names = FALSE)

# Check out which complexes have the largest proportion of pairs changing

save %>% 
  group_by(Complex_name) %>% 
  mutate(Total = length(Pair_comparison_ID)) %>% 
  ungroup() %>% 
  filter(adjusted_p < 0.1) %>% 
  group_by(Complex_name) %>% 
  mutate(Significant = length(Pair_comparison_ID)) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(Prop = (Significant/Total)*100) %>% 
  ungroup() %>% 
  select(Complex_name, Prop) %>% 
  unique() %>% 
  arrange(desc(Prop)) %>% 
  View()

rm(save)

# Plot the heatmaps for the age trend for each pair

age_trend_real_gene <- bind_rows(age_trend_real_gene_list,
                                 .id = "complex_name")

pdf("Downstream_Analysis/Results/proteincomplex_age_trend_transcript.pdf",
    width = 8)

for (i in unique(age_trend_real_gene$complex_name)){
  
  d <- age_trend_real_gene %>% 
    filter(complex_name == i) %>% 
    mutate(signif = ifelse(pair_comparison %in% 
                             complex_gene_age_trend_final_significant$pair_comparison[
                               complex_gene_age_trend_final_significant$complex_name == i], 
                           "TRUE", "FALSE")) %>% 
    mutate(signif = as.logical(signif)) %>% 
    select(pair_comparison,CorrTrend,signif) %>% 
    separate(pair_comparison, into = c("gene_1","gene_2"),sep = "_")
  
  d <- d %>% 
    mutate(
      gene_1_old = factor(gene_1, levels = levels(order_rule[[i]]$gene_symbol)),
      gene_2_old = factor(gene_2, levels = levels(order_rule[[i]]$gene_symbol))
    ) %>% 
    mutate(
      gene_1 = ifelse(
        as.numeric(gene_1_old) < as.numeric(gene_2_old),
        as.character(gene_1_old), as.character(gene_2_old)
      ),
      gene_2 = ifelse(
        as.numeric(gene_1_old) < as.numeric(gene_2_old),
        as.character(gene_2_old), as.character(gene_1_old)
      )
    ) %>% 
    mutate(
      gene_1 = factor(gene_1, levels = levels(order_rule[[i]]$gene_symbol)),
      gene_2 = factor(gene_2, levels = levels(order_rule[[i]]$gene_symbol))
    )
  
  SIZE <- 32/length(unique(d$gene_1))
  
  #quartz()
  plot <- d %>% 
    ggplot(aes(gene_1, gene_2)) +
    geom_hline(yintercept = seq(2, length(levels(d$gene_2))) - .5, 
               color="#eeeeee", size = .25) + 
    geom_vline(xintercept = seq(2, length(levels(d$gene_1))) - .5, 
               color="#eeeeee", size = .25) +
    ggtitle(paste0('Age trend estimates - ',i)) +
    theme_bw() +
    xlab('Gene symbol') +
    ylab('Gene symbol') +
    geom_tile(aes(fill = CorrTrend), color='#eeeeee') +
    geom_point(data = d %>% filter(signif), color='black', shape = 8, 
               size = SIZE, stroke = SIZE/2) +
    scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b" , 
                         name = "Age estimate") +
    scale_x_discrete(
      drop = FALSE, 
      limits = levels(d$gene_1)[-length(levels(d$gene_1))]
    ) +
    scale_y_discrete(
      drop = FALSE, 
      limits = levels(d$gene_2)[-1]
    ) +
    theme(
      axis.text.x=element_text(angle=90, size = 10, colour = "black"),
      axis.text.y=element_text(size = 10, colour = 'black'),
      axis.ticks=element_blank(),
      axis.line=element_blank(),
      panel.border=element_blank(),
      panel.grid.major = element_blank()
    )
  
  print(plot)
  
}

dev.off()

rm(age_trend_real_gene_list, d, empirical_pvals, i, SIZE,plot)


test_gene <- complex_gene_age_trend_final_significant %>% 
  separate(pair_comparison, 
           into = c("gene_1","gene_2"),
           sep = "_") %>% 
  select(complex_name, gene_1, gene_2, CorrTrend)

test_protein <- complex_protein_age_trend_final_significant %>% 
  separate(pair_comparison, 
           into = c("protein_1","protein_2"),
           sep = "_") %>% 
  select(complex_name, protein_1, protein_2, CorrTrend)

test_gene %>% 
  inner_join(test_protein, by = "complex_name") %>% 
  filter(gene_1 == protein_1 & gene_2 == protein_2 | 
           gene_1 == protein_2 & gene_2 == protein_1) 

# There are 5 pair comparisons in which RNA cor increase with age and protein
# cor decrease with age: 

# A tibble: 5 x 7
# complex_name                           gene_1 gene_2     CorrTrend.x   protein_1 protein_2     CorrTrend.y
# <chr>                                 <chr>   <chr>        <dbl>         <chr>     <chr>         <dbl>
#   1 COP9_signalosome                    Cops3  Cops8        0.501       Cops3     Cops8          -0.508
# 2 Mitochondrial_complex_V               Atp5h  Atp5pb       0.255       Atp5h     Atp5pb         -0.552
# 3 Mitochondrial_complex_IV              Cox5a  Cox7c        0.299       Cox5a     Cox7c          -0.628
# 4 multi-eIF_complex                     Eif3b  Eif3i        0.811       Eif3b     Eif3i          -0.307
# 5 SMG-1-Upf1-eRF1-eRF3_complex_(SURF)   Etf1   Gspt2        0.485       Etf1      Gspt2          -0.423

rm(test_gene, order_rule, test_protein)

# Now, plotting only the proteasome for publication and ordering by 
# categories. 

pdf("Downstream_Analysis/Results/proteincomplex_age_trend_proteasome.pdf",
    width = 8, height = 14)
  
  d <- age_trend_real_gene %>% 
    filter(complex_name == "26S_Proteasome") %>% 
    mutate(signif = ifelse(pair_comparison %in% 
                                    complex_gene_age_trend_final_significant$pair_comparison[
                                      complex_gene_age_trend_final_significant$complex_name == "26S_Proteasome"],
                           "TRUE", "FALSE")) %>% 
    mutate(signif = as.logical(signif)) %>% 
    select(pair_comparison,CorrTrend,signif) %>% 
    separate(pair_comparison, into = c("gene_1","gene_2"),sep = "_")
  
  
  ORDER <- unique(
               c(
                 unique(d$gene_1),
                   unique(d$gene_2)))
  
  ORDER <- factor(ORDER, levels = rev(
    c(
    paste0("Psma",seq(1:7)),
    paste0("Psmb",seq(1:10)),
    paste0("Psmc",seq(1:6)),
    paste0("Psmd",seq(1:4)),
    "Psmd4.1",
    paste0("Psmd",c(5:9)),
    paste0("Psmd",c(11:14)),
    paste0("Psme",seq(1:4)),
    "Psmf1","Psmg1"
  )))
  
  ORDER <- as.character(
    sort(
      ORDER)
  )
  
  
  d <- d %>% 
    mutate(
      gene_1_old = factor(gene_1, 
                          levels = ORDER),
      gene_2_old = factor(gene_2, levels = ORDER)
    ) %>% 
    mutate(
      gene_1 = ifelse(
        as.numeric(gene_1_old) < as.numeric(gene_2_old),
        as.character(gene_1_old), as.character(gene_2_old)
      ),
      gene_2 = ifelse(
        as.numeric(gene_1_old) < as.numeric(gene_2_old),
        as.character(gene_2_old), as.character(gene_1_old)
      )
    ) %>% 
    mutate(
      gene_1 = factor(gene_1, levels = rev(ORDER)),
      gene_2 = factor(gene_2, levels = ORDER)
    )
  
  SIZE <- 32/length(unique(d$gene_1))
  
  #quartz()
  plot1 <- d %>% 
    ggplot(aes(gene_1, gene_2)) +
    geom_hline(yintercept = seq(2, length(levels(d$gene_2))) - .5, 
               color="#eeeeee", size = .25) + 
    geom_vline(xintercept = seq(2, length(levels(d$gene_1))) - .5, 
               color="#eeeeee", size = .25) +
    ggtitle("Transcript level") +
    theme_bw() +
    xlab('Gene symbol') +
    ylab('Gene symbol') +
    geom_tile(aes(fill = CorrTrend), color='#eeeeee') +
    geom_point(data = d %>% filter(signif), color='black', shape = 8, 
               size = SIZE, stroke = SIZE/2) +
    scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b" , 
                         name = "Age effect") +
    scale_x_discrete(
      drop = FALSE, 
      limits = levels(d$gene_1)[-1], position = "top"
    ) +
    scale_y_discrete(
      drop = FALSE, 
      limits = levels(d$gene_2)[-1]
    ) +
    theme(
      axis.text.x=element_text(angle=90, size = 10, colour = "black"),
      axis.text.y=element_text(size = 10, colour = 'black'),
      axis.ticks=element_blank(),
      axis.line=element_blank(),
      panel.border=element_blank(),
      panel.grid.major = element_blank()
    )

d <- age_trend_real_protein %>% 
  filter(complex_name == "26S_Proteasome") %>% 
  mutate(signif = ifelse(pair_comparison %in% 
                           complex_protein_age_trend_final_significant$pair_comparison[
                             complex_protein_age_trend_final_significant$complex_name == "26S_Proteasome"],
                         "TRUE", "FALSE")) %>% 
  mutate(signif = as.logical(signif)) %>% 
  select(pair_comparison,CorrTrend, signif) %>% 
  separate(pair_comparison, into = c("protein_1","protein_2"),sep = "_")


d <- d %>% 
  mutate(
    protein_1_old = factor(protein_1, 
                           levels = ORDER),
    protein_2_old = factor(protein_2, 
                           levels = ORDER)
  ) %>% 
  mutate(
    protein_1 = ifelse(
      as.numeric(protein_1_old) < as.numeric(protein_2_old),
      as.character(protein_1_old), as.character(protein_2_old)
    ),
    protein_2 = ifelse(
      as.numeric(protein_1_old) < as.numeric(protein_2_old),
      as.character(protein_2_old), as.character(protein_1_old)
    )
  ) %>% 
  mutate(
    protein_1 = factor(protein_1, 
                       levels = rev(ORDER)),
    protein_2 = factor(protein_2, 
                       levels = ORDER)
  )

SIZE <- 32/length(unique(d$protein_1))

#quartz()
plot2 <- d %>% 
  ggplot(aes(protein_1, protein_2)) +
  geom_hline(yintercept = seq(2, length(levels(d$protein_2))) - .5, 
             color="#eeeeee", size = .25) + 
  geom_vline(xintercept = seq(2, length(levels(d$protein_1))) - .5, 
             color="#eeeeee", size = .25) +
  ggtitle("Protein level") +
  theme_bw() +
  xlab('Gene symbol') +
  ylab('Gene symbol') +
  geom_tile(aes(fill = CorrTrend), color='#eeeeee') +
  geom_point(data = d %>% filter(signif), color='black', shape = 8, 
             size = SIZE, stroke = SIZE/2) +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b" , 
                       name = "Age effect") +
  scale_x_discrete(
    drop = FALSE, 
    limits = levels(d$protein_1)[-1], position = "top"
  ) +
  scale_y_discrete(
    drop = FALSE, 
    limits = levels(d$protein_2)[-1]
  ) +
  theme(
    axis.text.x=element_text(angle=90, size = 10, colour = "black"),
    axis.text.y=element_text(size = 10, colour = 'black'),
    axis.ticks=element_blank(),
    axis.line=element_blank(),
    panel.border=element_blank(),
    panel.grid.major = element_blank()
  )

cowplot::plot_grid(plotlist = list(plot1,plot2),ncol = 1)

dev.off()

rm(d, plot1, plot2, ORDER, SIZE)



# Now, plotting supplemental figures for Large_Drosha_complex, Nuclear_pore_complex_(NPC),
# cytoplasmic_ribosomal_large_subunit, chaperonin-containing_T-complex, 
# Mitochondrial_complex_I, Mitochondrial_complex_II, Mitochondrial_complex_III, 
# Mitochondrial_complex_IV, Mitochondrial_complex_V

pdf("Downstream_Analysis/Results/proteincomplex_supplemental_age_effects.pdf",
    width = 18, height = 8)

for (i in c(
           'chaperonin-containing_T-complex', 
           'Suv39h1',
           'dynactin_complex',
           'multi-eIF_complex',
           'Large_Drosha_complex', 
           'Ubiquilin-proteasome_complex',
           'TNF-alpha/NF-kappa_B_signaling_complex_7',
           'SMG-1-Upf1-eRF1-eRF3_complex_(SURF)',
           'COP9_signalosome',
            'cytoplasmic_ribosomal_large_subunit',
            'Nuclear_pore_complex_(NPC)',
            'cytoplasmic_ribosomal_large_subunit', 
            'Mitochondrial_complex_I', 
            'Mitochondrial_complex_II', 
            'Mitochondrial_complex_III', 
            'Mitochondrial_complex_IV', 
            'Mitochondrial_complex_V')) {

  d <- age_trend_real_gene %>% 
    filter(complex_name == i) %>% 
    mutate(signif = ifelse(pair_comparison %in% 
                             complex_gene_age_trend_final_significant$pair_comparison[
                               complex_gene_age_trend_final_significant$complex_name == i],
                           "TRUE", "FALSE")) %>% 
    mutate(signif = as.logical(signif)) %>% 
    select(pair_comparison,CorrTrend,signif) %>% 
    separate(pair_comparison, into = c("gene_1","gene_2"),sep = "_")
  
  d <- d %>% 
    mutate(
      gene_1_old = factor(gene_1, levels = levels(order_rule[[i]]$gene_symbol)),
      gene_2_old = factor(gene_2, levels = levels(order_rule[[i]]$gene_symbol))
    ) %>% 
    mutate(
      gene_1 = ifelse(
        as.numeric(gene_1_old) < as.numeric(gene_2_old),
        as.character(gene_1_old), as.character(gene_2_old)
      ),
      gene_2 = ifelse(
        as.numeric(gene_1_old) < as.numeric(gene_2_old),
        as.character(gene_2_old), as.character(gene_1_old)
      )
    ) %>% 
    mutate(
      gene_1 = factor(gene_1, levels = levels(order_rule[[i]]$gene_symbol)),
      gene_2 = factor(gene_2, levels = levels(order_rule[[i]]$gene_symbol))
    )
  
  SIZE <- 32/length(unique(d$gene_1))
  
  #quartz()
  plot1 <- d %>% 
    ggplot(aes(gene_1, gene_2)) +
    geom_hline(yintercept = seq(2, length(levels(d$gene_2))) - .5, 
               color="#eeeeee", size = .25) + 
    geom_vline(xintercept = seq(2, length(levels(d$gene_1))) - .5, 
               color="#eeeeee", size = .25) +
    ggtitle("Transcript") +
    theme_bw() +
    theme(title = element_text(size = 15)) +
    xlab('') +
    ylab('') +
    geom_tile(aes(fill = CorrTrend), color='#eeeeee') +
    geom_point(data = d %>% filter(signif), color='black', shape = 8, 
               size = SIZE, stroke = SIZE/2) +
    scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b" , 
                         name = "Age effect") +
    scale_x_discrete(
      drop = FALSE, 
      limits = levels(d$gene_1)[-length(levels(d$gene_1))]
    ) +
    scale_y_discrete(
      drop = FALSE, 
      limits = levels(d$gene_2)[-1]
    ) +
    theme(
      axis.text.x=element_text(angle=90, size = 10, colour = "black"),
      axis.text.y=element_text(size = 10, colour = 'black'),
      axis.ticks=element_blank(),
      axis.line=element_blank(),
      panel.border=element_blank(),
      panel.grid.major = element_blank()
    )
  
  
  d <- age_trend_real_protein %>% 
    filter(complex_name == i) %>% 
    mutate(signif = ifelse(pair_comparison %in% 
                             complex_protein_age_trend_final_significant$pair_comparison[
                               complex_protein_age_trend_final_significant$complex_name == i],
                           "TRUE", "FALSE")) %>% 
    mutate(signif = as.logical(signif)) %>% 
    select(pair_comparison,CorrTrend, signif) %>% 
    separate(pair_comparison, into = c("protein_1","protein_2"),sep = "_")
  
  
  d <- d %>% 
    mutate(
      protein_1_old = factor(protein_1, 
                             levels = levels(order_rule[[i]]$gene_symbol)),
      protein_2_old = factor(protein_2, 
                             levels = levels(order_rule[[i]]$gene_symbol))
    ) %>% 
    mutate(
      protein_1 = ifelse(
        as.numeric(protein_1_old) < as.numeric(protein_2_old),
        as.character(protein_1_old), as.character(protein_2_old)
      ),
      protein_2 = ifelse(
        as.numeric(protein_1_old) < as.numeric(protein_2_old),
        as.character(protein_2_old), as.character(protein_1_old)
      )
    ) %>% 
    mutate(
      protein_1 = factor(protein_1, 
                         levels = levels(order_rule[[i]]$gene_symbol)),
      protein_2 = factor(protein_2, 
                         levels = levels(order_rule[[i]]$gene_symbol))
    )
  
  SIZE <- 32/length(unique(d$protein_1))
  
  #quartz()
  plot2 <- d %>% 
    ggplot(aes(protein_1, protein_2)) +
    geom_hline(yintercept = seq(2, length(levels(d$protein_2))) - .5, 
               color="#eeeeee", size = .25) + 
    geom_vline(xintercept = seq(2, length(levels(d$protein_1))) - .5, 
               color="#eeeeee", size = .25) +
    ggtitle("Protein") +
    theme_bw() +
    theme(title = element_text(size = 15)) +
    xlab('') +
    ylab('') +
    geom_tile(aes(fill = CorrTrend), color='#eeeeee') +
    geom_point(data = d %>% filter(signif), color='black', shape = 8, 
               size = SIZE, stroke = SIZE/2) +
    scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b" , 
                         name = "Age effect") +
    scale_x_discrete(
      drop = FALSE, 
      limits = levels(d$protein_1)[-length(levels(d$protein_1))]
    ) +
    scale_y_discrete(
      drop = FALSE, 
      limits = levels(d$protein_2)[-1]
    ) +
    theme(
      axis.text.x=element_text(angle=90, size = 10, colour = "black"),
      axis.text.y=element_text(size = 10, colour = 'black'),
      axis.ticks=element_blank(),
      axis.line=element_blank(),
      panel.border=element_blank(),
      panel.grid.major = element_blank()
    )
  
  print(cowplot::plot_grid(plot1, plot2, 
                           ncol = 2, 
                           labels = i,
                           label_size = 20,
                           vjust = 1.3,
                           scale = 0.9))
  
}

dev.off()
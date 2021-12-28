################## Protein complexes analysis - correlation ####################

# KIDNEY

# This script uses the information on Greg's table 
# to gather the genes and proteins of each complex
# and check how their correlation is affected by age

# Author: Isabela Gerdes Gyuricza
# Date: 07_26_2021

################################################################################
############ clear workspace

rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

################################################################################
############ libraries

library(corrplot)
library(tidyverse)

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

genes_info <- ensimplR::batchGenes(c("ENSMUSG00000021577","ENSMUSG00000009863","ENSMUSG00000058076","ENSMUSG00000000171"))


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
              select(Mouse.ID, Age, Sex), by = "Mouse.ID") 

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

######## Testing change in correlation among proteins from the same complex and 
######## how they change with age 

col <- colorRampPalette(c("#2166ac", "#67a9cf", "#d1e5f0", "#f7f7f7", "#fddbc7",
                          "#ef8a62", "#b2182b"))

pdf("Downstream_Analysis/Results/KIDNEY_proteincomplex_corplots_protein.pdf",
    width = 14, height = 5.5)

order_rule <- vector('list', length = length(unique(df$complex_name)))
names(order_rule) <- unique(df$complex_name)

for (i in unique(df$complex_name)){
  
  df_test_6 <- df %>% 
    filter(complex_name == i,
           Age == 6) %>% 
    select(mouse.id,gene_symbol,protein_expression) %>% 
    spread(gene_symbol,protein_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_12 <- df %>% 
    filter(complex_name == i,
           Age == 12) %>% 
    select(mouse.id, gene_symbol,protein_expression) %>% 
    spread(gene_symbol,protein_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_18 <- df %>% 
    filter(complex_name == i,
           Age == 18) %>% 
    select(mouse.id, gene_symbol,protein_expression) %>% 
    spread(gene_symbol,protein_expression) %>% 
    column_to_rownames("mouse.id")
  
  par(mfrow=c(1,3))
  
  COR_06 <- cor(df_test_6, use = "pairwise.complete.obs")
  COR_12 <- cor(df_test_12, use = "pairwise.complete.obs")
  COR_18 <- cor(df_test_18, use = "pairwise.complete.obs")
  
  CLUST <- hclust(dist(abs(COR_06)), method = 'ward.D2')
  ORDER <- CLUST$order
  
  order_rule[[i]] <- CLUST$labels[ORDER]
  
  corrplot(COR_06[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=1,tl.col = "black",col = col(200),
           title = paste0(i, " - Age 6"),mar=c(0,0,1.5,0))
  
  corrplot(COR_12[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=1,tl.col = "black",col = col(200),
           title = paste0(i, " - Age 12"),mar=c(0,0,1.5,0))
  
  corrplot(COR_18[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=1,tl.col = "black",col = col(200),
           title = paste0(i, " - Age 18"),mar=c(0,0,1.5,0))
}

dev.off()

rm(df_test_6, df_test_12, df_test_18,i, COR_06, COR_12, COR_18, CLUST, ORDER)


pdf("Downstream_Analysis/Results/KIDNEY_proteincomplex_corplots_transcript.pdf",
    width = 14, height = 5.5)

for (i in unique(df$complex_name)){
  
  df_test_6 <- df %>% 
    filter(complex_name == i,
           Age == 6) %>% 
    select(mouse.id,gene_symbol,gene_expression) %>% 
    unique() %>% 
    spread(gene_symbol,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_6 <- df_test_6[,!is.na(colSums(df_test_6))]
  
  df_test_12 <- df %>% 
    filter(complex_name == i,
           Age == 12) %>% 
    select(mouse.id, gene_symbol,gene_expression) %>% 
    unique() %>% 
    spread(gene_symbol,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_12 <- df_test_12[,!is.na(colSums(df_test_12))]
  
  df_test_18 <- df %>% 
    filter(complex_name == i,
           Age == 18) %>% 
    select(mouse.id, gene_symbol,gene_expression) %>% 
    unique() %>% 
    spread(gene_symbol,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_18 <- df_test_18[,!is.na(colSums(df_test_18))]
  
  par(mfrow=c(1,3))
  #quartz()
  
  COR_06 <- cor(df_test_6, use = "pairwise.complete.obs")
  COR_12 <- cor(df_test_12, use = "pairwise.complete.obs")
  COR_18 <- cor(df_test_18, use = "pairwise.complete.obs")
  
  CLUST <- hclust(dist(abs(COR_06)), method = 'ward.D2')
  
  ORDER <- CLUST$order
  
  corrplot(COR_06[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=1,tl.col = "black",col = col(200),
           title = paste0(i, " - Age 6"),mar=c(0,0,1.5,0))
  
  corrplot(COR_12[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=1,tl.col = "black",col = col(200),
           title = paste0(i, " - Age 12"),mar=c(0,0,1.5,0))
  
  corrplot(COR_18[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=1,tl.col = "black",col = col(200),
           title = paste0(i, " - Age 18"),mar=c(0,0,1.5,0))
  
}

dev.off()

rm(df_test_6, df_test_12, df_test_18, COR_06, COR_12, COR_18, CLUST, ORDER,i)

######## Now checking the changes in correlation by sex

pdf("Downstream_Analysis/Results/KIDNEY_proteincomplex_corplots_protein_BYSEX.pdf",
    width = 14, height = 10.5)

for (i in unique(df$complex_name)){
  
  df_test_6_M <- df %>% 
    filter(complex_name == i,
           Age == 6,
           Sex == "M") %>% 
    select(mouse.id,gene_symbol,protein_expression) %>% 
    spread(gene_symbol,protein_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_12_M <- df %>% 
    filter(complex_name == i,
           Age == 12,
           Sex == "M") %>% 
    select(mouse.id, gene_symbol,protein_expression) %>% 
    spread(gene_symbol,protein_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_18_M <- df %>% 
    filter(complex_name == i,
           Age == 18,
           Sex == "M") %>% 
    select(mouse.id, gene_symbol,protein_expression) %>% 
    spread(gene_symbol,protein_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_6_F <- df %>% 
    filter(complex_name == i,
           Age == 6,
           Sex == "F") %>% 
    select(mouse.id,gene_symbol,protein_expression) %>% 
    spread(gene_symbol,protein_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_12_F <- df %>% 
    filter(complex_name == i,
           Age == 12,
           Sex == "F") %>% 
    select(mouse.id, gene_symbol,protein_expression) %>% 
    spread(gene_symbol,protein_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_18_F <- df %>% 
    filter(complex_name == i,
           Age == 18,
           Sex == "F") %>% 
    select(mouse.id, gene_symbol,protein_expression) %>% 
    spread(gene_symbol,protein_expression) %>% 
    column_to_rownames("mouse.id")
  
  par(mfrow=c(2,3))
  
  COR_06_M <- cor(df_test_6_M, use = "pairwise.complete.obs")
  COR_12_M <- cor(df_test_12_M, use = "pairwise.complete.obs")
  COR_18_M <- cor(df_test_18_M, use = "pairwise.complete.obs")
  
  COR_06_F <- cor(df_test_6_F, use = "pairwise.complete.obs")
  COR_12_F <- cor(df_test_12_F, use = "pairwise.complete.obs")
  COR_18_F <- cor(df_test_18_F, use = "pairwise.complete.obs")
  
  CLUST <- hclust(dist(abs(COR_06_M)), method = 'ward.D2')
  ORDER <- CLUST$order
  
  order_rule[[i]] <- CLUST$labels[ORDER]
  
  corrplot(COR_06_M[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=1,tl.col = "black",col = col(200),
           title = paste0(i, " - Age 6 - Male"),mar=c(0,0,1.5,0))
  
  corrplot(COR_12_M[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=1,tl.col = "black",col = col(200),
           title = paste0(i, " - Age 12 - Male"),mar=c(0,0,1.5,0))
  
  corrplot(COR_18_M[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=1,tl.col = "black",col = col(200),
           title = paste0(i, " - Age 18 - Male"),mar=c(0,0,1.5,0))
  
  
  
  corrplot(COR_06_F[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=1,tl.col = "black",col = col(200),
           title = paste0(i, " - Age 6 - Female"),mar=c(0,0,1.5,0))
  
  corrplot(COR_12_F[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=1,tl.col = "black",col = col(200),
           title = paste0(i, " - Age 12 - Female"),mar=c(0,0,1.5,0))
  
  corrplot(COR_18_F[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=1,tl.col = "black",col = col(200),
           title = paste0(i, " - Age 18 - Female"),mar=c(0,0,1.5,0))
}

dev.off()

rm(df_test_6_M, df_test_12_M, df_test_18_M,
   df_test_6_F, df_test_12_F, df_test_18_F,
   COR_06_M, COR_12_M, COR_18_M,
   COR_06_F, COR_12_F, COR_18_F, CLUST, ORDER)


pdf("Downstream_Analysis/Results/KIDNEY_proteincomplex_corplots_transcript_BYSEX.pdf",
    width = 14, height = 10.5)

order_rule <- vector('list', length = length(unique(df$complex_name)))
names(order_rule) <- unique(df$complex_name)

for (i in unique(df$complex_name)){
  
  df_test_6_M <- df %>% 
    filter(complex_name == i,
           Age == 6,
           Sex == "M") %>% 
    select(mouse.id,gene_symbol,gene_expression) %>% 
    unique() %>% 
    spread(gene_symbol,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_6_M <- df_test_6_M[,!is.na(colSums(df_test_6_M))]
  
  df_test_12_M <- df %>% 
    filter(complex_name == i,
           Age == 12,
           Sex == "M") %>% 
    select(mouse.id, gene_symbol,gene_expression) %>% 
    unique() %>% 
    spread(gene_symbol,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_12_M <- df_test_12_M[,!is.na(colSums(df_test_12_M))]
  
  
  df_test_18_M <- df %>% 
    filter(complex_name == i,
           Age == 18,
           Sex == "M") %>% 
    select(mouse.id, gene_symbol,gene_expression) %>% 
    unique() %>% 
    spread(gene_symbol,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_18_M <- df_test_18_M[,!is.na(colSums(df_test_18_M))]
  
  
  df_test_6_F <- df %>% 
    filter(complex_name == i,
           Age == 6,
           Sex == "F") %>% 
    select(mouse.id,gene_symbol,gene_expression) %>% 
    unique() %>% 
    spread(gene_symbol,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_6_F <- df_test_6_F[,!is.na(colSums(df_test_6_F))]
  
  
  df_test_12_F <- df %>% 
    filter(complex_name == i,
           Age == 12,
           Sex == "F") %>% 
    select(mouse.id, gene_symbol,gene_expression) %>% 
    unique() %>% 
    spread(gene_symbol,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_12_F <- df_test_12_F[,!is.na(colSums(df_test_12_F))]
  
  
  df_test_18_F <- df %>% 
    filter(complex_name == i,
           Age == 18,
           Sex == "F") %>% 
    select(mouse.id, gene_symbol,gene_expression) %>% 
    unique() %>% 
    spread(gene_symbol,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_18_F <- df_test_18_F[,!is.na(colSums(df_test_18_F))]
  
  
  par(mfrow=c(2,3))
  
  COR_06_M <- cor(df_test_6_M, use = "pairwise.complete.obs")
  COR_12_M <- cor(df_test_12_M, use = "pairwise.complete.obs")
  COR_18_M <- cor(df_test_18_M, use = "pairwise.complete.obs")
  
  COR_06_F <- cor(df_test_6_F, use = "pairwise.complete.obs")
  COR_12_F <- cor(df_test_12_F, use = "pairwise.complete.obs")
  COR_18_F <- cor(df_test_18_F, use = "pairwise.complete.obs")
  
  CLUST <- hclust(dist(abs(COR_06_M)), method = 'ward.D2')
  ORDER <- CLUST$order
  
  order_rule[[i]] <- CLUST$labels[ORDER]
  
  corrplot(COR_06_M[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=1,tl.col = "black",col = col(200),
           title = paste0(i, " - Age 6 - Male"),mar=c(0,0,1.5,0))
  
  corrplot(COR_12_M[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=1,tl.col = "black",col = col(200),
           title = paste0(i, " - Age 12 - Male"),mar=c(0,0,1.5,0))
  
  corrplot(COR_18_M[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=1,tl.col = "black",col = col(200),
           title = paste0(i, " - Age 18 - Male"),mar=c(0,0,1.5,0))
  
  
  
  corrplot(COR_06_F[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=1,tl.col = "black",col = col(200),
           title = paste0(i, " - Age 6 - Female"),mar=c(0,0,1.5,0))
  
  corrplot(COR_12_F[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=1,tl.col = "black",col = col(200),
           title = paste0(i, " - Age 12 - Female"),mar=c(0,0,1.5,0))
  
  corrplot(COR_18_F[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=1,tl.col = "black",col = col(200),
           title = paste0(i, " - Age 18 - Female"),mar=c(0,0,1.5,0))
}

dev.off()


rm(df_test_6_M, df_test_12_M, df_test_18_M,
   df_test_6_F, df_test_12_F, df_test_18_F,
   COR_06_M, COR_12_M, COR_18_M,
   COR_06_F, COR_12_F, COR_18_F, CLUST, ORDER,i)

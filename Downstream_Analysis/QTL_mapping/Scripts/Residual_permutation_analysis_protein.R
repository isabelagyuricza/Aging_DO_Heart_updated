#################### QTL mapping analysis on the DO heart ######################

# Gathering permutations tests for proteins that I have run on the cluster 
# and gathering them.

################################################################################

# Author: Isabela Gerdes Gyuricza - Churchill lab 
# Date: 08-19-2021

################################################################################
############ clear workspace

rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# loading library

library(tidyverse)
library(stringr)

# loading data

load("QTLviewer_dataset/JAC_DO_heart_v9.gz")

rm(dataset.mrna, genoprobs, K, map, ensembl.version)

# Reading the results scped from Sumner:

result.list <- list.files(
  "Downstream_Analysis/QTL_mapping/Results/_rslurm_permutation_apply_protein", 
  pattern = "results_")

genes_list <- list()
for (i in result.list){
  genes_list[[i]] <- readRDS(
    paste0(
      "Downstream_Analysis/QTL_mapping/Results/_rslurm_permutation_apply_protein/",
      i))
}; rm(i)

# Checking the number of genes, it should be 2731 in total.

check_genes <- c()
for (i in 1:length(genes_list)){
  check_genes[i] <- length(genes_list[[i]])
}; rm(i)

sum(check_genes) #[1] 2731 Perfect, all proteins are here.

# Setting up the list:

lods_list <- unlist(genes_list,recursive = FALSE)

names(lods_list) <- str_split(names(lods_list),"RDS.",simplify = TRUE)[,2]

# Taking the 95% quantile for each gene.

for (i in 1:length(lods_list)){
  if (class(lods_list[[i]]) != "numeric") {print(i)}
}
rm(i, check_genes)

#Nothing, no problem.

# Using one lod threshold for all the genes based on a 95% threhsold. 

quantile(unlist(lods_list), c(0.37, 0.8, 0.9, 0.95, 0.99))
#    37%      80%      90%      95%      99% 
# 5.823390 6.846700 7.316771 7.752120 8.697736 

# Setting up lod peaks table with more information 

pqtls <- dataset.protein$lod.peaks$age_int %>%
  left_join(markers, by = 'marker.id') %>%
  rename(qtl.chr = chr,
         qtl.pos = pos) %>%
  select(protein.id, marker.id, lod, qtl.chr, qtl.pos) %>% 
  left_join(dataset.protein$annot.protein %>% 
              dplyr::select(protein.id, gene.id, chr, start, end, symbol), 
            by='protein.id') %>% 
  rename(gene.start = start,
         gene.end = end,
         gene.chr = chr,
         gene.symbol = symbol) %>% 
  mutate(cis = (gene.chr == qtl.chr) & abs(qtl.pos-gene.start) <= 4) %>%  # Considering 
# a 4Mb window for cis QTL.
  select(protein.id, gene.id, everything())

pqtls %>% 
  filter(lod > 7.75) %>% 
  dim()
# 603 significants pQTLs

pqtls %>% 
  filter(lod > 7.75, 
         cis) %>% 
  dim()

# None of these are local 

sig_QTLs <- pqtls %>% 
  filter(lod > 7.75)

saveRDS(sig_QTLs, 
        file = "Downstream_Analysis/QTL_mapping/Results/significant_pQTLs.RDS")
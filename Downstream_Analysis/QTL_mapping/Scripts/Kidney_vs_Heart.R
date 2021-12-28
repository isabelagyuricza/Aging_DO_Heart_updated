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

load("Downstream_analysis/Data/JAC_DO_kidney_v5_03_21_2020.RData")

# Cleaning up what I don't need

rm(genoprobs, K, map, ensembl.version,dataset.protein.zam)

lod.peaks_rna_kidney <- dataset.mrna$lod.peaks$Age_Int %>%
  left_join(markers, by = 'marker.id') %>%  # Get marker ids
  rename(qtl.chr   = chr,
                qtl.pos   = pos) %>%
  select(gene.id, marker.id, lod, qtl.chr, qtl.pos) %>% 
  left_join(dataset.mrna$annot.mrna %>% 
              dplyr::select(gene.id,chr,
                            start,end,symbol), 
            by='gene.id') %>% 
  rename(gene.start=start,
                gene.end=end,
                gene.chr=chr,
                gene.symbol=symbol) %>% 
  mutate(cis = (gene.chr == qtl.chr) & abs(qtl.pos-gene.start) <= 4)

lod.peaks_protein_kidney <- dataset.protein.zaz$lod.peaks$Age_Int %>%
  left_join(markers, by = 'marker.id') %>%  # Get marker ids
  rename(qtl.chr   = chr,
         qtl.pos   = pos) %>%
  select(protein.id, marker.id, lod, qtl.chr, qtl.pos) %>% 
  left_join(dataset.protein.zaz$annot.protein %>% 
              dplyr::select(protein.id,chr,
                            start,end,symbol), 
            by='protein.id') %>% 
  rename(gene.start=start,
         gene.end=end,
         gene.chr=chr,
         gene.symbol=symbol) %>% 
  mutate(cis = (gene.chr == qtl.chr) & abs(qtl.pos-gene.start) <= 4)

rm(dataset.mrna, dataset.protein.zaz, markers)


load("QTLviewer_dataset/JAC_DO_heart_v9.gz")

rm(genoprobs, K, map, ensembl.version)

lod.peaks_rna_heart <- dataset.mrna$lod.peaks$age_int %>%
  left_join(markers, by = 'marker.id') %>%  # Get marker ids
  rename(qtl.chr   = chr,
         qtl.pos   = pos) %>%
  select(gene.id, marker.id, lod, qtl.chr, qtl.pos) %>% 
  left_join(dataset.mrna$annot.mrna %>% 
              dplyr::select(gene.id,chr,
                            start,end,symbol), 
            by='gene.id') %>% 
  rename(gene.start=start,
         gene.end=end,
         gene.chr=chr,
         gene.symbol=symbol) %>% 
  mutate(cis = (gene.chr == qtl.chr) & abs(qtl.pos-gene.start) <= 4)

lod.peaks_protein_heart <- dataset.protein$lod.peaks$age_int %>%
  left_join(markers, by = 'marker.id') %>%  # Get marker ids
  rename(qtl.chr   = chr,
         qtl.pos   = pos) %>%
  select(protein.id, marker.id, lod, qtl.chr, qtl.pos) %>% 
  left_join(dataset.protein$annot.protein %>% 
              dplyr::select(protein.id,chr,
                            start,end,symbol), 
            by='protein.id') %>% 
  rename(gene.start=start,
         gene.end=end,
         gene.chr=chr,
         gene.symbol=symbol) %>% 
  mutate(cis = (gene.chr == qtl.chr) & abs(qtl.pos-gene.start) <= 4)

rm(dataset.mrna, dataset.protein, markers)

lod.peaks_protein_kidney %>% 
  filter(cis) %>% 
  inner_join(lod.peaks_protein_heart %>% 
               filter(cis), by = "protein.id")

lod.peaks_rna_kidney %>% 
  filter(cis) %>% 
  inner_join(lod.peaks_rna_heart %>% 
               filter(cis), by = "gene.id")

protein_distal_joint <- lod.peaks_protein_kidney %>% 
  filter(!cis) %>% 
  select(protein.id, qtl.chr, qtl.pos, gene.symbol, lod) %>% 
  rename(qtl.pos_kidney = qtl.pos,
         gene_symbol_kidney = gene.symbol,
         lod_kidney = lod) %>% 
  inner_join(
    lod.peaks_protein_heart %>% 
      filter(!cis) %>% 
      select(protein.id, qtl.chr, qtl.pos, gene.symbol,lod) %>% 
      rename(qtl.pos_heart = qtl.pos,
             gene_symbol_heart = gene.symbol,
             lod_heart = lod), by = c("protein.id", "qtl.chr")) %>% 
  filter(abs(qtl.pos_kidney - qtl.pos_heart) < 4)

#There are 41 common out of 10533 distal pQTLs in the kidney and 8007 in the heart

rna_distal_joint <- lod.peaks_rna_kidney %>% 
  filter(!cis) %>% 
  select(gene.id, qtl.chr, qtl.pos, gene.symbol,lod) %>% 
  rename(qtl.pos_kidney = qtl.pos,
         gene_symbol_kidney = gene.symbol,
         lod_kidney = lod) %>% 
  inner_join(
    lod.peaks_rna_heart %>% 
      filter(!cis) %>% 
      select(gene.id, qtl.chr, qtl.pos, gene.symbol,lod) %>% 
      rename(qtl.pos_heart = qtl.pos,
             gene_symbol_heart = gene.symbol,
             lod_heart = lod), by = c("gene.id", "qtl.chr")) %>% 
  filter(abs(qtl.pos_kidney - qtl.pos_heart) < 4)

# There are 189 common out of 27386 for kidney and 18668 for heart


###################### QTL analysis on the DO heart ############################

# Creating supplemental table

################################################################################

# Author: Isabela Gerdes Gyuricza - Churchill lab 
# Date: 11-15-2021

################################################################################
############ clear workspace

rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

################################################################################
############ libraries

library(tidyverse)

################################################################################
############ load data

load("QTLviewer_dataset/JAC_DO_heart_v9.gz")

chr3 <- readRDS("Downstream_Analysis/QTL_mapping/Results/proteins_hotspot_chr3.RDS")

chr12 <- readRDS("Downstream_Analysis/QTL_mapping/Results/proteins_hotspot_chr12.RDS")


eQTLs <- dataset.mrna$lod.peaks$age_int %>% 
         left_join(markers, by = 'marker.id') %>%
         rename(Gene_ID = gene.id, 
                Marker_ID = marker.id,
                Lod_score = lod,
                QTL_chr   = chr,
                QTL_pos   = pos) %>%
         select(Gene_ID, Marker_ID, Lod_score, QTL_chr, QTL_pos) %>% 
         left_join(dataset.mrna$annot.mrna %>% 
                     rename(Gene_ID = gene.id) %>% 
                     select(Gene_ID,chr,
                            start,end,symbol), 
            by='Gene_ID') %>% 
        rename(Gene_start=start,
               Gene_end=end,
               Gene_chr=chr,
               Gene_symbol=symbol) %>% 
        mutate(Cis_QTL = (Gene_chr == QTL_chr) & abs(QTL_pos-Gene_start) <= 4) %>% 
        mutate(Significant = ifelse(Lod_score > 7.75, TRUE, FALSE),
               Gene_symbol = paste0("\"", Gene_symbol, "\"")) %>% 
  as.data.frame()

write.csv(eQTLs, "Downstream_Analysis/QTL_mapping/Results/eQTLs_table.csv",
          row.names = FALSE)

hotspot <- rbind(
  chr3 %>% 
    select(protein.id, qtl.chr, qtl.pos, 
           lod, gene.id, symbol, chr,
           start, end, mean.corr) %>% 
    mutate(Hotspot_chr = 3,
           qtl.chr = as.character(qtl.chr)) %>% 
    rename(Protein_ID = protein.id, 
           QTL_chr = qtl.chr,
           QTL_pos = qtl.pos, 
           Lod_score = lod,
           Gene_ID = gene.id,
           Gene_symbol = symbol, 
           Gene_chr = chr,
           Gene_start = start,
           Gene_end = end,
           Mean_correlation = mean.corr),
  chr12 %>% 
    select(protein.id, qtl.chr, qtl.pos, 
           lod, gene.id, symbol, chr,
           start, end, mean.corr) %>% 
    mutate(Hotspot_chr = 12,
           qtl.chr = as.character(qtl.chr)) %>% 
    rename(Protein_ID = protein.id, 
           QTL_chr = qtl.chr,
           QTL_pos = qtl.pos, 
           Lod_score = lod,
           Gene_ID = gene.id,
           Gene_symbol = symbol, 
           Gene_chr = chr,
           Gene_start = start,
           Gene_end = end,
           Mean_correlation = mean.corr)
)

Hotspot_fun <- function(Hotspot_chr){

  if(
    sum(Hotspot_chr == c(3,12)) > 1 |
     sum(Hotspot_chr == c(12,3)) > 1) {
    
    r <- TRUE
    
    } else {
      
      r <- FALSE
    }
  
  return(r)
  
  }    

hotspot <- hotspot %>%
  group_by(Protein_ID) %>% 
  mutate(Both_hotspots = Hotspot_fun(Hotspot_chr)) %>% 
  ungroup()

pQTLs <- dataset.protein$lod.peaks$age_int %>% 
  left_join(markers, by = 'marker.id') %>%
  rename(Protein_ID = protein.id, 
         Marker_ID = marker.id,
         Lod_score = lod,
         QTL_chr   = chr,
         QTL_pos   = pos) %>%
  select(Protein_ID, Marker_ID, Lod_score, QTL_chr, QTL_pos) %>% 
  left_join(dataset.protein$annot.protein %>% 
              rename(Protein_ID = protein.id,
                     Gene_ID = gene.id) %>% 
              select(Protein_ID, Gene_ID,chr,
                     start,end,symbol), 
            by='Protein_ID') %>% 
  rename(Gene_start=start,
         Gene_end=end,
         Gene_chr=chr,
         Gene_symbol=symbol) %>% 
  mutate(Cis_QTL = (Gene_chr == QTL_chr) & abs(QTL_pos-Gene_start) <= 4) %>% 
  mutate(Significant = ifelse(Lod_score > 7.75, TRUE, FALSE)) %>% 
  select(Protein_ID, Gene_ID, Marker_ID, Lod_score, QTL_chr, QTL_pos, 
         Gene_chr, Gene_start, Gene_end, Gene_symbol, Cis_QTL, Significant) %>% 
  left_join(hotspot) %>% 
  mutate(Gene_symbol = paste0("\"", Gene_symbol, "\"")) %>% 
  as.data.frame()

write.csv(pQTLs, "Downstream_Analysis/QTL_mapping/Results/pQTLs_table.csv",
          row.names = FALSE)

annot <- readRDS("Downstream_Analysis/Data/all_gene_sets.RDS")

tmp <- annot %>% 
  filter(gs_subcat == "GO:BP") %>% 
  mutate(gs_name = factor(gs_name)) %>%
  dplyr::select(gs_name, ensembl_gene)

GO_Pathways <- split(tmp$ensembl_gene, tmp$gs_name)

GO_Pathways <- plyr::ldply(GO_Pathways, rbind) %>% 
  rename(pathway = .id) %>%
  gather("test","Gene_ID", -pathway) %>%
  select(-test) %>%
  filter(!is.na(Gene_ID))

GO_Pathways <- GO_Pathways[!duplicated(GO_Pathways$Gene_ID),]


hotspot <- hotspot %>% 
  left_join(GO_Pathways)

common <- intersect(chr3$gene.id, chr12$gene.id)
  
hotspot %>% filter(Gene_ID %in% common,
                   Hotspot_chr == 3 |
                     Hotspot_chr == 12) %>% View()

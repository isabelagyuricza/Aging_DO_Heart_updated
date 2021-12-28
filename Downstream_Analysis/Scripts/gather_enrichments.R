############## Enrichment analysis comparison on the DO heart ##################

# Comparing the enrichment analyis results between protein and transcripts data

################################################################################

# Author: Isabela Gerdes Gyuricza - Churchill lab 
# Date: 09-08-2021

################################################################################
############ clear workspace

rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

################################################################################
############ loading libraries

library(tidyverse)

# Loading enrichment results

transcript <- read.csv("Downstream_Analysis/Results/enrichment_age_transcripts.csv") %>% 
  mutate(Gene_IDs = as.list(Gene_IDs))

protein <- read.csv("Downstream_Analysis/Results/enrichment_age_proteins.csv") %>% 
  mutate(Gene_IDs = as.list(Gene_IDs))

intersect(transcript$Pathway_name, protein$Pathway_name)
#"GOMF_UNFOLDED_PROTEIN_BINDING"            

# There is only one common category between transcript and protein. 

# Let's try combining the pathways based on gene symbols.

all_gene_sets <- readRDS(file="Downstream_Analysis/Data/all_gene_sets.RDS")

# For proteins

DEA_proteins <- read_csv("Downstream_Analysis/Results/DEA_age_proteins.csv")

enrichment <- protein$Gene_IDs

names(enrichment) <-  protein$Pathway_name

# Transforming these into character vector

for (i in names(enrichment)) {
  
  tmp <- gsub("[[:punct:]]"," ",enrichment[[i]])
  
  tmp <- str_split(tmp,"   ", simplify = TRUE) %>% as.vector()
  
  tmp <- str_split(tmp,"  ", simplify = TRUE) %>% as.vector()
  
  tmp <- tmp[grepl("ENS",tmp)]
  
  enrichment[[i]] <- tmp
}

rm(i, tmp)

temp <- plyr::ldply(enrichment, rbind) %>%
  rename(Pathway_name = .id) %>%
  gather("test","ensembl_gene", -Pathway_name) %>%
  dplyr::select(-test) %>%
  filter(!is.na(ensembl_gene)) %>%
  left_join(all_gene_sets %>% select(ensembl_gene, gene_symbol)) %>%
  unique()

enrich_protein <- temp %>% 
  left_join(protein %>% 
              select(Pathway_name, adjusted_p)) %>% 
  arrange(adjusted_p) %>% 
  rename(adjusted_p_pathway = adjusted_p,
         Gene_ID = ensembl_gene, 
         Gene_symbol = gene_symbol) %>% 
  left_join(DEA_proteins) %>% 
  mutate(STD_Effect = Effect/Effect_SE) %>% 
  select(Pathway_name, Protein_ID, Gene_ID, Gene_symbol, STD_Effect, 
         adjusted_p, adjusted_p_pathway)

enrich_protein <- enrich_protein[!duplicated(enrich_protein$Gene_symbol),]

rm(enrichment)

# For transcripts

DEA_transcript <- read_csv("Downstream_Analysis/Results/DEA_age_transcripts.csv")

enrichment <- transcript$Gene_IDs

names(enrichment) <-  transcript$Pathway_name

# Transforming these into character vector

for (i in names(enrichment)) {
  
  tmp <- gsub("[[:punct:]]"," ",enrichment[[i]])
  
  tmp <- str_split(tmp,"   ", simplify = TRUE) %>% as.vector()
  
  tmp <- str_split(tmp,"  ", simplify = TRUE) %>% as.vector()
  
  tmp <- tmp[grepl("ENS",tmp)]
  
  enrichment[[i]] <- tmp
}

rm(i, tmp)

temp <- plyr::ldply(enrichment, rbind) %>%
  rename(Pathway_name = .id) %>%
  gather("test","ensembl_gene", -Pathway_name) %>%
  dplyr::select(-test) %>%
  filter(!is.na(ensembl_gene)) %>%
  left_join(all_gene_sets %>% select(ensembl_gene, gene_symbol)) %>%
  unique()

enrich_transcript <- temp %>% 
  left_join(transcript %>% 
              select(Pathway_name, adjusted_p)) %>% 
  arrange(adjusted_p) %>% 
  rename(adjusted_p_pathway = adjusted_p,
         Gene_ID = ensembl_gene, 
         Gene_symbol = gene_symbol) %>% 
  left_join(DEA_transcript) %>% 
  mutate(STD_Effect = Effect/Effect_SE) %>% 
  select(Pathway_name, Gene_ID, Gene_symbol, STD_Effect, 
         adjusted_p, adjusted_p_pathway)

enrich_transcript <- enrich_transcript[!duplicated(enrich_transcript$Gene_symbol),]

rm(enrichment)

# Joining the datasets

enrich <- enrich_protein %>% 
  rename(Pathway_name_protein = Pathway_name,
         STD_Effect_protein = STD_Effect,
         adjusted_p_protein = adjusted_p,
         adjusted_p_pathway_protein = adjusted_p_pathway) %>% 
  inner_join(enrich_transcript %>% 
               rename(Pathway_name_transcript = Pathway_name,
                      STD_Effect_transcript = STD_Effect,
                      adjusted_p_transcript = adjusted_p,
                      adjusted_p_pathway_transcript = adjusted_p_pathway)) %>% 
  select(Pathway_name_protein, Pathway_name_transcript, Protein_ID, Gene_ID, 
         Gene_symbol, STD_Effect_protein, STD_Effect_transcript, adjusted_p_protein,
         adjusted_p_transcript, adjusted_p_pathway_protein, adjusted_p_pathway_transcript)

# Saving for supplemental table 
# 
# write.csv(enrich, file = "Downstream_Analysis/Results/gathered_enrichment.csv",
#           row.names = FALSE)

# enrich %>% 
#   ggplot(aes(x = STD_Effect_transcript, y = STD_Effect_protein)) +
#   geom_point(aes(color = Pathway_name_transcript), size = 5, alpha = 0.6) +
#   geom_text(aes(label = Gene_symbol)) +
#   #scale_color_manual(values = c('#e41a1c','#377eb8','#4daf4a','#984ea3')) +
#   geom_vline(xintercept = 0, col="black", linetype = "dashed", size = 1) +
#   geom_hline(yintercept = 0, col="black", linetype = "dashed", size = 1) +
#   labs(x = "STD transcript age effect",
#        y = "STD protein age effect") +
#   #scale_x_continuous(limits = c(-6,6)) +
#   #scale_y_continuous(limits = c(-11,11)) +
#   theme_bw()


# There are too many categories to plot. Choosing only the ones for genes
# on the edge. 


pdf("Downstream_Analysis/Results/common_enrichment.pdf", 
    width = 22, height = 7)

path <- enrich %>% 
  filter(Gene_symbol %in% c("Rab21","Atp2a2","Lman1",
                            "Cd63","Ctsz","S100a9","Psmb9",
                            "Hspd1","Mical1","Cct8")) %>% 
  select(Pathway_name_transcript) %>% 
  unlist()

plot1 <- enrich %>% 
  mutate(Pathway_name_transcript = ifelse(Pathway_name_transcript %in% path,
                                          Pathway_name_transcript,
                                          "Other category")) %>% 
  ggplot(aes(x = STD_Effect_transcript, y = STD_Effect_protein)) +
  geom_point(aes(color = Pathway_name_transcript), size = 5, alpha = 0.6) +
  geom_text(aes(label = ifelse(STD_Effect_transcript > 1.5 &
                                 STD_Effect_transcript < 4.5 &
                                 STD_Effect_protein > 0 &
                                 STD_Effect_protein < 8.5 |
                                 STD_Effect_transcript < -1.5 &
                                 STD_Effect_transcript > -3.5 &
                                 STD_Effect_protein > 3 &
                                 STD_Effect_protein < 8, "",
                               Gene_symbol)),
                size = 6) +
  scale_color_manual(values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c',
                                '#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
                                '#cab2d6','#6a3d9a','grey')) +
  geom_vline(xintercept = 0, col="black", linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, col="black", linetype = "dashed", size = 1) +
  labs(x = "Standardized age effect - Transcript",
       y = "Standardized age effect - Protein") +
  scale_x_continuous(limits = c(-11,11)) +
  scale_y_continuous(limits = c(-11,11)) +
  theme_bw() +
  theme(legend.text = element_text(size = 2))


path <- enrich %>% 
  filter(Gene_symbol %in% c("Rab21","Atp2a2","Lman1",
                            "Cd63","Ctsz","S100a9","Psmb9",
                            "Hspd1","Mical1","Cct8")) %>% 
  select(Pathway_name_protein) %>% 
  unlist()


plot2 <- enrich %>% 
  mutate(Pathway_name_protein = ifelse(Pathway_name_protein %in% path,
                                          Pathway_name_protein,
                                          "Other category")) %>% 
  ggplot(aes(x = STD_Effect_transcript, y = STD_Effect_protein)) +
  geom_point(aes(color = Pathway_name_protein), size = 5, alpha = 0.6) +
  geom_text(aes(label = ifelse(STD_Effect_transcript > 1.5 &
                                 STD_Effect_transcript < 4.5 &
                                 STD_Effect_protein > 0 &
                                 STD_Effect_protein < 8.5 |
                                 STD_Effect_transcript < -1.5 &
                                 STD_Effect_transcript > -3.5 &
                                 STD_Effect_protein > 3 &
                                 STD_Effect_protein < 8, "",
                               Gene_symbol)),
                size = 6) +
  scale_color_manual(values = c('#e41a1c','#377eb8','#4daf4a','#984ea3',
                                '#ff7f00','#e7298a','grey')) +
  geom_vline(xintercept = 0, col="black", linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, col="black", linetype = "dashed", size = 1) +
  labs(x = "Standardized age effect - Transcript",
       y = "Standardized age effect - Protein") +
  scale_x_continuous(limits = c(-11,11)) +
  scale_y_continuous(limits = c(-11,11)) +
  theme_bw() +
  theme(legend.text = element_text(size = 2))

cowplot::plot_grid(plot1,plot2,ncol = 2)

dev.off()
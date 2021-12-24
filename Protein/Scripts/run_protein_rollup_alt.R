###################### Processing peptide data #################################

# Processing peptide raw data by removing polymorphic peptides adapted from 
# Greg Keele (Churchill lab). Also removing the two problematic mice 
# (DO.1077 and DO.1157)

################################################################################

# Author: Isabela Gerdes Gyuricza - Churchill lab 
# Date: 06_30_2021


################################################################################
################################################################################
######
######                Title: run_protein_rollup.R
######                Description: Rolls up protein abundances from summarized
######                peptide quantifications to be used for further analysis
######                
######                Manuscript: Keele & Zhang et al. 
######
######                Author: Greg Keele
######
################################################################################
################################################################################


################################################################################
############ clear workspace

rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

################################################################################
############ loading libraries

library(ensimplR)
library(tidyverse)

options(stringsAsFactors = FALSE)

## Set working directory
setwd("~/Box/JAC_Heart_Data/Heart_Data_June2021")

## Support functions
source("Protein/Scripts/protein_rollup.R")
library(tidyverse)

## Polymorphic peptides
polymorph_dat <- read.csv("Protein/Data/peptides_w_polymorphisms.csv", 
                          header = TRUE) %>%
  mutate(peptide.id = gsub(x = peptide.id, pattern = "^[-|A-Z]+\\.", 
                           replacement = "", perl = TRUE),
         ## Remove second period and post-string
         peptide.id = gsub(x = peptide.id, pattern = "\\.[A-Z|-]+$", 
                           replacement = "", perl = TRUE),
         ## Remove * representing differently charged methyl groups
         peptide.id = gsub(x = peptide.id, pattern = "*", 
                           replacement = "", fixed = TRUE),
         ## Remove -, always at beginning or end of AA sequence
         peptide.id = gsub(x = peptide.id, pattern = "-", 
                           replacement = "", fixed = TRUE))
polymorph_dat$balance <- rowSums(polymorph_dat[,LETTERS[1:8]])

## Grab variant-to-expectation correlations from CC multi-tissue data
cc_peptide_cor <- readRDS("Protein/Data/cc_multitissue_cor_dat.RDS")

## Filtering out peptides based on the data matching our expectations (cor > 0.7)
polymorph_remove <- cc_peptide_cor %>%
  filter(!is.na(heart_cor) & heart_cor > 0.7) %>%
  pull(peptide.id)

################################################################################ 
## 
##                             DO roll-up
##
################################################################################
## DO annotations
do_annot <- read.csv("Protein/Data/SampleAnnotations.csv", header = TRUE) %>%
  dplyr::rename(mouse.id = Mouse.ID,
                tag = Tag) %>%
  mutate(Batch = gsub(x = Batch, pattern = "Set", replacement = "TMT")) %>%
  mutate(sample.id = paste(Batch, tag, sep = "~")) %>%
  dplyr::select(-c(Sample.Number, NYGC.ID))

## Raw peptide data in CC
raw_do_peptide_dat <- read_tsv("Protein/Data/heart_samples_peptides.tsv")

## Process DO peptides
do_peptide_dat <- raw_do_peptide_dat %>%
  dplyr::select(ProteinId, PeptideSequence, 
                starts_with("rq"), Class, -Description) %>%
  # Filter out reverse proteins
  filter(!grepl(x = ProteinId, pattern = "##", fixed = TRUE),
         # Filter out contaminants
         !grepl(x = ProteinId, pattern = "contaminant", fixed = TRUE),
         # Filter to only things in Ensembl
         grepl(x = ProteinId, pattern = "ENSMUSP", fixed = TRUE)) %>%
  # Convert TMT-channels from columns to rows
  gather(starts_with("rq"), key = "sample_id", value = "Intensity") %>%
  # Process peptide sequence
  ## Remove first period and pre-string
  mutate(PeptideSequence = gsub(x = PeptideSequence, pattern = "^[-|A-Z]+\\.", 
                                replacement = "", perl = TRUE),
         ## Remove second period and post-string
         PeptideSequence = gsub(x = PeptideSequence, pattern = "\\.[A-Z|-]+$", 
                                replacement = "", perl = TRUE),
         ## Remove * representing differently charged methyl groups
         PeptideSequence = gsub(x = PeptideSequence, pattern = "*", 
                                replacement = "", fixed = TRUE),
         ## Remove -, always at beginning or end of AA sequence
         PeptideSequence = gsub(x = PeptideSequence, pattern = "-", 
                                replacement = "", fixed = TRUE),
         ## Remove "." and number from ProteinId
         protein.id = gsub(ProteinId, pattern = "\\.[0-9]*$", 
                           replacement = "", perl = TRUE),
         ## Pulling TMT tag
         tag = gsub(x = sample_id, pattern = "rq_", 
                    replacement = "", perl = TRUE),
         tag = gsub(x = tag, pattern = "_sn", replacement = "", fixed = TRUE),
         Class = gsub(x = Class, pattern = "HS",
                      replacement = "TMT"),
         sample.id = paste0(Class, "~", tag)) %>%
  ## Rename PeptideSequence as peptide.id
  dplyr::rename(peptide.id = PeptideSequence,
                Batch = Class) %>%
  # Sum up multiple observations of a peptide within an experiment (per Tian)
  group_by(protein.id, peptide.id, sample.id, Batch, tag) %>%
  summarize(Intensity = sum(Intensity, na.rm = TRUE)) %>%
  ungroup

## Merging in polymorphic data to the CC
do_peptide_dat <- do_peptide_dat %>%
  left_join(polymorph_dat %>%
              dplyr::select(peptide.id, balance) %>%
              mutate(polymorphic = TRUE)) %>%
  mutate(balance = ifelse(is.na(balance), 8, balance),
         polymorphic = ifelse(is.na(polymorphic), FALSE, polymorphic))

# ## Number of polymorphic peptides originally being filtered
do_peptide_dat %>%
  dplyr::select(peptide.id, polymorphic) %>%
  distinct %>%
  filter(polymorphic) %>%
  nrow
# 7299
# 
# ## Number of polymorphic peptides that are filtered based on cor > 0.7
do_peptide_dat %>%
  dplyr::select(peptide.id) %>%
  distinct %>%
  inner_join(data.frame(peptide.id = polymorph_remove)) %>%
  nrow
# 936

## Filtering all polymorphic peptides
# do_protein_dat_old <- protein_rollup(long_peptide_dat = do_peptide_dat, 
#                                  filter_out_peptides = polymorph_dat$peptide.id, 
#                                  bridge_sample = NULL, 
#                                  zero_as_na = FALSE, 
#                                  na_as_zero = FALSE) %>%
#   left_join(do_annot %>%
#               dplyr::select(sample.id, mouse.id, Sex, Age, Generation)) %>% 
#   # Removing the two mice that were targetted in the genoprobs QC (DO.1077 and DO.1157)
#   dplyr::filter(!mouse.id %in% c("DO-1077","DO-1157")) 

## Filtering only polymorphic peptides that matched expectation in CC heart
do_protein_dat_new <- protein_rollup(long_peptide_dat = do_peptide_dat, 
                                     filter_out_peptides = polymorph_remove, 
                                     bridge_sample = NULL, 
                                     zero_as_na = FALSE, 
                                     na_as_zero = FALSE) %>%
  left_join(do_annot %>%
              dplyr::select(sample.id, mouse.id, Sex, Age, Generation)) %>% 
  # Removing the two mice that were targetted in the genoprobs QC (DO.1077 and DO.1157)
  dplyr::filter(!mouse.id %in% c("DO-1077","DO-1157")) 

# 
# do_protein_dat_old %>% 
#   group_by(protein.id) %>%
#   dplyr::summarize(na_perc = mean(is.na(Intensity))) %>%
#   pull(na_perc) %>%
#   hist
# do_protein_dat_old %>% 
#   group_by(protein.id) %>%
#   dplyr::summarize(zero_perc = sum(Intensity == 0, 
#                                    na.rm = T)/length(Intensity)) %>%
#   pull(zero_perc) %>%
#   hist

## Batch correction
source("Protein/Scripts/debatch_functions.R")

# debatched_do_protein_dat_old <- debatch_protein_lmer(long_dat = do_protein_dat_old %>%
#                                                    mutate(Age = factor(Age, levels = c(6, 12, 18))), 
#                                                  regress_out_sex = FALSE, 
#                                                  filter = 0.5, 
#                                                  filter_value = "NA") %>% 
#   gather(key = protein.id, value = Intensity, 
#          -c(mouse.id, Age, Sex, Generation, Batch, tag)) 

debatched_do_protein_dat_new <- debatch_protein_lmer(long_dat = do_protein_dat_new %>%
                                                       mutate(Age = factor(Age, levels = c(6, 12, 18))), 
                                                     regress_out_sex = FALSE, 
                                                     filter = 0.5, 
                                                     filter_value = "NA") %>% 
  gather(key = protein.id, value = Intensity, 
         -c(mouse.id, Age, Sex, Generation, Batch, tag)) 


# debatched_do_protein_dat_old %>% pull(protein.id) %>% unique() %>% length()
# # [1] 3960 proteins in total

debatched_do_protein_dat_new %>% pull(protein.id) %>% unique() %>% length()
# [1] 4223 proteins in total

## Correlate old and new data
# debatched_do_protein_mat_old <- debatched_do_protein_dat_old %>%
#   dplyr::select(mouse.id, protein.id, Intensity) %>%
#   pivot_wider(names_from = protein.id, values_from = Intensity) %>%
#   column_to_rownames("mouse.id") %>%
#   as.matrix
# debatched_do_protein_mat_new <- debatched_do_protein_dat_new %>%
#   dplyr::select(mouse.id, protein.id, Intensity) %>%
#   pivot_wider(names_from = protein.id, values_from = Intensity) %>%
#   column_to_rownames("mouse.id") %>%
#   as.matrix

# saveRDS(debatched_do_protein_mat_new, "~/Box/Heart_Data_June2021/Protein/Data/debatched_do_protein_mat_new.rds")
# 
# overlapping_proteins <- intersect(colnames(debatched_do_protein_mat_old), colnames(debatched_do_protein_mat_new))
# 
# overlapping_cor <- sapply(1:length(overlapping_proteins), function(i) cor(debatched_do_protein_mat_old[,overlapping_proteins[i]], 
#                                                                           debatched_do_protein_mat_new[,overlapping_proteins[i]]))
# hist(overlapping_cor)


# Adjusting mouse.id to match the rna and genoprobs

debatched_do_protein_dat <- debatched_do_protein_dat_new %>% 
  dplyr::mutate(mouse.id = gsub("-",".",mouse.id))

# Saving debatched protein data

write.csv(debatched_do_protein_dat, 
          file = "Protein/Results/debatch_nopoly_protein_heart.csv", 
          quote = FALSE, row.names = FALSE)

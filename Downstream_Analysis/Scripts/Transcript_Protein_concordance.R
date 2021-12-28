################ Transcript-protein correlation analysis #######################

################################################################################

# Author: Isabela Gerdes Gyuricza - Churchill lab 
# Date: 11-22-2021

################################################################################
############ clear workspace

rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

################################################################################
############ loading libraries

library(tidyverse)

# Load QTLviewer data 

load("QTLviewer_dataset/JAC_DO_heart_v9.gz")

rm(genoprobs, K, map, markers, ensembl.version)

result_age_transcript <- read_csv(
  "Downstream_Analysis/Results/DEA_age_transcripts.csv"
  ) %>% 
  mutate(STD_Effect_Transcript = Effect/Effect_SE) %>% 
  rename(adjusted_p_transcript = adjusted_p)

result_age_protein <- read_csv(
  "Downstream_Analysis/Results/DEA_age_proteins.csv"
  ) %>% 
  mutate(STD_Effect_Protein = Effect/Effect_SE) %>% 
  rename(adjusted_p_protein = adjusted_p)


result_age <- result_age_transcript %>% 
  select(Gene_ID, 
         Gene_symbol, 
         STD_Effect_Transcript, 
         adjusted_p_transcript) %>% 
  inner_join(result_age_protein %>% 
               select(Gene_ID, Protein_ID, Gene_symbol, 
                      STD_Effect_Protein, adjusted_p_protein)) %>% 
  mutate(group.age="Non significant") %>%
  mutate(group.age = ifelse(
    (
      (adjusted_p_transcript < 0.1) & (adjusted_p_protein < 0.1) &
         (STD_Effect_Transcript > 0) & (STD_Effect_Protein > 0)
      ), 
    "A", 
    group.age)) %>%
  mutate(group.age = ifelse(
    (
      (adjusted_p_transcript < 0.1) & (adjusted_p_protein < 0.1) &
         (STD_Effect_Transcript < 0) & (STD_Effect_Protein > 0)
      ), 
    "B", 
    group.age)) %>%
  mutate(group.age = ifelse(
    (
      (adjusted_p_transcript < 0.1) & (adjusted_p_protein < 0.1) &
         (STD_Effect_Transcript > 0) & (STD_Effect_Protein < 0)
      ), 
    "C",
    group.age)) %>%
  mutate(group.age = ifelse(
    (
      (adjusted_p_transcript < 0.1) & (adjusted_p_protein < 0.1) &
         (STD_Effect_Transcript < 0) & (STD_Effect_Protein < 0)
      ), 
    "D", 
    group.age)) %>%
  mutate(group.age = factor(
    group.age, levels=c("Non significant",
                        "A",
                        "B",
                        "C",
                        "D"))
    )

table(result_age$group.age)
# Non significant               A               B               C               D 
#       3524                  290             131              39              74 


plot1 <- result_age %>% 
ggplot(
  aes(
    x=STD_Effect_Transcript, y=STD_Effect_Protein, color = group.age)) +
  geom_point(alpha = 0.5, size = 3) +
  scale_color_manual(values = c('grey','#377eb8','#4daf4a','#984ea3','#e41a1c')) +
  scale_x_continuous(limits = c(-10,10), breaks = c(-10,-5,0,-5,10)) +
  scale_y_continuous(limits = c(-20,20)) +
  guides(color = "none") +
  theme_bw()

# Now for sex 

result_sex_transcript <- read_csv(
  "Downstream_Analysis/Results/DEA_sex_transcripts.csv"
) %>% 
  mutate(STD_Effect_Transcript = Effect/Effect_SE) %>% 
  rename(adjusted_p_transcript = adjusted_p)

result_sex_protein <- read_csv(
  "Downstream_Analysis/Results/DEA_sex_proteins.csv"
) %>% 
  mutate(STD_Effect_Protein = Effect/Effect_SE) %>% 
  rename(adjusted_p_protein = adjusted_p)

result_sex <- result_sex_transcript %>% 
  select(Gene_ID, 
         Gene_symbol, 
         STD_Effect_Transcript, 
         adjusted_p_transcript) %>% 
  inner_join(result_sex_protein %>% 
               select(Gene_ID, Protein_ID, Gene_symbol, 
                      STD_Effect_Protein, adjusted_p_protein)) %>% 
  left_join(dataset.mrna$annot.mrna %>% 
              select(gene.id, symbol, chr) %>% 
              rename(Gene_ID = gene.id, 
                     Gene_symbol = symbol)) %>% 
  mutate(group.sex="Non significant") %>%
  mutate(group.sex = ifelse(
    (
      (adjusted_p_transcript < 0.1) & (adjusted_p_protein < 0.1) &
        (STD_Effect_Transcript > 0) & (STD_Effect_Protein > 0)
    ), 
    "A", 
    group.sex)) %>%
  mutate(group.sex = ifelse(
    (
      (adjusted_p_transcript < 0.1) & (adjusted_p_protein < 0.1) &
        (STD_Effect_Transcript < 0) & (STD_Effect_Protein > 0)
    ), 
    "B", 
    group.sex)) %>%
  mutate(group.sex = ifelse(
    (
      (adjusted_p_transcript < 0.1) & (adjusted_p_protein < 0.1) &
        (STD_Effect_Transcript > 0) & (STD_Effect_Protein < 0)
    ), 
    "C",
    group.sex)) %>%
  mutate(group.sex = ifelse(
    (
      (adjusted_p_transcript < 0.1) & (adjusted_p_protein < 0.1) &
        (STD_Effect_Transcript < 0) & (STD_Effect_Protein < 0)
    ), 
    "D", 
    group.sex)) %>%
  mutate(group.sex = factor(
    group.sex, levels=c("Non significant",
                        "A",
                        "B",
                        "C",
                        "D"))
  )

table(result_sex$group.sex)
# Non significant               A               B               C               D 
#       3549                   171              16              48             274

plot2 <- result_sex %>% 
  filter(!chr %in% c("X","Y")) %>%  #Removing sexual genes for plotting.
  ggplot(
    aes(
      x=STD_Effect_Transcript, y=STD_Effect_Protein, color = group.sex)) +
  geom_point(alpha = 0.5, size = 3) +
  scale_color_manual(values = c('grey','#377eb8','#4daf4a','#984ea3','#e41a1c')) +
  scale_x_continuous(limits = c(-25,25)) +
  scale_y_continuous(limits = c(-20,20)) +
  guides(color = "none") +
  theme_bw()


pdf("Downstream_Analysis/Results/age_sex_effects.pdf", 
    width = 14, height = 6)

cowplot::plot_grid(plot1, plot2, ncol = 2)

dev.off()

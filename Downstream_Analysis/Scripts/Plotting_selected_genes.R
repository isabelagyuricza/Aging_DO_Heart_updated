############################ Plotting genes ####################################

## Isabela Gerdes Gyuricza - 11_23_2021

################################################################################
############ clear workspace

rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

################################################################################
############ libraries

library(tidyverse)

################################################################################   
# Load data ####

load("QTLviewer_dataset/JAC_DO_heart_v9.gz")

rm(genoprobs, K, map, markers, ensembl.version)

# Age effect 

transcript <- dataset.mrna$data$norm %>% 
  data.frame() %>% 
  rownames_to_column("mouse.id") %>% 
  gather("gene.id","Expression",-mouse.id) %>% 
  left_join(dataset.mrna$annot.mrna %>% 
              select(gene.id, symbol)
  ) %>% 
  filter(symbol %in% c("Ryr1","S100a9","Hspa1b")) %>% #Selecting a few genes with high age effects
  left_join(dataset.mrna$annot.samples) %>% 
  mutate(Sex = ifelse(Sex == "F","Female","Male")) %>% 
  mutate(symbol = factor(symbol, levels = c("Ryr1","S100a9","Hspa1b")))

plot1 <- transcript %>% 
  ggplot(aes(x = factor(Age, levels = c(6,12,18)), y = Expression)) +
  geom_boxplot(outlier.shape = NA, aes(color = Sex, fill = Sex), alpha = 0.5) +
  geom_point(alpha = 0.5, size = 2,
             aes(color = Sex), 
                 position = position_jitterdodge(jitter.width = 0.2)) +
  labs(x = "Age", y = "Normalized expression") +
  theme_bw() +
  facet_wrap(~ symbol , scales = "free",ncol = 3) +
  ggtitle("Transcripts - age effect")

# Sex effect

transcript <- dataset.mrna$data$norm %>% 
  data.frame() %>% 
  rownames_to_column("mouse.id") %>% 
  gather("gene.id","Expression",-mouse.id) %>% 
  left_join(dataset.mrna$annot.mrna %>% 
              select(gene.id, symbol)
  ) %>% 
  filter(symbol %in% c("C7","Pbdc1","Ace")) %>% #Selecting a few genes with high sex effects
  left_join(dataset.mrna$annot.samples) %>% 
  mutate(Sex = ifelse(Sex == "F","Female","Male")) %>% 
  mutate(symbol = factor(symbol, levels = c("C7","Pbdc1","Ace")))

plot2 <- transcript %>% 
  mutate(Age = factor(Age, levels = c(6, 12, 18))) %>% 
  ggplot(aes(x = Sex, y = Expression)) +
  geom_boxplot(outlier.shape = NA,
               aes(color = Age, 
                   fill = Age), 
              alpha = 0.5) +
  geom_point(alpha = 0.5, 
             size = 2,
             aes(color = Age), 
             position = position_jitterdodge(jitter.width = 0.2)) +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  scale_fill_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  labs(x = "Sex", y = "Normalized expression") +
  theme_bw() +
  facet_wrap(~ symbol , scales = "free",ncol = 3) +
  ggtitle("Transcripts - sex effect")

# Age by sex effect 

transcript <- dataset.mrna$data$norm %>% 
  data.frame() %>% 
  rownames_to_column("mouse.id") %>% 
  gather("gene.id","Expression",-mouse.id) %>% 
  left_join(dataset.mrna$annot.mrna %>% 
              select(gene.id, symbol)
  ) %>% 
  filter(symbol %in% c("Hspa1b","Smpx","Gm4841")) %>% #Selecting a few genes with high age by sex effects
  left_join(dataset.mrna$annot.samples) %>% 
  mutate(Sex = ifelse(Sex == "F","Female","Male")) %>% 
  mutate(symbol = factor(symbol, levels = c("Hspa1b","Smpx","Gm4841")))

plot3 <- transcript %>% 
  ggplot(aes(x = factor(Age, levels = c(6,12,18)), y = Expression)) +
  geom_boxplot(outlier.shape = NA, aes(color = Sex, fill = Sex), alpha = 0.5) +
  geom_point(alpha = 0.5, size = 2,
             aes(color = Sex), 
             position = position_jitterdodge(jitter.width = 0.2)) +
  labs(x = "Age", y = "Normalized expression") +
  theme_bw() +
  facet_wrap(~ symbol , scales = "free",ncol = 3) +
  ggtitle("Transcripts - age by sex interaction effect")


# Age effect 

protein <- dataset.protein$data$norm %>% 
  data.frame() %>% 
  rownames_to_column("mouse.id") %>% 
  gather("protein.id","Expression",-mouse.id) %>% 
  left_join(dataset.protein$annot.protein %>% 
              select(protein.id, symbol)
  ) %>% 
  filter(symbol %in% c("Sepsecs","Rtn1","Actn3")) %>% #Selecting a few proteins with high age effects
  left_join(dataset.protein$annot.samples) %>% 
  mutate(Sex = ifelse(Sex == "F","Female","Male")) %>% 
  mutate(symbol = factor(symbol, levels = c("Sepsecs","Rtn1","Actn3")))

plot4 <- protein %>% 
  ggplot(aes(x = factor(Age, levels = c(6,12,18)), y = Expression)) +
  geom_boxplot(outlier.shape = NA, aes(color = Sex, fill = Sex), alpha = 0.5) +
  geom_point(alpha = 0.5, size = 2,
             aes(color = Sex), 
             position = position_jitterdodge(jitter.width = 0.2)) +
  labs(x = "Age", y = "Normalized abundance") +
  theme_bw() +
  facet_wrap(~ symbol , scales = "free",ncol = 3) +
  ggtitle("Proteins - age effect")

# Sex effect 

protein <- dataset.protein$data$norm %>% 
  data.frame() %>% 
  rownames_to_column("mouse.id") %>% 
  gather("protein.id","Expression",-mouse.id) %>% 
  left_join(dataset.protein$annot.protein %>% 
              select(protein.id, symbol)
  ) %>% 
  filter(symbol %in% c("C8g","A1bg","Mug1")) %>% #Selecting a few proteins with high sex effects
  left_join(dataset.protein$annot.samples) %>% 
  mutate(Sex = ifelse(Sex == "F","Female","Male")) %>% 
  mutate(symbol = factor(symbol, levels = c("C8g","A1bg","Mug1")))

plot5 <- protein %>% 
  mutate(Age = factor(Age, levels = c(6, 12, 18))) %>% 
  ggplot(aes(x = Sex, y = Expression)) +
  geom_boxplot(outlier.shape = NA,
               aes(color = Age, 
                   fill = Age), 
               alpha = 0.5) +
  geom_point(alpha = 0.5, 
             size = 2,
             aes(color = Age), 
             position = position_jitterdodge(jitter.width = 0.2)) +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  scale_fill_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  labs(x = "Sex", y = "Normalized abundance") +
  theme_bw() +
  facet_wrap(~ symbol , scales = "free",ncol = 3) +
  ggtitle("Proteins - sex effect")

# Age by sex effect 

protein <- dataset.protein$data$norm %>% 
  data.frame() %>% 
  rownames_to_column("mouse.id") %>% 
  gather("protein.id","Expression",-mouse.id) %>% 
  left_join(dataset.protein$annot.protein %>% 
              select(protein.id, symbol)
  ) %>% 
  filter(symbol %in% c("Fth1","Hnrnpl","Actg2")) %>% #Selecting a few proteins with high age by sex effects
  left_join(dataset.protein$annot.samples) %>% 
  mutate(Sex = ifelse(Sex == "F","Female","Male")) %>% 
  mutate(symbol = factor(symbol, levels = c("Fth1","Hnrnpl","Actg2")))

plot6 <- protein %>% 
  ggplot(aes(x = factor(Age, levels = c(6,12,18)), y = Expression)) +
  geom_boxplot(outlier.shape = NA, aes(color = Sex, fill = Sex), alpha = 0.5) +
  geom_point(alpha = 0.5, size = 2,
             aes(color = Sex), 
             position = position_jitterdodge(jitter.width = 0.2)) +
  labs(x = "Age", y = "Normalized abundance") +
  theme_bw() +
  facet_wrap(~ symbol , scales = "free",ncol = 3) +
  ggtitle("Proteins - age by sex interaction effect")


pdf("Downstream_Analysis/Results/selected_expression.pdf", 
    width = 12, 
    height = 10)

cowplot::plot_grid(plot1, plot2, plot3, ncol = 1)
cowplot::plot_grid(plot4, plot5, plot6, ncol = 1)


dev.off()


############## Now, for transcripts and proteins for figure 4 ##################

transcript <- dataset.mrna$data$norm %>% 
  data.frame() %>% 
  rownames_to_column("mouse.id") %>% 
  gather("gene.id","Expression",-mouse.id) %>% 
  left_join(dataset.mrna$annot.mrna %>% 
              select(gene.id, symbol)
  ) %>% 
  filter(symbol %in% c("Rab8a","Cyb5r3","Hspd1","Akt2")) %>% 
  left_join(dataset.mrna$annot.samples) %>% 
  mutate(Sex = ifelse(Sex == "F","Female","Male"),
         symbol = factor(symbol, levels = c("Rab8a","Cyb5r3","Hspd1","Akt2")),
         Type = "Transcript") %>% 
  select(mouse.id, symbol, Sex, Age, Expression, Type)


protein <- dataset.protein$data$norm %>% 
  data.frame() %>% 
  rownames_to_column("mouse.id") %>% 
  gather("protein.id","Expression",-mouse.id) %>% 
  left_join(dataset.protein$annot.protein %>% 
              select(protein.id, symbol)
  ) %>% 
  filter(symbol %in% c("Rab8a","Cyb5r3","Hspd1","Akt2")) %>%
  left_join(dataset.protein$annot.samples) %>% 
  mutate(Sex = ifelse(Sex == "F","Female","Male"),
         symbol = factor(symbol, levels = c("Rab8a","Cyb5r3","Hspd1","Akt2")),
         Type = "Protein") %>% 
  select(mouse.id, symbol, Sex, Age, Expression, Type)

df <- rbind(
  transcript,
  protein
) %>% 
  mutate(Age = factor(Age, levels = c(6, 12, 18)),
         Type = factor(Type, levels = c("Transcript","Protein")))

# Akt2
plot1 <-  df %>% 
  filter(symbol == "Akt2") %>% 
  ggplot(aes(x = Age, y = Expression)) +
  geom_boxplot(outlier.shape = NA, aes(color = Age, fill = Age), alpha = 0.3) +
  geom_point(aes(color = Age), 
             position = position_jitter(height = 0, width = 0.1), 
             alpha = 0.5) +
  geom_smooth(aes(x = as.integer(Age), y = Expression), method="lm",
              color="black", linetype="dashed") +
  labs(x = "Age", y = "Normalized expression") +
  theme_bw() +
  facet_grid(Type ~ Sex, scales = "free") +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  scale_fill_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  ggtitle("Akt2")

subset <- df %>% 
  filter(symbol == "Akt2") %>% 
  spread("Type","Expression")

plot2 <- subset %>% 
  ggplot(aes(x = Transcript, y = Protein, color = Age)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  ggtitle("")

# Rab8a
 plot3 <-  df %>% 
  filter(symbol == "Rab8a") %>% 
  ggplot(aes(x = Age, y = Expression)) +
  geom_boxplot(outlier.shape = NA, aes(color = Age, fill = Age), alpha = 0.3) +
  geom_point(aes(color = Age), 
             position = position_jitter(height = 0, width = 0.1), 
             alpha = 0.5) +
  geom_smooth(aes(x = as.integer(Age), y = Expression), method="lm",
              color="black", linetype="dashed") +
  labs(x = "Age", y = "Normalized expression") +
  theme_bw() +
  facet_grid(Type ~ Sex, scales = "free") +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  scale_fill_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  ggtitle("Rab8a")

subset <- df %>% 
  filter(symbol == "Rab8a") %>% 
  spread("Type","Expression")
  
plot4 <- subset %>% 
  ggplot(aes(x = Transcript, y = Protein, color = Age)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  ggtitle("")

# Cyb5r3 
plot5 <-  df %>% 
  filter(symbol == "Cyb5r3") %>% 
  ggplot(aes(x = Age, y = Expression)) +
  geom_boxplot(outlier.shape = NA, aes(color = Age, fill = Age), alpha = 0.3) +
  geom_point(aes(color = Age), 
             position = position_jitter(height = 0, width = 0.1), 
             alpha = 0.5) +
  geom_smooth(aes(x = as.integer(Age), y = Expression), method="lm",
              color="black", linetype="dashed") +
  labs(x = "Age", y = "Normalized expression") +
  theme_bw() +
  facet_grid(Type ~ Sex, scales = "free") +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  scale_fill_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  ggtitle("Cyb5r3")

subset <- df %>% 
  filter(symbol == "Cyb5r3") %>% 
  spread("Type","Expression")

plot6 <- subset %>% 
  ggplot(aes(x = Transcript, y = Protein, color = Age)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  ggtitle("")


# Hspd1 
plot7 <-  df %>% 
  filter(symbol == "Hspd1") %>% 
  ggplot(aes(x = Age, y = Expression)) +
  geom_boxplot(outlier.shape = NA, aes(color = Age, fill = Age), alpha = 0.3) +
  geom_point(aes(color = Age), 
             position = position_jitter(height = 0, width = 0.1), 
             alpha = 0.5) +
  geom_smooth(aes(x = as.integer(Age), y = Expression), method="lm",
              color="black", linetype="dashed") +
  labs(x = "Age", y = "Normalized expression") +
  theme_bw() +
  facet_grid(Type ~ Sex, scales = "free") +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  scale_fill_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  ggtitle("Hspd1")

subset <- df %>% 
  filter(symbol == "Hspd1") %>% 
  spread("Type","Expression")

plot8 <- subset %>% 
  ggplot(aes(x = Transcript, y = Protein, color = Age)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  ggtitle("")

pdf("Downstream_Analysis/Results/Multi_genes_expression.pdf", width = 9.5,
    height = 13)

cowplot::plot_grid(plot1, plot2, 
                   plot3, plot4,
                   plot5, plot6,
                   plot7, plot8,
                   ncol = 2)
dev.off()

################################################################################
# Supplemental figure for protein complex analysis

transcript <- dataset.mrna$data$norm %>% 
  data.frame() %>% 
  rownames_to_column("mouse.id") %>% 
  gather("gene.id","Expression",-mouse.id) %>% 
  left_join(dataset.mrna$annot.mrna %>% 
              select(gene.id, symbol)
  ) %>% 
  filter(symbol %in% c("Psmb8","Psmb9","Psmb3","Psmd7")) %>% 
  left_join(dataset.mrna$annot.samples) %>% 
  mutate(Sex = ifelse(Sex == "F","Female","Male"),
         symbol = factor(symbol, levels = c("Psmb8","Psmb9","Psmb3","Psmd7")),
         Type = "Transcript") %>% 
  select(mouse.id, symbol, Sex, Age, Expression, Type)


protein <- dataset.protein$data$norm %>% 
  data.frame() %>% 
  rownames_to_column("mouse.id") %>% 
  gather("protein.id","Expression",-mouse.id) %>% 
  left_join(dataset.protein$annot.protein %>% 
              select(protein.id, symbol)
  ) %>% 
  filter(symbol %in% c("Psmb8","Psmb9","Psmb3","Psmd7")) %>%
  left_join(dataset.protein$annot.samples) %>% 
  mutate(Sex = ifelse(Sex == "F","Female","Male"),
         symbol = factor(symbol, levels = c("Psmb8","Psmb9","Psmb3","Psmd7")),
         Type = "Protein") %>% 
  select(mouse.id, symbol, Sex, Age, Expression, Type)

df <- rbind(
  transcript,
  protein
) %>% 
  mutate(Age = factor(Age, levels = c(6, 12, 18)),
         Type = factor(Type, levels = c("Transcript","Protein")))

# Transcripts

plot1 <- df %>% 
  filter(symbol %in% c("Psmb8","Psmb9"),
         Type == "Transcript") %>% 
  ggplot(aes(x = Age, y = Expression)) +
  geom_boxplot(outlier.shape = NA, aes(color = Age, fill = Age), alpha = 0.3) +
  geom_point(aes(color = Age), 
             position = position_jitter(height = 0, width = 0.1), 
             alpha = 0.5) +
  geom_smooth(aes(x = as.integer(Age), y = Expression), method="lm",
              color="black", linetype="dashed") +
  labs(x = "Age", y = "Normalized expression") +
  theme_bw() +
  facet_wrap(~ symbol, scales = "free") +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  scale_fill_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  ggtitle("Psmb8 and Psmb9 - Transcript")

plot2 <- df %>% 
  filter(symbol %in% c("Psmb3","Psmd7"),
         Type == "Transcript") %>% 
  ggplot(aes(x = Age, y = Expression)) +
  geom_boxplot(outlier.shape = NA, aes(color = Age, fill = Age), alpha = 0.3) +
  geom_point(aes(color = Age), 
             position = position_jitter(height = 0, width = 0.1), 
             alpha = 0.5) +
  geom_smooth(aes(x = as.integer(Age), y = Expression), method="lm",
              color="black", linetype="dashed") +
  labs(x = "Age", y = "Normalized expression") +
  theme_bw() +
  facet_wrap(~ symbol, scales = "free") +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  scale_fill_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  ggtitle("Psmb3 and Psmd7 - Transcript")


subset <- df %>% 
  filter(symbol %in% c("Psmb8","Psmb9"),
         Type == "Transcript") %>% 
  spread("symbol","Expression")

plot3 <-  subset %>% 
  ggplot(aes(x = Psmb8, y = Psmb9, color = Age)) +
  geom_point(alpha = 0.3, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  ggtitle("")

subset <- df %>% 
  filter(symbol %in% c("Psmb3","Psmd7"),
         Type == "Transcript") %>% 
  spread("symbol","Expression")

plot4 <- subset %>% 
  ggplot(aes(x = Psmb3, y = Psmd7, color = Age)) +
  geom_point(alpha = 0.3, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  ggtitle("")

# Proteins

plot5 <- df %>% 
  filter(symbol %in% c("Psmb8","Psmb9"),
         Type == "Protein") %>% 
  ggplot(aes(x = Age, y = Expression)) +
  geom_boxplot(outlier.shape = NA, aes(color = Age, fill = Age), alpha = 0.3) +
  geom_point(aes(color = Age), 
             position = position_jitter(height = 0, width = 0.1), 
             alpha = 0.5) +
  geom_smooth(aes(x = as.integer(Age), y = Expression), method="lm",
              color="black", linetype="dashed") +
  labs(x = "Age", y = "Normalized expression") +
  theme_bw() +
  facet_wrap(~ symbol, scales = "free") +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  scale_fill_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  ggtitle("Psmb8 and Psmb9 - Protein")

plot6 <- df %>% 
  filter(symbol %in% c("Psmb3","Psmd7"),
         Type == "Protein") %>% 
  ggplot(aes(x = Age, y = Expression)) +
  geom_boxplot(outlier.shape = NA, aes(color = Age, fill = Age), alpha = 0.3) +
  geom_point(aes(color = Age), 
             position = position_jitter(height = 0, width = 0.1), 
             alpha = 0.5) +
  geom_smooth(aes(x = as.integer(Age), y = Expression), method="lm",
              color="black", linetype="dashed") +
  labs(x = "Age", y = "Normalized expression") +
  theme_bw() +
  facet_wrap(~ symbol, scales = "free") +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  scale_fill_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  ggtitle("Psmb3 and Psmd7 - Protein")


subset <- df %>% 
  filter(symbol %in% c("Psmb8","Psmb9"),
         Type == "Protein") %>% 
  spread("symbol","Expression")

plot7 <-  subset %>% 
  ggplot(aes(x = Psmb8, y = Psmb9, color = Age)) +
  geom_point(alpha = 0.3, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  ggtitle("")

subset <- df %>% 
  filter(symbol %in% c("Psmb3","Psmd7"),
         Type == "Protein") %>% 
  spread("symbol","Expression")

plot8 <- subset %>% 
  ggplot(aes(x = Psmb3, y = Psmd7, color = Age)) +
  geom_point(alpha = 0.3, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) +
  ggtitle("")

pdf("Downstream_Analysis/Results/proteincomplex_proteasome_cor.pdf",
    width = 10, height = 10)


cowplot::plot_grid(plot1,
                   plot3,
                   plot2,
                   plot4,
                   plot5,
                   plot7,
                   plot6,
                   plot8,
                   ncol = 2,
                   rel_widths = c(1.5,1))
dev.off()
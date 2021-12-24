#################### QTL mapping analysis on the DO heart ######################

# Generating a transcriptome map for the LOD scores across the genome.
# Adapting Duy Pham' script for that. 

################################################################################

# Author: Isabela Gerdes Gyuricza - Churchill lab 
# Date: 09-30-2021





################################################################################
#
#   This script is used to generate a 3-panel plot including transcriptome, 
# cis QTL, and distal QTL density plots.
#
#   Note*:
#       Install 'cowplots' package (install.packages('cowplot')) but do not load!
#
#
#   Input:
#       1: viewer_data:   Path to qtl viewer .RData 
#       2: dataset:       Which dataset to pull lod.peaks table
#       3: type_peak:     Which dataframe to pull from lod.peaks list
#       4: slide:         How much to slide across genome
#       5: window:        Window to count QTLs
#       6: lod_thres:     Threshold to filter LOD scores in lod.peaks table
#       7: density_thres: Y-axis value for horizontal line in density plots      
#       8: cis_color:     Color for local QTLs
#       9: dis_color:     Color for distal QTLs
#
#
#   Output: 
#       1: 3-panel plot showing transcriptome of LOD scores 
# (gene location by QTL location), 
#           distal QTL density, and cis QTL density
#
#
#
#
#   Authors: Duy Pham
#   Date:    April 9, 2019
#   E-mails: duy.pham@jax.org
################################################################################


################################################################################
############ clear workspace

rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)


### Options and Libraries
options(stringsAsFactors = FALSE)
library(tidyverse)
library(grid)
library(gridExtra)
library(plyr)
library(data.table)

# Loading dataset

load("QTLviewer_dataset/JAC_DO_heart_v9.gz")


pdf("Downstream_Analysis/QTL_mapping/Results/genome_scans.pdf",
    width = 14, height = 8)

################################ eQTLs #########################################

dataset     <- 'dataset.mrna'
#dataset     <- 'dataset.protein'
type_peak   <- 'age_int'
protein <- FALSE
#protein <- TRUE
slide  <- 1
window <- 4
lod.thres <- 7.75 #Using the lod thresholds from the permutations
density.thres <- 40
cis_color <- "gray30"
dis_color <- 'dodgerblue4'


dataset <- get(dataset)

if(protein == FALSE) {
  
  annots <- dataset$annot.mrna
  
  peaks <- dataset[['lod.peaks']][[type_peak]] %>%
    left_join(markers, by = 'marker.id') %>%  # Get marker ids
    dplyr::rename(qtl.chr   = chr,
                  qtl.pos   = pos) %>%
    dplyr::select(gene.id, marker.id, lod, qtl.chr, qtl.pos) %>% 
    left_join(annots %>% 
                dplyr::select(gene.id,chr,
                              start,end,symbol), 
              by='gene.id')
  
  
} else {
  
  annots <- dataset$annot.protein
  
  peaks <- dataset[['lod.peaks']][[type_peak]] %>%
    left_join(markers, by = 'marker.id') %>%  # Get marker ids
    dplyr::rename(qtl.chr   = chr,
                  qtl.pos   = pos) %>%
    dplyr::select(protein.id, marker.id, lod, qtl.chr, qtl.pos) %>% 
    left_join(annots %>% 
                dplyr::select(protein.id, gene.id,chr,
                              start,end,symbol), 
              by='protein.id')
  
}

peaks <- peaks %>% 
  dplyr::rename(gene.start=start,
                gene.end=end,
                gene.chr=chr,
                gene.symbol=symbol) %>% 
  mutate(cis = (gene.chr == qtl.chr) & abs(qtl.pos-gene.start) <= 4)


dataset[['lod.peaks']][[type_peak]] <- peaks


### Get data and filter lod.peaks table

lod.peaks    <- dataset[['lod.peaks']][[type_peak]]
lod.peaks    <- lod.peaks[lod.peaks$lod > lod.thres,]
lod.peaks    <- lod.peaks[complete.cases(lod.peaks),]
lod.peaks    <- lod.peaks[lod.peaks$gene.chr %in% c(1:19,'X'),]
local.peaks  <- lod.peaks[lod.peaks$cis,]
distal.peaks <- lod.peaks[!lod.peaks$cis,]
colnames(lod.peaks) <- gsub('.','_',colnames(lod.peaks), fixed = TRUE)



### Format chromosomes for transcriptome plot
all.chr = lod.peaks %>%
  select(qtl_chr, gene_chr) %>%
  gather(k, v) %>%
  select(v) %>%
  distinct() %>%
  arrange(v)

all.chr <- all.chr$v[!is.na(all.chr$v)]

if(length(grep("M", all.chr)) > 0){
  wh <- grep("M", all.chr)
  all.chr <- all.chr[c(1:(wh-1), (wh+1):length(all.chr), wh)]
}




### Format lod.peaks for transcriptome plot
data = lod.peaks %>% 
  mutate(cis      = (gene_chr == qtl_chr) & (abs(gene_start - qtl_pos) <= 4),
         qtl_chr  = factor(qtl_chr, levels   = all.chr[order(as.numeric(all.chr))]),
         gene_chr = factor(gene_chr,levels   = rev(all.chr[order(as.numeric(all.chr))])),
         gene_pos = (gene_end + gene_start) * 0.5)
cis.colors        = c(dis_color, cis_color)
names(cis.colors) = c("FALSE", "TRUE")


## Cis LOD counts
lod_df_cis <- list()
for(i in unique(markers$chr)){
  
  # Finding floor of minimum marker position and ceiling of maximum marker position
  min <- round_any(min(map[[i]]), 1, f = floor)
  max <- round_any(max(map[[i]]), 4, f = ceiling)
  
  # Creating x-axis scale. min to max with slide (or 'by' in seq function)
  x_axis_scale <- seq(min, max, slide)
  chr <- rep(i, length(x_axis_scale))
  
  # Getting LOD peaks from chromosome i
  sub <- subset(local.peaks, qtl.chr == i)
  
  
  # Creating dataframe of counts of lod peaks above threshold 
  count <- vector()
  pos <- ((x_axis_scale+window)+x_axis_scale)/2
  for(j in 1:length(pos)){
    count[j] <- sum((abs(sub$qtl.pos - pos[j])  <= window/2))
  }
  
  lod_df_cis[[i]] <- data.frame(chr = chr, pos = pos, count = count)
}

lod_df_cis     <- rbindlist(lod_df_cis)
lod_df_cis$chr <- factor(lod_df_cis$chr, levels = c(1:19,'X'))



### Distal LOD counts
lod_df_dis <- list()
for(i in unique(markers$chr)){
  
  # Finding floor of minimum marker position and ceiling of maximum marker position
  min <- round_any(min(map[[i]]), 1, f = floor)
  max <- round_any(max(map[[i]]), 4, f = ceiling)
  
  # Creating x-axis scale. min to max with slide (or 'by' in seq function)
  x_axis_scale <- seq(min, max, slide)
  chr <- rep(i, length(x_axis_scale))
  
  # Getting LOD peaks from chromosome i
  sub <- subset(distal.peaks, qtl.chr == i)
  
  
  # Creating dataframe of counts of lod peaks above threshold 
  count <- vector()
  pos <- ((x_axis_scale+window)+x_axis_scale)/2
  for(j in 1:length(pos)){
    count[j] <- sum((abs(sub$qtl.pos - pos[j])  <= window/2))
  }
  
  lod_df_dis[[i]] <- data.frame(chr = chr, pos = pos, count = count)
}

lod_df_dis <- rbindlist(lod_df_dis)
lod_df_dis$chr <- factor(lod_df_dis$chr, levels = c(1:19,'X'))

### Plot transcriptome
g_rna <- ggplot(data, aes(x = qtl_pos, y = gene_pos, size = log(lod)), alpha = 0.4) +
  geom_point(aes(color = cis), alpha = 0.4) + 
  scale_color_manual(values = cis.colors, labels = c('Distal-pQTL','Local-pQTL')) +
  facet_grid(gene_chr ~ qtl_chr, scales = "free",switch = 'y', shrink = TRUE) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = 0, color = "black"),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.0, "line"),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(angle = 180),
        axis.text = element_blank(),
        axis.line.y.left = element_blank(),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_blank(),
        legend.key = element_blank()) +
  xlab('eQTL position (Chr)') +
  ylab('Gene position (Chr)')

### Plot distal density plot
d_rna <- ggplot(lod_df_dis, aes(x = pos, y = count)) +
  geom_line(col = dis_color, size = 1) +
  #geom_hline(yintercept = density.thres, linetype = "dotted", color = 'black', size = 1) +
  ylab(paste('No. distal eQTL/',window,'Mbp')) +
  xlab('Chromosome') +
  theme(panel.spacing.x = unit(0, "lines"),
        plot.margin = unit(c(0,0,1.3,0),"cm"),
        panel.background = element_blank(),
        panel.border = element_rect(fill = 0, color = "black"),
        axis.text.x = element_blank(),
        strip.text = element_text(size= 23),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12)) + 
  scale_y_continuous(limits = c(0,45), breaks = seq(0,45,15)) +  
  #I first ran without this line and then adjusted this values based on each plot.
  facet_grid(.~chr, scales = "free",switch = 'x')

################################ pQTLs #########################################

#dataset     <- 'dataset.mrna'
dataset     <- 'dataset.protein'
type_peak   <- 'age_int'
#protein <- FALSE
protein <- TRUE
slide  <- 1
window <- 4
lod.thres <- 7.75 #Using the lod thresholds from the permutations
density.thres <- 40
cis_color <- "gray30"
dis_color <- 'dodgerblue4'


dataset <- get(dataset)

if(protein == FALSE) {
  
  annots <- dataset$annot.mrna
  
  peaks <- dataset[['lod.peaks']][[type_peak]] %>%
    left_join(markers, by = 'marker.id') %>%  # Get marker ids
    dplyr::rename(qtl.chr   = chr,
                  qtl.pos   = pos) %>%
    dplyr::select(gene.id, marker.id, lod, qtl.chr, qtl.pos) %>% 
    left_join(annots %>% 
                dplyr::select(gene.id,chr,
                              start,end,symbol), 
              by='gene.id')
  
  
} else {
  
  annots <- dataset$annot.protein
  
  peaks <- dataset[['lod.peaks']][[type_peak]] %>%
    left_join(markers, by = 'marker.id') %>%  # Get marker ids
    dplyr::rename(qtl.chr   = chr,
                  qtl.pos   = pos) %>%
    dplyr::select(protein.id, marker.id, lod, qtl.chr, qtl.pos) %>% 
    left_join(annots %>% 
                dplyr::select(protein.id, gene.id,chr,
                              start,end,symbol), 
              by='protein.id')
  
}

peaks <- peaks %>% 
  dplyr::rename(gene.start=start,
                gene.end=end,
                gene.chr=chr,
                gene.symbol=symbol) %>% 
  mutate(cis = (gene.chr == qtl.chr) & abs(qtl.pos-gene.start) <= 4)


dataset[['lod.peaks']][[type_peak]] <- peaks


### Get data and filter lod.peaks table

lod.peaks    <- dataset[['lod.peaks']][[type_peak]]
lod.peaks    <- lod.peaks[lod.peaks$lod > lod.thres,]
lod.peaks    <- lod.peaks[complete.cases(lod.peaks),]
lod.peaks    <- lod.peaks[lod.peaks$gene.chr %in% c(1:19,'X'),]
local.peaks  <- lod.peaks[lod.peaks$cis,]
distal.peaks <- lod.peaks[!lod.peaks$cis,]
colnames(lod.peaks) <- gsub('.','_',colnames(lod.peaks), fixed = TRUE)



### Format chromosomes for transcriptome plot
all.chr = lod.peaks %>%
  select(qtl_chr, gene_chr) %>%
  gather(k, v) %>%
  select(v) %>%
  distinct() %>%
  arrange(v)

all.chr <- all.chr$v[!is.na(all.chr$v)]

if(length(grep("M", all.chr)) > 0){
  wh <- grep("M", all.chr)
  all.chr <- all.chr[c(1:(wh-1), (wh+1):length(all.chr), wh)]
}




### Format lod.peaks for transcriptome plot
data = lod.peaks %>% 
  mutate(cis      = (gene_chr == qtl_chr) & (abs(gene_start - qtl_pos) <= 4),
         qtl_chr  = factor(qtl_chr, levels   = all.chr[order(as.numeric(all.chr))]),
         gene_chr = factor(gene_chr,levels   = rev(all.chr[order(as.numeric(all.chr))])),
         gene_pos = (gene_end + gene_start) * 0.5)
cis.colors        = c(dis_color, cis_color)
names(cis.colors) = c("FALSE", "TRUE")




## Cis LOD counts
lod_df_cis <- list()
for(i in unique(markers$chr)){
  
  # Finding floor of minimum marker position and ceiling of maximum marker position
  min <- round_any(min(map[[i]]), 1, f = floor)
  max <- round_any(max(map[[i]]), 4, f = ceiling)
  
  # Creating x-axis scale. min to max with slide (or 'by' in seq function)
  x_axis_scale <- seq(min, max, slide)
  chr <- rep(i, length(x_axis_scale))
  
  # Getting LOD peaks from chromosome i
  sub <- subset(local.peaks, qtl.chr == i)
  
  
  # Creating dataframe of counts of lod peaks above threshold 
  count <- vector()
  pos <- ((x_axis_scale+window)+x_axis_scale)/2
  for(j in 1:length(pos)){
    count[j] <- sum((abs(sub$qtl.pos - pos[j])  <= window/2))
  }
  
  lod_df_cis[[i]] <- data.frame(chr = chr, pos = pos, count = count)
}

lod_df_cis     <- rbindlist(lod_df_cis)
lod_df_cis$chr <- factor(lod_df_cis$chr, levels = c(1:19,'X'))



### Distal LOD counts
lod_df_dis <- list()
for(i in unique(markers$chr)){
  
  # Finding floor of minimum marker position and ceiling of maximum marker position
  min <- round_any(min(map[[i]]), 1, f = floor)
  max <- round_any(max(map[[i]]), 4, f = ceiling)
  
  # Creating x-axis scale. min to max with slide (or 'by' in seq function)
  x_axis_scale <- seq(min, max, slide)
  chr <- rep(i, length(x_axis_scale))
  
  # Getting LOD peaks from chromosome i
  sub <- subset(distal.peaks, qtl.chr == i)
  
  
  # Creating dataframe of counts of lod peaks above threshold 
  count <- vector()
  pos <- ((x_axis_scale+window)+x_axis_scale)/2
  for(j in 1:length(pos)){
    count[j] <- sum((abs(sub$qtl.pos - pos[j])  <= window/2))
  }
  
  lod_df_dis[[i]] <- data.frame(chr = chr, pos = pos, count = count)
}

lod_df_dis <- rbindlist(lod_df_dis)
lod_df_dis$chr <- factor(lod_df_dis$chr, levels = c(1:19,'X'))

### Plot transcriptome
g_protein <- ggplot(data, aes(x = qtl_pos, y = gene_pos, size = log(lod)), alpha = 0.4) +
  geom_point(aes(color = cis), alpha = 0.4) + 
  scale_color_manual(values = cis.colors, labels = c('Distal-pQTL','Local-pQTL')) +
  facet_grid(gene_chr ~ qtl_chr, scales = "free",switch = 'y', shrink = TRUE) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = 0, color = "black"),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.0, "line"),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(angle = 180),
        axis.text = element_blank(),
        axis.line.y.left = element_blank(),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_blank(),
        legend.key = element_blank()) +
  xlab('pQTL position (Chr)') +
  ylab('Gene position (Chr)')

### Plot distal density plot
d_protein <- ggplot(lod_df_dis, aes(x = pos, y = count)) +
  geom_line(col = dis_color, size = 1) +
  #geom_hline(yintercept = density.thres, linetype = "dotted", color = 'black', size = 1) +
  ylab(paste('No. distal pQTL/',window,'Mbp')) +
  xlab('Chromosome') +
  theme(panel.spacing.x = unit(0, "lines"),
        plot.margin = unit(c(0,0,1.3,0),"cm"),
        panel.background = element_blank(),
        panel.border = element_rect(fill = 0, color = "black"),
        axis.text.x = element_blank(),
        strip.text = element_text(size= 23),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12)) + 
  scale_y_continuous(limits = c(0,45), breaks = seq(0,45,15)) + 
  #I first ran without this line and then adjusted this values based on each plot.
  facet_grid(.~chr, scales = "free",switch = 'x')



### Put plots in 3 panels vertically. DONT LOAD COWPLOT LIBRARY!!
cowplot::plot_grid(g_rna, g_protein, d_rna, d_protein,
                   axis = 'lr', ncol = 2, align = 'v', 
                   rel_heights = c(1,.3,.3), 
                   vjust = c(1.5,6.0,.85))

dev.off()

################################################################################

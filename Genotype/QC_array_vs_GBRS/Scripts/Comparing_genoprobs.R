############################ Genotypes QC ######################################

# Comparing array genoprobs with GBRS genoprobs to look for samples mixups.

################################################################################

# Author: Isabela Gerdes Gyuricza - Churchill lab 
# Date: 06_27_2021

################################################################################
############ clear workspace

rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

################################################################################
############ loading libraries

library(qtl2)
library(tidyverse)

################################################################################
############ loading data

setwd("~/Box/JAC_Heart_Data/Heart_Data_June2021")

load("Genotype/Data/GBRS_genoprobs_gathered.RData")
# Getting genoprobs with 192 mice

load("Genotype/Data/GBRS_map.RData")
# map data

genoprobs_array <- readRDS("Genotype/Data/CS_AlleleProbs.rds")
# array genoprobs

cross2 <- readRDS("Genotype/Data/QTL2_cross2_LongAndCS.rds")
# object that contains array map

array_map <- cross2$pmap

rm(cross2)

# Array genoprobs and map have Y, MT and pseudo-autosomal regions. Removing those 

genoprobs_array$Y <- NULL

genoprobs_array$M <- NULL

genoprobs_array$P <- NULL


array_map$Y <- NULL

array_map$M <- NULL

array_map$P <- NULL

# Array genoprobs is not interpolated...

array_map <- insert_pseudomarkers(array_map, pseudomarker_map = map)

genoprobs_array_interp <- interp_genoprob(genoprobs_array, array_map)

genoprobs_array_interp <- interp_genoprob(genoprobs_array_interp, map)

# Now it is!

rm(genoprobs_array, array_map)

genoprobs_array <- genoprobs_array_interp

rm(genoprobs_array_interp)

# Renaming GBRS genoprobs to make it easier

genoprobs <- genoprobs_GBRS


# Keeping same mouse.id pattern

# Vector of names you want to change to 

vector <- rownames(genoprobs_array[[1]])

vector <- gsub("DO-","DO.",vector)

# Names are the original IDs.

names(vector) <- rownames(genoprobs_array[[1]])

# Applying the qtl2 function..

genoprobs_array <- replace_ids(genoprobs_array,vector)

rm(vector)

# No need to subset only common mice, functions below already do that

# Confirming that array and GBRS genoprobs have same number of markers

for (i in 1:length(genoprobs_array)) {
  
  print(dim(genoprobs_array[[i]]) == dim(genoprobs[[i]]))
  print(i)
  
} #True for haplotypes and markers

rm(i)


# Using Greg's function for comparing genoprobs

source("Genotype/QC_array_vs_GBRS/Scripts/compare_genoprobs_functions.R")

dist_genoprobs <- compare_dist_genoprobs(genoprobs1 = genoprobs_array, 
                                         genoprobs2 = genoprobs,
                                         map_list = map, 
                                         use_dosages = FALSE)


mean_dist <- rowMeans(dist_genoprobs$dist_mat)


# Plotting all them togehter..

pdf("Genotype/QC_array_vs_GBRS/histogram_distances.pdf",
    width = 10, height = 8)

hist(mean_dist, main = "GBRS genoprobs vs Array genoprobs")

dev.off()

# There are some high values.. 

mean_dist[mean_dist > 0.8]
# DO.1077   DO.1157 
# 0.8171449 0.8624699 

# Apparently, mice DO.1077 and DO.1157 are swapped in the RNAseq..

########## Trying to fix them. 

vector <- rownames(genoprobs[[1]])

names(vector) <- rownames(genoprobs[[1]])

vector[which(names(vector) == "DO.1077")] <- "DO.1157"

vector[which(names(vector) == "DO.1157")] <- "DO.1077"

genoprobs_swap <- replace_ids(genoprobs,vector)


# Comparison again:

dist_genoprobs <- compare_dist_genoprobs(genoprobs1 = genoprobs_array, 
                                         genoprobs2 = genoprobs_swap,
                                         map_list = map, 
                                         use_dosages = FALSE)


mean_dist <- rowMeans(dist_genoprobs$dist_mat)

# Plotting all them togehter..

pdf("Genotype/QC_array_vs_GBRS/histogram_distances_swap.pdf",
    width = 10, height = 8)

hist(mean_dist, main = "swapped GBRS genoprobs vs Array genoprobs")

dev.off()

# There is still something there

mean_dist[mean_dist > 0.8]
# DO.1157   DO.1077 
# 0.8639815 0.8540533 

#They are not swapped. Are they duplications of some other mice?

# Using Greg's function

dup_geno_dist <- check_duplicate_genoprobs(genoprobs,"DO.1157",map_list = map)

mean <- rowMeans(dup_geno_dist$dist_mat)

mean[mean < 0.2]
#DO.1157 
#0

# DO.1157 is only identical to itself


dup_geno_dist <- check_duplicate_genoprobs(genoprobs,"DO.1077",map_list = map)

mean <- rowMeans(dup_geno_dist$dist_mat)

mean[mean < 0.2]
# DO.0938    DO.1077 
# 0.09469223 0.00000000 

# DO.1077 might be a duplication of DO.0938 in the RNA-seq, but DO.1157 doesn't
# match anything. Try excluding DO.1077 and DO.1157 and see if the comparison looks good.


# Creating the vector to exclude mice

vector <- rownames(genoprobs[[1]])

 names(vector) <- rownames(genoprobs[[1]])
 
vector <- vector[-which(names(vector) %in% c("DO.1157","DO.1077"))]

genoprobs <- replace_ids(genoprobs,vector)

dist_genoprobs <- compare_dist_genoprobs(genoprobs1 = genoprobs_array, 
                                         genoprobs2 = genoprobs,
                                         map_list = map, 
                                         use_dosages = FALSE)

mean_dist <- rowMeans(dist_genoprobs$dist_mat)

# Plotting it

pdf("Genotype/QC_array_vs_GBRS/histogram_distances_removing.pdf",
    idth = 10, height = 8)

hist(mean_dist, main = "GBRS genoprobs vs Array genoprobs")

dev.off()


mean_dist[mean_dist > 0.8]
# named numeric(0)

# Nice! Keeping these two out on the array genoprobs

vector <- rownames(genoprobs_array[[1]])

names(vector) <- rownames(genoprobs_array[[1]])

vector <- vector[which(names(vector) %in% rownames(genoprobs[[1]]))]

genoprobs_array <- replace_ids(genoprobs_array,vector)

rownames(genoprobs_array[[1]]) %>% length()
# Now we have only 187 mice. Saving new genoprobs

saveRDS(genoprobs_array, file = "Genotype/QC_array_vs_GBRS/Results/genoprobs.RDS")

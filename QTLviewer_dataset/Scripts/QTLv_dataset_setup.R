########################## QTLviewer dataset ###################################

# Generating QTLviewer with genotype data (genoprobs, K, markers and map) and 
# expression matrices (rna and protein).

################################################################################

# Author: Isabela Gerdes Gyuricza - Churchill lab 
# Date: 06_30_2021

################################################################################
############ clear workspace

rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

################################################################################
############ loading libraries

library(ensimplR)
library(qtl2)
library(tidyverse)

################################################################################
############ loading data

setwd("~/Box/JAC_Heart_Data/Heart_Data_June2021")

genoprobs <- readRDS("Genotype/QC_array_vs_GBRS/Results/genoprobs.RDS")
# Loading new array genoprobs

load("Genotype/Data/GBRS_map.RData")
# Loading map data

counts <- readRDS("RNA/Results/rna_counts.RDS")
# Loading raw counts

protein_data <- read_csv("Protein/Results/debatch_nopoly_protein_heart.csv")
# Loading processed protein data

annot.samples <- protein_data %>%
  dplyr::rename(Batch_MS = Batch,
                Tag_MS = tag) %>%
  dplyr::select(mouse.id, Sex, Age, Generation, Batch_MS, Tag_MS) %>% 
  unique() %>%   
  column_to_rownames("mouse.id")
# Loading samples annotation that will be used for both rna and protein datasets


# First, grabbing only common mice between genoprobs array, RNA and protein.

mice <- intersect(
  intersect(
            rownames(genoprobs[[1]]),
                  rownames(counts)
    ),
    protein_data$mouse.id %>% unique()
  )

# 185 mice in common 



###################### Generating genotype dataset #############################
################################################################################

#Fisrt, putting genoprobs in the same order as "mice"

for (i in 1:length(genoprobs)) {

  genoprobs[[i]] <- genoprobs[[i]][mice,,]

}

rm(i)

# Transforming genoprobs into the qtl2 object

class(genoprobs) <- c("calc_genoprob","list")

# Generating Kinship matrix using qlt2

K <- calc_kinship(genoprobs, type = "loco", cores = 2)

identical(rownames(K$`1`), rownames(genoprobs[[1]]))
#[1] TRUE

# Generating markers variable using qtl2

markers <- qtl2convert::map_list_to_df(map)  

# Adding ENSEMBL version for QTLviewer

ensembl.version <- 84

# Saving genotype QTLviewer dataset

 save(genoprobs, K, map, markers,ensembl.version, 
        file = "QTLviewer_dataset/Data/Genotype/dataset.genotype.RData")

################################################################################
################################################################################



################################################################################
# Setting up common variables, such as covar.matrix, covar.info and annot.samples

annot.samples <- annot.samples[mice,] %>% 
  rownames_to_column("mouse.id") %>% 
  mutate(Age = as.numeric(
    as.character(Age)
  )
  ) %>% 
  column_to_rownames("mouse.id")

covar.matrix <- model.matrix(~ Sex + Age, data = annot.samples)[,-1]

class(covar.matrix) #matrix, done. 

annot.samples <- annot.samples %>% 
  rownames_to_column("mouse.id")

covar.info <- tibble(sample.column = c("Sex","Age"),
                     covar.column = c("Sex","Age"),
                     display.name = c("Sex","Age"),
                     interactive = as.factor(c("TRUE","TRUE")),
                     primary = as.factor(c("TRUE","TRUE")),
                     lod.peaks = c("sex_int","age_int"))

covar.info <- covar.info %>% 
  mutate(interactive = as.logical(interactive),
         primary = as.logical(primary))

################################################################################
################################################################################



######################## Generating rna dataset ################################
################################################################################

# Filtering RNA dataset keeping only genes that have more than 1 read for at 
# least half of the samples and a median of at least 1.

# Computing the median per gene
med <- apply(counts,2,median)

# Computing number of samples with reads >=1
samples <- apply(counts >= 1, 2, sum)  

# index for subsetting genes
indx.keep <- which(samples >= 95 & med >= 1)

counts <- counts[mice, indx.keep] 
#reordering mice to match annots and keeping only the filtered genes

dim(counts)
# 185 x 21016

##Check if the Cdk2na is in the dataset! 

counts[,"ENSMUSG00000044303"] # - great! It is there :)

rm(med, samples, indx.keep)

# Generating data variable with raw, vst and rankzed expression

norm <- t(
  DESeq2::vst(
  t(
  floor(
    counts
  )
)
)
)

# Setting up RANKZ function:
rankZ <- function(x) {
  x <- rank(x, na.last = "keep") / (length(x) - sum(is.na(x)) + 1)
  return(qnorm(x))
}


data <- list(raw = counts, 
             norm = norm,
             rz = apply(
               norm, 2, rankZ
               )
)

rm(norm,counts)

# Creating annot.mrna variable

annot.mrna <- batchGenes(colnames(data$raw), release = 84)

entrez_id <- str_split(annot.mrna$external_ids,",",simplify = TRUE)[,1]

entrez_id <- regmatches(entrez_id, gregexpr("[[:digit:]]+", entrez_id))

names(entrez_id) <- annot.mrna$id

test_list <- function(x) {
  if (length(x) > 1 | length(x) == 0) {
    x <- "NA"
  } else {
    x
  }
}

entrez_id <- lapply(entrez_id, test_list)

entrez_id <- bind_rows(entrez_id)

annot.mrna <- annot.mrna %>% 
  select(id,symbol,chromosome,start,end,strand) %>% 
  mutate(entrez_id = as.character(entrez_id),
         middle = (((end-start)/2)+start)/1000000,
         start = start/1000000,
         end = end/1000000) %>% 
  rename(gene.id = id,
         chr = chromosome) %>% 
  select(gene.id, symbol, entrez_id,chr,start,end,middle,strand)


nearest_marker <- character(length = nrow(annot.mrna))
for(i in 1:nrow(annot.mrna)){
  cm <- annot.mrna$chr[i]
  if (cm %in% c(1:19,"X")){
    sub <- subset(markers, chr == cm)
    nearest_marker[i] <- sub$marker[which.min(abs(sub$pos - annot.mrna$start[i]))]
    
  }};rm(i,sub,cm)

annot.mrna$nearest.marker.id <- nearest_marker 


assign("dataset.mrna",
       list(annot.mrna = as_tibble(annot.mrna),
            covar.matrix = covar.matrix,
            covar.info = covar.info,
            datatype = "mrna",
            display.name = "DO heart transcriptomics",
            data = data,
            annot.samples = as_tibble(annot.samples)))

saveRDS(dataset.mrna, file = "QTLviewer_dataset/Data/RNA/dataset.mrna.RDS")

################################################################################
################################################################################


###################### Generating protein dataset ##############################
################################################################################

protein_data <- protein_data %>%
  spread(key = protein.id, value = Intensity) %>%
  column_to_rownames("mouse.id") %>%
  select(contains("ENSMUS")) %>%
  as.matrix()

# Filtering to keep only the common mice
protein_data <- protein_data[mice,]

dim(protein_data)
#[1]  185 4223

data <- list(# No raw data for protein 
             norm = protein_data,
             rz = apply(
               protein_data, 2, rankZ
             )
)

rm(protein_data)

# Creating annot.mrna variable

annot.protein <- batchGenes(colnames(data$norm), release = 75, species = "Mm")
#We coudn't find annotations for two proteins ("ENSMUSP00000052262" "ENSMUSP00000111624")
# Keeping these out.

entrez_id <- str_split(annot.protein$external_ids,",",simplify = TRUE)[,1]

entrez_id <- regmatches(entrez_id, gregexpr("[[:digit:]]+", entrez_id))

names(entrez_id) <- annot.protein$id

test_list <- function(x) {
  if (length(x) > 1 | length(x) == 0) {
    x <- "NA"
  } else {
    x
  }
}

entrez_id <- lapply(entrez_id, test_list)

entrez_id <- bind_rows(entrez_id)

annot.protein <- annot.protein %>% 
  select(id,gene_id,symbol,chromosome,start,end,strand) %>% 
  mutate(entrez_id = as.character(entrez_id),
         middle = (((end-start)/2)+start)/1000000,
         start = start/1000000,
         end = end/1000000) %>% 
  rename(protein.id = id,
         gene.id = gene_id,
         chr = chromosome) %>% 
  select(protein.id, gene.id, symbol, entrez_id,chr,start,end,middle,strand)


nearest_marker <- character(length = nrow(annot.protein))
for(i in 1:nrow(annot.protein)){
  cm <- annot.protein$chr[i]
  if (cm %in% c(1:19,"X")){
    sub <- subset(markers, chr == cm)
    nearest_marker[i] <- sub$marker[which.min(abs(sub$pos - annot.protein$start[i]))]
    
  }};rm(i,sub,cm)

annot.protein$nearest.marker.id <- nearest_marker 

data$norm <- data$norm[,annot.protein$protein.id]

data$rz <- data$rz[,annot.protein$protein.id]


assign("dataset.protein",
       list(annot.protein = as_tibble(annot.protein),
            covar.matrix = covar.matrix,
            covar.info = covar.info,
            datatype = "protein",
            display.name = "DO heart proteomics",
            data = data,
            annot.samples = as_tibble(annot.samples)))

saveRDS(dataset.protein, file = "QTLviewer_dataset/Data/Protein/dataset.protein.RDS")



###################### Generating combined dataset #############################
################################################################################

# Clear workspace
rm(list = ls())

# Loading genotype, rna and protein data
load("QTLviewer_dataset/Data/Genotype/dataset.genotype.RData")

dataset.mrna <- readRDS("QTLviewer_dataset/Data/RNA/dataset.mrna.RDS")

dataset.protein <- readRDS("QTLviewer_dataset/Data/Protein/dataset.protein.RDS")

# Creating temporary QTLvfile: 

#save(list = ls(), file = gzfile("~/Desktop/JAC_DO_heart_v8.gz"))


# Everything looks fine =)

################################################################################

# NOTE: These datasets were used as input for the QTLviewer software, which 
# generates lod peaks for the QTL mapping. These files were then combined, 
# generating a RData object (used as input data for
# all of our downstream analysis).
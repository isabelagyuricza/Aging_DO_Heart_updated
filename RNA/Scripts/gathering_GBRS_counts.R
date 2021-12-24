########################## Gathering GBRS counts ###############################

# Gathering multiple files from GBRS counts for each sample into one large 
# counts matrix. Also, excluding the two problematic samples (DO.1077 and DO.1157)

################################################################################

# Author: Isabela Gerdes Gyuricza - Churchill lab 
# Date: 06_29_2021

################################################################################
############ clear workspace

rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

################################################################################
############ loading libraries

library(tidyverse)

# Loading counts data from GBRS

setwd("~/Box/JAC_Heart_Data/Heart_Data_June2021")

file_list <- list.files("RNA/Data", 
                        pattern="*.txt")

# create list of all .gene.counts.txt files in folder

list_counts <- list()
for (i in file_list){
  x <- read.table(
    paste(
      "RNA/Data/", 
      i, 
      sep=""),
    header=TRUE,
    sep = "\t")
  
  colnames(x) <- c("Target_ID","A","B","C","D","E","F","G","H","Total","Notes")
  
  list_counts[[i]] <- x
  
}; rm(x,i)
# read in each .txt file in file_list and create a data frame with the same name as the file

# Exclude samples "1077" and "1157"

list_counts[grepl("_1077", names(list_counts))] <- NULL

list_counts[grepl("_1157", names(list_counts))] <- NULL

#Checking if all the genes for all the samples are the same

ls <- list()

for (i in 1:length(list_counts)){
  
  ls <- list_counts[[i]]$Target_ID
  
}
#Factor with 47643 levels = taking all the genes for all the samples there are 
# 47643 levels (genes). There is no difference of the set of the genes among samples.

# Creating input data with ALL the genes (list_counts)

df_counts <- data.frame(
  row.names=list_counts[[1]]$Target_ID,
                        stringsAsFactors = FALSE) 
#I took the genes from the first sample, because they are all the same.

for (i in names(list_counts)){
  
  x <- str_split(i,pattern = "-",simplify = TRUE)
  
  x <- x[1]
  
  df_counts[[x]] <- list_counts[[i]]$Total
  
}

# Changing the colnames

for (i in 1:length(colnames(df_counts))) {
  
  x <- str_split(colnames(df_counts)[i],"_",simplify = TRUE)
  
  x <- x[2]
  
  x <- as.character(x)
  
  x <- as.numeric(x)
  
  x <- formatC(x,width = 4,format = "d",flag = "0")
  
  colnames(df_counts)[i] <- x
  
}

for (i in 1:length(colnames(df_counts))){
  
  y <- colnames(df_counts)[i]
  
  x <- paste0("DO.",y)
  
  colnames(df_counts)[i] <- x
  
}

rna_counts <- t(df_counts)

saveRDS(rna_counts,file = "RNA/Results/rna_counts.RDS")

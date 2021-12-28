#########################################
# JAC Heart
#  compare heart and kidney DE results. 
#
#########################################

#########################################
# Gary A Churchill, 2021

# setwd("~/Work/Projects/Aging_Center_2021/Aging_Heart_v3")


####
library(fgsea)
library(org.Mm.eg.db)
library(tidyverse)
library(stringr)
library(cowplot)
library(readxl)
library(modelr)
library(broom)
library(ggbeeswarm)

#########################################
# load DEA Results for heart and kidney

#####
# load heart DE results
# Note that this requires the Supplemental Data files

Heart_RNATest <- read_csv("Downstream_Analysis/Results/DEA_age_transcripts.csv") %>%
  mutate(Gene_symbol = str_sub(Gene_symbol, 2, -2)) %>%
  dplyr::rename(gene.id = Gene_ID, 
                # symbol = Gene_symbol,
                b.age.mrna = Effect,
                se.age.mrna = Effect_SE,
                p.age.mrna = p,
                q.age.mrna = adjusted_p) %>%
  mutate(z.age.mrna = b.age.mrna/se.age.mrna)

Heart_ProteinTest <- read_csv("Downstream_Analysis/Results/DEA_age_proteins.csv") %>%
  mutate(Gene_symbol = str_sub(Gene_symbol, 2, -2)) %>%
  dplyr::rename(protein.id = Protein_ID,
                gene.id = Gene_ID, 
                symbol = Gene_symbol,
                b.age.prot = Effect,
                se.age.prot = Effect_SE,
                p.age.prot = p,
                q.age.prot = adjusted_p) %>%
  mutate(z.age.prot = b.age.prot/se.age.prot)
#
common.genes <- intersect(Heart_RNATest$gene.id, Heart_ProteinTest$gene.id)
Heart_GeneTest <- left_join(filter(Heart_ProteinTest, gene.id %in% common.genes),
          filter(Heart_RNATest, gene.id %in% common.genes), by="gene.id")

#####
# load kidney DE result - from Takemon et al. 
Kidney_GeneTest <- read_csv("Downstream_Analysis/Data/Kidney_GeneTest.csv")


#####
# how big is the common set after merging RNA and proteins?  
dim(Heart_GeneTest)    #4117
length(unique(Heart_GeneTest$gene.id))
dim(Kidney_GeneTest)   #6514
length(unique(Kidney_GeneTest$gene.id))
#
common.genes <- unique(intersect(Kidney_GeneTest$protein.id, 
                                 Heart_GeneTest$protein.id))
length(common.genes)  
# there are 3514 genes in common across the two merged sets


#######
# create a merged set of proteins and transcripts that are in both tissues
Heart_GeneTest <- filter(Heart_GeneTest, protein.id %in% common.genes)
Kidney_GeneTest <- filter(Kidney_GeneTest, protein.id %in% common.genes)

####
# select variables
tmp.heart <- select(Heart_GeneTest,
                         protein.id, gene.id, symbol, 
                         z.age.mrna, q.age.mrna, z.age.prot, q.age.prot)  %>%
  dplyr::rename(z.heart.mrna = z.age.mrna,  
                q.heart.mrna = q.age.mrna, 
                z.heart.prot = z.age.prot, 
                q.heart.prot = q.age.prot)
#
tmp.kidney <- select(Kidney_GeneTest,
                          protein.id, # gene.id, symbol, 
                          z.age.mrna, q.age.mrna, z.age.prot, q.age.prot)  %>%
  dplyr::rename(z.kidney.mrna = z.age.mrna,  
                q.kidney.mrna = q.age.mrna, 
                z.kidney.prot = z.age.prot, 
                q.kidney.prot = q.age.prot)

####
# combine 
GeneTest <- left_join(tmp.heart, tmp.kidney)
#
length(unique(GeneTest$gene.id))
length(unique(GeneTest$protein.id))
# 3514 protein and 3499 genes

#####
# clean up
rm(common.genes, Kidney_GeneTest, Heart_GeneTest, tmp.heart, tmp.kidney)
ls()

#######################################################




#######################################################
# look at the significant DE genes in common to heart and kidney
#
# suggest to use 0.1 for RNA and 0.01 for protein
#######################################################


####
# number of DE genes in each tissue - includes duplicate gene.id counts
with(GeneTest, sum(q.heart.mrna < 0.01))  #  173
with(GeneTest, sum(q.kidney.mrna < 0.01)) #  268 
with(GeneTest, sum(q.heart.prot < 0.01))  # 1782
with(GeneTest, sum(q.kidney.prot < 0.01)) # 1414



####
# look at genes significant in both tissues
(sig.rna <- sort(filter(GeneTest, (q.heart.mrna < 0.1 & q.kidney.mrna < 0.1))$symbol))
#  45 rna in common at q <  0.01
# 245 at q < 0.1
#
(sig.prot <- sort(filter(GeneTest, (q.heart.prot < 0.01 & q.kidney.prot < 0.01))$symbol))
# 745 proteins in common at q < 0.01



####
# which genes are significantly changed in both tissues and for both RNA & Protien?
intersect(sig.rna, sig.prot)  # 19 genes
#

#####
# function print a list of gene.id's as symbols
GID2Symbol <- function(list){
  sort(unique(filter(GeneTest, (gene.id %in% list))$symbol))
}


############
# make gene lists for enrichments analysis

# ####
# # genes signficant in both rna and protein = 19
# mrna.protein <- unique(filter(GeneTest, 
#                               symbol %in% intersect(sig.rna, sig.prot))$gene.id)
# GID2Symbol(mrna.protein)
# 
# ######
# # all genes with significant rna changes = 47 genes short 
# mrna.all <- unique(filter(GeneTest, 
#                           symbol %in% sig.rna)$gene.id)
# GID2Symbol(mrna.all)
# 
# 
# ######
# # all genes with significant protein changes = 742
# prot.all <- unique(filter(GeneTest, 
#                           symbol %in% sig.prot)$gene.id)
# GID2Symbol(prot.all)



####
# mrna that go up in both = 25 unique genes
mrna.up.up <- unique(filter(GeneTest, (symbol %in% sig.rna) & 
                              (z.heart.mrna > 0) & (z.kidney.mrna > 0))$gene.id)
GID2Symbol(mrna.up.up)
#
# "C1qc"      "C4b"       "C6"        "Cfi"       "Cnn2"      "Cygb"     
# "Ddx39"     "Efhd2"     "Eno2"      "Ephx1"     "Gnpda2"    "Icam1"    
# "Igha"      "Igkc"      "Lcp1"      "Nudt16"    "Plek"      "Rpl3"     
# "S100a6"    "Serpina3n" "Sh3bgrl3"  "Sh3kbp1"   "Sncg"      "Vat1"     
# "Vcam1"   

####
# mrna that go down in both = 19
mrna.dn.dn <- unique(filter(GeneTest, (symbol %in% sig.rna) & 
                              (z.heart.mrna < 0) & (z.kidney.mrna < 0))$gene.id)
GID2Symbol(mrna.dn.dn)
# 
#  "Ahsa1"    "Akap1"    "Atp2a2"   "Atxn2"    "Chordc1"  "Cox15"   
# "Dnaja1"   "Erp44"    "Fkbp4"    "Hdhd2"    "Hsp90ab1" "Hspa1b"  
# "Hspd1"    "Hsph1"    "Lman1"    "Slc16a10" "St13"     "Stip1"   
# "Tmem135" 



####
# mrna that go up in heart & down in kidney = 1
mrna.up.dn <- unique(filter(GeneTest, (symbol %in% sig.rna) & 
                              (z.heart.mrna > 0) & (z.kidney.mrna < 0))$gene.id)
GID2Symbol(mrna.up.dn)
# "Rbp1"

####
# mrna that go down in heart & up in kidney,  = 0
mrna.dn.up <- unique(filter(GeneTest, (symbol %in% sig.rna) & 
                              (z.heart.mrna < 0) & (z.kidney.mrna > 0))$gene.id)
GID2Symbol(mrna.dn.up)
# none

####
# prot that go up in both = 319 unique genes
prot.up.up <- unique(filter(GeneTest, (symbol %in% sig.prot) & 
              (z.heart.prot > 0) & (z.kidney.prot > 0))$gene.id)
GID2Symbol(prot.up.up)

####
# prot that go down in both = 89
prot.dn.dn <- unique(filter(GeneTest, (symbol %in% sig.prot) & 
                              (z.heart.prot < 0) & (z.kidney.prot < 0))$gene.id)
GID2Symbol(prot.dn.dn)

####
# prot that go up in heart, down in kidney = 290
prot.up.dn <- unique(filter(GeneTest, (symbol %in% sig.prot) & 
                              (z.heart.prot > 0) & (z.kidney.prot < 0))$gene.id)
GID2Symbol(prot.up.dn)

####
# prot that go up in kidney, down in heart = 45
prot.dn.up <- unique(filter(GeneTest, (symbol %in% sig.prot) & 
                              (z.heart.prot < 0) & (z.kidney.prot > 0))$gene.id)
GID2Symbol(prot.dn.up)



#########################################
# plot the heart v. kidney age-trends for RNA & Proteins
# Suppl Fig 7B

####
# slope correlations
with(GeneTest, cor.test(z.heart.mrna, z.kidney.mrna))
# r = 0.450 p < 2.2x10-16
#
with(GeneTest, cor.test(z.heart.prot, z.kidney.prot))
# r = 0.027  p = 0.104

######
# scatter plot age slopes heart v kidney for RNA and for Protein

# plots std-slopes Heart v Kidney
p.mrna <- mutate(GeneTest, 
                 Sig_RNA = (q.heart.mrna < 0.1 & q.kidney.mrna < 0.1)) %>%
  arrange(Sig_RNA) %>%
  ggplot(aes(x=z.kidney.mrna, y=z.heart.mrna)) +
  geom_text(aes(label=symbol, color=Sig_RNA)) +
  scale_color_manual(values=c("darkgray", "black")) +
  guides(color = "none") +
  theme_bw(base_size = 18) + 
  annotate("text", x = -5, y = 7.5, label="r = 0.450", color="blue", size = 6) +
  annotate("text", x = -4.5, y = 7.0, label="p < 2.2x10-16", color="blue", size = 6) +
  # geom_smooth(method = "lm", color="red") +
  ggtitle("Kidney v. Heart mRNA")
#
p.prot <- mutate(GeneTest, 
                 Sig_Protein = (q.heart.prot < 0.01 & q.kidney.prot < 0.01)) %>%
  arrange(Sig_Protein) %>%
  ggplot(aes(x=z.kidney.prot, y=z.heart.prot)) +
  geom_text(aes(label=symbol, color=Sig_Protein)) +
  scale_color_manual(values=c("darkgray", "black")) +
  guides(color = "none") +
  theme_bw(base_size = 18) + 
  annotate("text", x = -9, y = 17, label="r = 0.027", color="blue", size = 6) +
  annotate("text", x = -8.9, y = 16, label="p = 0.104", color="blue", size = 6)  +
  # geom_smooth(method = "lm", color="red") +
  ggtitle("Kidney v. Heart Protein")

pdf("Downstream_Analysis/Results/HeartKidneyCompare_v2.pdf",height=10, width=20)
# quartz()
plot_grid(p.mrna, p.prot) 
          # labels = c("mRNA", "Protein"), label_size = 18)
dev.off()

#######################################################







#######################################################
#######################################################
#######################################################
# enrichments version 2 use ClusterProfiler
  
  
###
# load wrapper functions
source("Downstream_Analysis/Scripts/Enrichment_Functions.R")
# note: library loads this will mask some dplyr functions
# e.g. select() and simplify()



#################
# ENRICHMENT CODE BLOCK


####
# create a universe for heart-kidney comparison
all.geneids <- unique(GeneTest[[2]])


####
# cnet plots set up 
#
#for plotting fold change on color scale
# colpal <- c('#4dac26','#f1b6da', '#d01c8b')
colpal <- c('#0072b2', '#f1b6da', '#e41a1c')
#
# standardized fold change vector 
fc_heart <- GeneTest$z.heart.prot
names(fc_heart) <- GeneTest$gene.id
# truncate range at +/- 5
fc_heart <- ifelse(fc_heart < -5, -5, fc_heart)
fc_heart <- ifelse(fc_heart > 5, 5, fc_heart)
#
fc_kidney <- GeneTest$z.heart.prot
names(fc_kidney) <- GeneTest$gene.id
#
# truncate range at +/- 5
fc_kidney <- ifelse(fc_kidney < -5, -5, fc_kidney)
fc_kidney <- ifelse(fc_kidney > 5, 5, fc_kidney)

#####
# pick a gene list to study
gene.list <- list(prot.up.up, mrna.up.up,
                  prot.dn.dn, mrna.dn.dn,
                  prot.dn.up, mrna.dn.up,
                  prot.up.dn, mrna.up.dn)

list.name <- c("Protein: up heart, up kidney",
               "RNA: up heart, up kidney",
               "Protein: down heart, down kidney",
               "RNA: down heart, down kidney",
               "Protein: down heart, up kidney",
               "RNA: down heart, up kidney",
               "Protein: up heart, down kidney",
               "RNA: up heart, down kidney")

# choose a category
ont.cat = c("BP","CC","MF")

#####
# run enrichment for the GO categories
pdf("Downstream_Analysis/Results/UpDownEnrich.pdf")
for(i in 1:length(gene.list)){
  for(j in 1:length(ont.cat)){
    
    enrich = WrapEnrichGO(gene.list[[i]], all.geneids, ont=ont.cat[j])
    enrich = enrich.filter(enrich, 0.01)
    
    if(dim(enrich)[1] > 0){
      print(barplot(enrich.trim(enrich,42), showCategory = 20) +
              ggtitle(str_c(list.name[i], ont.cat[j], sep=" ")))
      
      print(cnetplot(enrich,colorEdge=FALSE,
                     showCategory = 8, node_label="all",
                     foldChange = fc_heart) + guides(size = "none") +
              scale_color_gradientn(colors=colpal, name="heart z score") +
              ggtitle(str_c(list.name[i], ont.cat[j], sep=" ")))  
    }
  }
}
dev.off()



#########################################
# 
# Aging Heart 
# Script to compare the RNA and protein data
#   Supplemental Figure 3
# 
# Gary A Churchill, 2021
#
#########################################

#########################################
# set up the environment
#setwd("~/Work/Projects/Aging_Center_2021/Aging_Heart_v3")

###
library(tidyverse)
library(cowplot)
library(ggbeeswarm)
library(modelr)
library(broom)




#########################################################
#  Create a common set of genes with both RNA and Protein data
#########################################################


####
# shortcut - load saved version of merged data
load(file="Downstream_Analysis/Data/Merged_RNA_Protein.Rdata")
ls()
# "annot.gene"    "annot.samples" "mrna.set"      "protein.set" 


# ####
# # commented code is not needed unless you are constructing combined data
# load("QTLviewer_dataset/JAC_DO_heart_v9.gz")
# ls()
#
# ####
# # rename the datasets
# dataset.mrna <- dataset.DOheart.mrna
# dataset.protein <- dataset.DOheart.protein
# rm(dataset.DOheart.mrna, dataset.DOheart.protein)
# 
# ####
# # copy annotation files up one level to reduce typing...
# annot.mrna <- dataset.mrna$annot.mrna
# annot.protein <- dataset.protein$annot.protein
# 
# ####
# # combine annotation on the common set
# 
# # select proteins that have mrna
# tmp1 <- filter(annot.protein, gene.id %in% annot.mrna$gene.id)
# dim(tmp1)
# length(unique(tmp1$gene.id))
# 
# # # lets have a look at the duplicated proteins
# # dup.genes <- tmp1$gene.id[duplicated(tmp1$gene.id)]
# # dup.prots <- filter(tmp1, gene.id %in% dup.genes)$protein.id
# # dup.data <- dataset.protein$data$norm[,dup.prots]
# # colnames(dup.data) <- filter(tmp1, gene.id %in% dup.genes)$symbol
# # #
# # quartz()
# # plot(hclust(dist(t(dup.data))))
# # # interesting that they don't pair up nicely...
# # rm(dup.genes, dup.prots, dup.data)
# 
# #select mrna that have proteins
# tmp2 <- filter(annot.mrna, gene.id %in% annot.protein$gene.id)
# dim(tmp2)
# 
# # create the merged gene annotations
# # include protein with multiple rna by duplicating the rna data
# annot.gene <- left_join(select(tmp1, protein.id, gene.id), tmp2)
# annot.gene
# # protein.id is unique
# 
# 
# ####
# # subset the data matrices
# 
# # subset the protein using protein.id colnames
# protein.set <- dataset.protein$data$norm[,annot.gene$protein.id]
# dim(protein.set)
# 
# #subset the rna using gene.id colnames
# mrna.set <- dataset.mrna$data$norm[,unique(annot.gene$gene.id)]
# dim(mrna.set)
# # # replace colnames with unique identifier
# # colnames(mrna.set) <- annot.gene$protein.id
# 
# # pull the sample annotations
# annot.samples <- dataset.mrna$annot.samples
# 
# ###
# #clean up
# rm(tmp1, tmp2)
# rm(annot.mrna, annot.protein)
# rm(ensembl.version)
# rm(dataset.mrna, dataset.protein)
# rm(genoprobs, K, map, markers)
# ls()
# 
# save(file="Data/Merged_RNA_Protein.Rdata", list=ls())

# end of Merge Data block
#########################################################



#########################################################
# compute within gene variance and correlations for RNA and Protein
#########################################################

#######
# create a nested data structure for gene-level modeling

# protein data
tmp1 <- protein.set %>% as_tibble() %>%
  mutate(mouse.id = row.names(protein.set)) %>%
  left_join(select(annot.samples, mouse.id)) %>%
  select(mouse.id, everything()) %>%
  pivot_longer(2:4118, names_to = "protein.id", values_to = "protein")

# rna data
tmp2 <- mrna.set %>% as_tibble() %>%
  mutate(mouse.id = row.names(protein.set)) %>%
  left_join(select(annot.samples, mouse.id, Sex, Age)) %>%
  select(mouse.id, Sex, Age, everything()) %>%
  mutate(Age = factor(Age)) %>%
  pivot_longer(4:4050, names_to = "gene.id", values_to = "mrna") %>%
  left_join(select(annot.gene, protein.id, gene.id, symbol, chr))

# join them into one structure and nest
by_protein <- left_join(tmp1, tmp2) %>%
  select(protein.id, gene.id, symbol, chr, mouse.id, Sex, Age, mrna, protein) %>%
  group_by(protein.id, gene.id, symbol, chr) %>%
  nest()
#
by_protein <- filter(by_protein, !(chr=="Y" | chr=="X"))

rm(tmp1, tmp2)


#####
# grab two gene-level datasets for testing 
# Cav2 has missing protein data
Apoh <- filter(by_protein, symbol=="Apoh")$data[[1]]
Cav2 <- filter(by_protein, symbol=="Cav2")$data[[1]]


####################
# does RNA-protein correlation change with age on a gene-by-gene basis?

# function returns p-val for interaction test
my.test1 <- function(df){
  anova(lm(protein ~ Sex + Age + mrna, data = df),
        lm(protein ~ Sex + Age + mrna + Age:mrna, data=df))[2,6]
}
#
my.test1(Cav2)


# apply the test to all genes
by_protein <- by_protein %>%
  # mutate(p.test1 = simplify(map(data, my.test1))) %>%
  mutate(p.test1 = map(data, my.test1)) %>%
  unnest(p.test1) 
#
p.adjust(by_protein$p.test1, method="BH") %>% sort() %>% head(20)
# there are 0 (FDR<0.01) or 3(FDR < 0.1) significant 
# changes in RNA to Protein correlation with age
#

# # pvalue histogram for "does rna-prot corr change with age?"
# quartz()
# ggplot(by_protein, aes(x=p.test1)) +
#   geom_histogram(bins=100)
# #
# # there is no evidence that protein-rna correlations change with age


##########
# compute age and sex adjusted correlation for each gene

####
# function to adjust for Sex and Age
my.adj.corr <- function(df){
  res.mrna <- residuals(lm(mrna ~ Age + Sex + Age:Sex, data = df, na.action=na.exclude))
  res.prot <- residuals(lm(protein ~ Age + Sex + Age:Sex, data = df, na.action=na.exclude))
  cor(res.mrna, res.prot, use="pair")
}
# 
my.adj.corr(Apoh)
my.adj.corr(Cav2)
#
by_protein <- by_protein %>%
  mutate(corr.adj = purrr::simplify(map(data, my.adj.corr)))
#
quartz()
ggplot(by_protein, aes(x=corr.adj)) +
  geom_histogram(bins=100)
# here we see that sex and age adjusted RNA-prot correlations are generally positive


##########
# compute variance & correlation by Age

my.vars <- function(df){
  
  df <-  group_by(df, Age, Sex)
  
  summarise(df,
            var.mrna = var(mrna, na.rm=TRUE),
            var.prot = var(protein, na.rm=TRUE),
            cor.mrna.prot = cor(mrna, protein, use="pair"),
            .groups = "keep")
}
#
my.vars(Apoh)
#
Vars <- by_protein %>% 
  mutate(vars = map(data, my.vars)) %>%
  unnest(vars, .drop=TRUE) %>%
  select(-ends_with("test1"), -data) %>%
  ungroup()

###
# look at the median of SD and r across genes
Vars %>% group_by(Age) %>% 
  summarize(m.sd.mrna = median(sqrt(var.mrna)),
            m.sd.prot = median(sqrt(var.prot)),
            m.corr = median(cor.mrna.prot),
            cor.var = cor(log(var.mrna), log(var.prot), use="pair"))


# break down by Age and Sex is interesting...
Vars %>% group_by(Age, Sex) %>%
  summarize(m.sd.mrna = median(sqrt(var.mrna)),
            m.sd.prot = median(sqrt(var.prot)),
            m.corr = median(cor.mrna.prot),
            cor.var = cor(log(var.mrna), log(var.prot), use="pair"))


#########
# Construct Suppl Figure 3

#
p1 <- ggplot(Vars, aes(x=Age, y=log2(var.mrna))) + 
  geom_beeswarm(cex=0.3)  +
  geom_boxplot(outlier.shape = NA, alpha = 0.1, aes(color=Age)) +
  # theme_bw(base_size = 18) + 
  guides(color="none") +
  # facet_wrap(~Sex) +
  ggtitle("mRNA log Variance by Age")
quartz(); p1
#
p2 <- ggplot(Vars, aes(x=Age, y=log2(var.prot))) + 
  geom_beeswarm(cex=0.3)  +
  geom_boxplot(outlier.shape = NA, alpha = 0.1, aes(color=Age)) +
  # theme_bw(base_size = 18) + 
  guides(color="none") +
  # facet_wrap(~Sex) +
  ggtitle("Protein log Variance by Age")
quartz(); p2
#
p3 <- ggplot(Vars, aes(x=Age, y=cor.mrna.prot)) + 
  geom_beeswarm(cex=0.35)  +
  geom_boxplot(outlier.shape = NA, alpha = 0.1, aes(color=Age)) +
  # theme_bw(base_size = 18) + 
  guides(color="none") +
  # facet_wrap(~Sex) +
  ggtitle("mRNA-Protein Correlation by Age")
quartz(); p3
# #
# p4 <- ggplot(Vars, aes(x=log2(var.mrna), y=log2(var.prot))) + 
#   geom_point()  +
#   geom_smooth(method="lm", se=FALSE, aes(color=Age)) + 
#   guides(color="none") +
#   # facet_wrap(~Sex) +
#   ggtitle("mRNA v Protein log Variance by Age")
# quartz(); p4

quartz()
plot_grid(p1, p2, p3, labels = c("A", "B", "C"), ncol=2)







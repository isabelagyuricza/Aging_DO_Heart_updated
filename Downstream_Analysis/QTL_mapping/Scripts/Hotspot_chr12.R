########### Chromosome 12 hotspot - QTL analysis on the DO heart ################

# Analyzing the correlation between proteins on the hotspot on chromosme 12 
# and doing enrichment analysis.

################################################################################

# Author: Isabela Gerdes Gyuricza - Churchill lab 
# Date: 11-10-2021

################################################################################
############ clear workspace

rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

################################################################################
############ libraries

library(RColorBrewer)
library(corrplot)
library(qtl2)
library(pcaMethods)
library(tidyverse)

source("~/Box/JAC_Heart_Data/GBRS/QTL_viewer/heart_newversion/Results/Transbands/Drivers_allele_effects/facet_CC_colors_function.R")


################################################################################
############ load data

load("QTLviewer_dataset/JAC_DO_heart_v9.gz")

################################################################################

### for transbands of protein interactive scan for age

### Q x Age interactions on Chr 12

#####################
###  extract the lod peaks and start trimming them to just the hotspot

lod.peaks <- dataset.protein$lod.peaks$age_int %>%
  data.frame() %>% 
  filter(substr(marker.id,1,3) == "12_") %>%
  separate(marker.id, into=c("qtl.chr","qtl.pos")) %>%
  mutate(qtl.chr=as.integer(qtl.chr), qtl.pos=as.double(qtl.pos)/(10^6)) %>%
  left_join(dataset.protein$annot.protein)

# use this graphic to fine tune the edges of the hotspot
#quartz()
ggplot(lod.peaks, aes(x=qtl.pos, y=lod)) +
  geom_point() +
  xlim(98,102) #Region with the highest lod scores ~145Mb to 149Mb
#
lod.peaks <- dplyr::filter(lod.peaks, qtl.pos>98 & qtl.pos<102)
dim(lod.peaks)
# 224 genes

############################
# grab normalized mRNA expr data, regress out Sex (adjusting for sex), then rankz it.

# expression data
expr.data <- dataset.protein$data$norm[,lod.peaks$protein.id]

# compute residuals to adjust x wrt covariate
myresids <- function(x, cov){
  residuals(lm(x~cov,na.action=na.exclude))
}

#
expr.data <- apply(expr.data,2,myresids,dataset.protein$annot.samples$Sex)

# Apply a normal scores transform to the data.

rankZ <- function(x) {
  x <- rank(x, na.last = "keep") / (length(x) - sum(is.na(x)) + 1)
  return(qnorm(x))
}

expr.data <- apply(expr.data,2, rankZ) #Across mice

rm(myresids, rankZ)

###
# order the data by QTL peak location
lod.peaks <- arrange(lod.peaks, by = qtl.pos)
expr.data <- expr.data[,lod.peaks$protein.id]
###

############################
# look at correlations

# grab average absolute correlations and filter the lod peaks
tmp.m <- apply(abs(cor(expr.data, use = "pairwise.complete.obs")),2,mean)
sum(tmp.m > 0.3)  #194
sum(tmp.m > 0.35) #177
sum(tmp.m > 0.4)  #156
#
lod.peaks <- mutate(lod.peaks, mean.corr = tmp.m)
#
lod.peaks <- filter(lod.peaks, mean.corr>0.3)

expr.data <- expr.data[,lod.peaks$protein.id]

################################################################################

# Saving data to compare with the proteins from the hotspot on chr 3 

saveRDS(lod.peaks, 
        file = "Downstream_Analysis/QTL_mapping/Results/proteins_hotspot_chr12.RDS")

################################################################################

# hierarchical clustering
quartz()
hm <- gplots::heatmap.2(cor(expr.data, use = "pairwise.complete.obs"),
                        Rowv = TRUE,
                        Colv=TRUE,dendrogram ="column",
                        trace="none",
                        sepcolor = "black",
                        keysize = 0.9, 
                        key.title = "Color_key",
                        margins = c(8,10),
                        col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))

#There are a group of genes that seem to be positively correlated.

# reorder genes in lod.peaks and expr.data
lod.peaks <- lod.peaks[hm$colInd,]
expr.data <- expr.data[,hm$colInd]

############# Generating supplemental figure for the correlations ##############

test <- expr.data %>% 
  data.frame() %>% 
  rownames_to_column("mouse.id") %>% 
  gather("protein.id","expr",-mouse.id) %>% 
  left_join(lod.peaks) %>% 
  select(mouse.id, protein.id, symbol, expr)

# Transforming duplicated symbols into unique

symbols <- test %>% 
  select(protein.id, symbol) %>% 
  unique() %>% 
  mutate(symbol = make.unique(symbol)) %>% 
  select(protein.id, symbol)

test <- test %>% 
  select(-symbol) %>% 
  left_join(symbols) %>% 
  select(-protein.id) %>% 
  spread(symbol, expr) %>% 
  column_to_rownames("mouse.id") %>% 
  as.matrix()

pdf("Downstream_Analysis/QTL_mapping/Results/Correlation_hotspot_chr12.pdf",
    width = 10,
    height = 10)

#quartz()
gplots::heatmap.2(cor(test, use = "pairwise.complete.obs"),
                        Rowv = TRUE,
                        Colv=TRUE,
                        dendrogram ="column",
                        trace="none",
                        keysize = 1, 
                        key.title = "Color_key",
                        col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(128))

dev.off()

rm(symbols, test)

################################################################################

# compute a PCA summary of the hotspot genes

# Each mouse will have one value for the hotspot (protein expression)
pca.chr12 <- pca(expr.data[,lod.peaks$protein.id], scale="uv", center=TRUE, nPcs=8)
summary(pca.chr12)
# first PC explains 50%

# save PC1 loadings
tmp <- loadings(pca.chr12)
tmp2 <- data_frame(protein.id=row.names(tmp),PC1 = tmp[,"PC1"])
lod.peaks <- left_join(lod.peaks, tmp2)
rm(tmp, tmp2)

# create a datframe for phenotypes and add PCs
pheno <- dataset.protein$annot.samples %>%
  dplyr::select(mouse.id, Sex, Age) %>%
  mutate(PC1 = scores(pca.chr12)[,1])

###
# map PC1, to confirm the signal on the hotspot!

scan.additive <- scan1(genoprobs, pheno$PC1, kinship = K,
                       addcovar = dataset.protein$covar.matrix)
scan.interaction <- scan1(genoprobs, pheno$PC1, kinship = K,
                          addcovar = dataset.protein$covar.matrix,
                          intcovar = dataset.protein$covar.matrix[,"Age"])
scan.out <- cbind(scan.additive, scan.interaction, scan.interaction-scan.additive)
colnames(scan.out) <- c("add","full","diff")
rm(scan.additive, scan.interaction)


which.max(scan.out[,3])
#12_99450809 
#44939  Using this marker for the model.

#Finding the intervals

bayes_int(scan1_output = scan.out, lodcolumn = 3, prob = 0.99,map = map, chr = 12)
# ci_lo      pos    ci_hi
# 1 98.74842 99.45081 100.3528

# Using fit_1 to get the allele effects

# 1) Using age 6 as baseline and getting the combined coefficients dataframe

annot.samples <- dataset.protein$annot.samples %>% 
  mutate(fAge = as.factor(Age))

ac <- model.matrix(~fAge+Sex, 
                   data=annot.samples)[,-1]

row.names(ac) <- annot.samples$mouse.id

ic <- model.matrix(~fAge, data=annot.samples)[,-1]

row.names(ic) <- annot.samples$mouse.id

fit <- fit1(genoprobs$`12`[,,"12_99450809"], pheno$PC1, kinship = K[12],
            addcovar = ac, intcovar = ic,zerosum = FALSE)

# 1) Using age 6 as baseline

fit_6 <- data.frame(estimate = fit$coef[LETTERS[1:8]], 
                    se = fit$SE[LETTERS[1:8]], Age = 6) %>% 
  rownames_to_column("Founder")


# 2) Using age 12 as baseline 

annot.samples$fAge <- factor(annot.samples$fAge, levels = c(12,6,18))

ac <- model.matrix(~fAge+Sex, 
                   data=annot.samples)[,-1]

row.names(ac) <- annot.samples$mouse.id

ic <- model.matrix(~fAge, data=annot.samples)[,-1]

row.names(ic) <- annot.samples$mouse.id

fit_2 <- fit1(genoprobs$`12`[,,"12_99450809"], pheno$PC1, kinship = K[12],
              addcovar = ac, intcovar = ic, zerosum = FALSE)

fit_12 <- data.frame(estimate = fit_2$coef[LETTERS[1:8]], 
                     se = fit_2$SE[LETTERS[1:8]], Age = 12) %>% 
  rownames_to_column("Founder")


# 3) Using age 18 as baseline 

annot.samples$fAge <- factor(annot.samples$fAge, levels = c(18,6,12))

ac <- model.matrix(~fAge+Sex, 
                   data=annot.samples)[,-1]

row.names(ac) <- annot.samples$mouse.id

ic <- model.matrix(~fAge, data=annot.samples)[,-1]

row.names(ic) <- annot.samples$mouse.id

fit_3 <- fit1(genoprobs$`12`[,,"12_99450809"], pheno$PC1, kinship = K[12],
              addcovar = ac, intcovar = ic, zerosum = FALSE)

fit_18 <- data.frame(estimate = fit_3$coef[LETTERS[1:8]], 
                     se = fit_3$SE[LETTERS[1:8]], Age = 18) %>% 
  rownames_to_column("Founder")


# Gathering them all in the same dataframe.

final_df <- bind_rows(fit_6,fit_12,fit_18)

final_df$Age <- factor(final_df$Age, levels = c(6, 12, 18))

###
# plot the regression coefs

pdf("Downstream_Analysis/QTL_mapping/Results/PC1_hotspot_chr12.pdf")

data(CCcolors)
final_df %>% 
  mutate(Founder = factor(Founder, levels = c(LETTERS[1:8]), 
                          labels = c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/LtJ", 
                                     "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ"))) %>% 
  ggplot(aes(x=Age, y=estimate)) +
  geom_point() +
  geom_line(aes(group=Founder), linetype = "dashed") +
  facet_wrap(~Founder, ncol=4) +
  geom_errorbar(aes(x=Age, ymin=estimate-se, ymax=estimate+se), width=0.2) +
  theme_bw() + ggtitle("PC1") -> p

p <- facet_CC_colour(p, ncol = 4)

#quartz()
grid.draw(p)

dev.off()

rm(final_df,fit,fit_12,fit_18,fit_2,fit_3,fit_6,ac,ic)

rm(annot.samples)

########################## Enrichment analysis #################################

library(org.Mm.eg.db)
library(clusterProfiler)

# library loads this will mask some dplyr functions

enrichment <- enrichGO(
  gene = unique(lod.peaks$gene.id),
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "all",
  universe = dataset.protein$annot.protein$gene.id,
  qvalueCutoff = 0.1,
  pvalueCutoff = 0.1,
  pAdjustMethod = "BH",
  readable = TRUE)

enrichment_df <- enrichment@result[enrichment@result$p.adjust < 0.05,]

#################### saving it for supplemental material #######################

save <- enrichment_df %>% 
  rename(Ontology_therm = ONTOLOGY, 
         Category_ID = ID, 
         Category_description = Description,
         Gene_ratio = GeneRatio,
         Bg_ratio = BgRatio, 
         p = pvalue, 
         adjusted_p = p.adjust,
         Gene_IDs = geneID) %>% 
  select(-qvalue) 

write.csv(save, "Downstream_Analysis/QTL_mapping/Results/Enrichment_hotspot_chr12.csv",
          row.names = FALSE)

# There are too many duplicated categories (with the same genes)

list_symbols <- list()

for(i in unique(enrichment_df$ID)) {
  
  list <- enrichment_df %>% filter(ID == i) %>%
    select(geneID) %>%
    unlist()
  
  list <- str_split(list, pattern = "/", simplify = TRUE) %>% as.vector()
  
  list_symbols[[i]] <- list
  
}

enrichment_df <- plyr::ldply(list_symbols, rbind) %>%
  gather("rm","Gene_symbol", -`.id`) %>%
  rename(ID = .id) %>%
  select(-rm) %>%
  filter(!is.na(Gene_symbol)) %>% 
  inner_join(enrichment_df) %>% 
  select(-geneID) %>% 
  arrange(p.adjust)

enrichment_df <- enrichment_df[!duplicated(enrichment_df$Gene_symbol),]

# Plotting these pathways.. 

categories <- enrichment[enrichment$ID %in% enrichment_df$ID] %>% 
  arrange(p.adjust)

tmp <- enrichment

tmp@result <- categories

pdf("Downstream_Analysis/QTL_mapping/Results/Enrichment_hotspot_chr12.pdf",
    width = 8, height = 6)

dotplot(tmp, color = "p.adjust", font.size = 1.2) +
  scale_color_gradientn(colours = rev(c('#b3cde3','#8c96c6','#8856a7','#810f7c'))) +
  scale_x_continuous(limits = c(0,0.15), breaks = c(0, 0.05, 0.1, 0.15), expand = c(0,0)) 

dev.off()

################################################################################

# Comparing proteins between hotspot chromosome 3 and 12. 

rm(list = ls())

load("QTLviewer_dataset/JAC_DO_heart_v9.gz")


chr3 <- readRDS("Downstream_Analysis/QTL_mapping/Results/proteins_hotspot_chr3.RDS")

chr12 <- readRDS("Downstream_Analysis/QTL_mapping/Results/proteins_hotspot_chr12.RDS")

common <- intersect(chr3$protein.id, chr12$protein.id)

dataset.protein$annot.protein %>% filter(protein.id %in% common) %>% View()


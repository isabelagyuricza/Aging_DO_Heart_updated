################ DEA and enrichment analysis on the DO heart ###################

# Running differential expression analysis (DEA) with age followed by enrichment
# analysis on the heart.

################################################################################

# Author: Isabela Gerdes Gyuricza - Churchill lab 
# Date: 09-16-2021

################################################################################
############ clear workspace

rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

################################################################################
############ loading libraries

library(DESeq2)
library(tidyverse)

# Load QTLviewer data 

load("QTLviewer_dataset/JAC_DO_heart_v9.gz")

# Deleting what I'm not going to use. 

rm(genoprobs, K, map, markers, ensembl.version)

# 1) First, running DEA for the transcripts

annot.samples <- dataset.mrna$annot.samples %>% 
  mutate(fAge = factor(Age, levels=c("6","12","18")), # Create factor Age to use
         # in different tests.
         Age = (Age/12)-1) # Rescale age to units per year


# Make sure annot.samples and expr matrix are in the same order 
identical(annot.samples$mouse.id, rownames(dataset.mrna$data$raw))
#[1] TRUE

# 1) Testing for Age

############################ FULL SET ##########################################

dds <- DESeqDataSetFromMatrix(countData = t(
                                floor(dataset.mrna$data$raw)
                                      ),
                              colData = annot.samples,
                              design = ~ Sex + Age)

# Controlling for sex while testing for age (as continous) by LRT.

dds_de <- DESeq(dds, test="LRT", reduced= ~ Sex)

# Generating data and saving it for supplemental table

result_age <- results(dds_de) %>% 
  data.frame() %>% 
  rownames_to_column("gene.id") %>% 
  left_join(dataset.mrna$annot.mrna %>% 
              dplyr::select(gene.id, symbol)) %>% 
  dplyr::select(gene.id, symbol, log2FoldChange, lfcSE, pvalue) %>% 
  # Recomputing the p-values. The adjustment of DESEq2 is a bit weird
  rename(Gene_ID = gene.id, 
         Gene_symbol = symbol,
         Effect = log2FoldChange,
         Effect_SE = lfcSE,
         p = pvalue) %>% 
  mutate(adjusted_p = p.adjust(p, method = "BH"),
         Gene_symbol = paste0("\"", Gene_symbol, "\""))

write.csv(result_age,
          file = "Downstream_Analysis/Results/DEA_age_FULL_transcripts.csv",
          row.names = FALSE)

################################################################################

# Running DEA only for transcript that have proteins.. 

common_set <- intersect(dataset.mrna$annot.mrna$gene.id, 
                        dataset.protein$annot.protein$gene.id)

data <- t(dataset.mrna$data$raw) %>% 
  data.frame() %>% 
  rownames_to_column("gene.id") %>% 
  filter(gene.id %in% common_set) %>% 
  column_to_rownames("gene.id") %>% 
  as.matrix() %>% 
  floor()
# We have 4047 transcripts that have proteins

identical(annot.samples$mouse.id, colnames(data))
#[1] TRUE

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = annot.samples,
                              design = ~ Sex + Age) 


# Controlling for sex while testing for age (as continous) by LRT.

dds_de <- DESeq(dds, test="LRT", reduced= ~ Sex)

# Generating data and saving it for supplemental table

result_age <- results(dds_de) %>% 
  data.frame() %>% 
  rownames_to_column("gene.id") %>% 
  left_join(dataset.mrna$annot.mrna %>% 
              dplyr::select(gene.id, symbol)) %>% 
  dplyr::select(gene.id, symbol, log2FoldChange, lfcSE, pvalue) %>% 
  # Recomputing the p-values. The adjustment of DESEq2 is a bit weird
  rename(Gene_ID = gene.id, 
         Gene_symbol = symbol,
         Effect = log2FoldChange,
         Effect_SE = lfcSE,
         p = pvalue) %>% 
  mutate(adjusted_p = p.adjust(p, method = "BH"),
         Gene_symbol = paste0("\"", Gene_symbol, "\""))

write.csv(result_age, file = "Downstream_Analysis/Results/DEA_age_transcripts.csv",
          row.names = FALSE)

result_age %>% filter(adjusted_p < 0.01) %>% dim()
# [1] 206   6


# 2) Testing for sex 

############################ FULL SET ##########################################

dds <- DESeqDataSetFromMatrix(countData = t(
                                     floor(dataset.mrna$data$raw)
                                      ),
                              colData = annot.samples,
                              design = ~ fAge + Sex)

dds_de <- DESeq(dds, test="LRT", reduced= ~ fAge)

# Generating data and saving it for supplemental table

result_sex <- results(dds_de) %>% 
  data.frame() %>% 
  rownames_to_column("gene.id") %>% 
  left_join(dataset.mrna$annot.mrna %>% 
              dplyr::select(gene.id, symbol)) %>% 
  dplyr::select(gene.id, symbol, log2FoldChange, lfcSE, pvalue) %>% 
  # Recomputing the p-values. The adjustment of DESEq2 is a bit weird
  rename(Gene_ID = gene.id, 
         Gene_symbol = symbol,
         Effect = log2FoldChange,
         Effect_SE = lfcSE,
         p = pvalue) %>% 
  mutate(adjusted_p = p.adjust(p, method = "BH"),
         Gene_symbol = paste0("\"", Gene_symbol, "\""))

write.csv(result_sex,
          file = "Downstream_Analysis/Results/DEA_sex_FULL_transcripts.csv",
           row.names = FALSE)

################################################################################

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = annot.samples,
                              design = ~ fAge + Sex) # Using factor Age here

dds_de <- DESeq(dds, test="LRT", reduced= ~ fAge)

# Generating data and saving it for supplemental table

result_sex <- results(dds_de) %>% 
  data.frame() %>% 
  rownames_to_column("gene.id") %>% 
  left_join(dataset.mrna$annot.mrna %>% 
              dplyr::select(gene.id, symbol)) %>% 
  dplyr::select(gene.id, symbol, log2FoldChange, lfcSE, pvalue) %>% 
  # Recomputing the p-values. The adjustment of DESEq2 is a bit weird
  rename(Gene_ID = gene.id, 
         Gene_symbol = symbol,
         Effect = log2FoldChange,
         Effect_SE = lfcSE,
         p = pvalue) %>% 
  mutate(adjusted_p = p.adjust(p, method = "BH"),
         Gene_symbol = paste0("\"", Gene_symbol, "\""))

write.csv(result_sex, file = "Downstream_Analysis/Results/DEA_sex_transcripts.csv",
           row.names = FALSE)

result_sex %>% filter(adjusted_p < 0.01) %>% dim()
#[1] 887   6

# 3) Testing for sex:age interaction

############################ FULL SET ##########################################

dds <- DESeqDataSetFromMatrix(countData = t(
                                         floor(dataset.mrna$data$raw)
                                          ),
                              colData = annot.samples,
                              design = ~ Age + Sex + Age:Sex)

dds_de <- DESeq(dds, test="LRT", reduced= ~ Age + Sex)

# Generating data and saving it for supplemental table

result_age_sex <- results(dds_de) %>% 
  data.frame() %>% 
  rownames_to_column("gene.id") %>% 
  left_join(dataset.mrna$annot.mrna %>% 
              dplyr::select(gene.id, symbol)) %>% 
  dplyr::select(gene.id, symbol, log2FoldChange, lfcSE, pvalue) %>% 
  # Recomputing the p-values. The adjustment of DESEq2 is a bit weird
  rename(Gene_ID = gene.id, 
         Gene_symbol = symbol,
         Effect = log2FoldChange,
         Effect_SE = lfcSE,
         p = pvalue) %>% 
  mutate(adjusted_p = p.adjust(p, method = "BH"),
         Gene_symbol = paste0("\"", Gene_symbol, "\""))

write.csv(result_age_sex,
          file = "Downstream_Analysis/Results/DEA_int_FULL_transcripts.csv",
           row.names = FALSE)

################################################################################

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = annot.samples,
                              design = ~ Age + Sex + Age:Sex)

# Controlling for sex while testing for age (as continous) by LRT.

dds_de <- DESeq(dds, test="LRT", reduced= ~ Age + Sex)

# Generating data and saving it for supplemental table

result_age_sex <- results(dds_de) %>% 
  data.frame() %>% 
  rownames_to_column("gene.id") %>% 
  left_join(dataset.mrna$annot.mrna %>% 
              dplyr::select(gene.id, symbol)) %>% 
  dplyr::select(gene.id, symbol, log2FoldChange, lfcSE, pvalue) %>% 
  # Recomputing the p-values. The adjustment of DESEq2 is a bit weird
  rename(Gene_ID = gene.id, 
         Gene_symbol = symbol,
         Effect = log2FoldChange,
         Effect_SE = lfcSE,
         p = pvalue) %>% 
  mutate(adjusted_p = p.adjust(p, method = "BH"),
         Gene_symbol = paste0("\"", Gene_symbol, "\""))

write.csv(result_age_sex, file = "Downstream_Analysis/Results/DEA_int_transcripts.csv",
           row.names = FALSE)

result_age_sex %>% filter(adjusted_p < 0.01) %>% dim()
#[1] 3 6

rm(dds, dds_de, data)

# Gathering all the results in a list 

results <- list(age_rna = result_age,
                sex_rna = result_sex, 
                int_rna = result_age_sex)

rm(result_age, result_age_sex, result_sex)

######################## Running DEA for proteins now 

# Setting up anova function

anova_comput <- function (data) {
  
  df_comput <- data$data$norm %>% 
    data.frame() %>% 
    rownames_to_column("mouse.id") %>% 
    left_join(dataset.protein$annot.samples) %>% 
    dplyr::select(mouse.id, Sex, Age,contains("ENSMUS")) %>% 
    mutate(Age = (Age/12) - 1)
  
  df_value <- df_comput %>% dplyr::select(contains("ENSMUS"))
  
  add <- function(x){
    
    fit_additive <- summary(lm(x ~ Sex + Age, data=df_comput))$coefficients
    
    df_final_additive <- tibble(pvalue.Age = fit_additive["Age","Pr(>|t|)"],
                                estimate.Age = fit_additive["Age","Estimate"],
                                se.Age = fit_additive["Age","Std. Error"],
                                pvalue.Sex = fit_additive["SexM","Pr(>|t|)"],
                                estimate.Sex = fit_additive["SexM","Estimate"],
                                se.Sex = fit_additive["SexM","Std. Error"])
    
    return(df_final_additive)}
  
  result_add <- apply(df_value, 2, add)
  
  int <- function(x){
    
    fit_interactive <- summary(
      lm(
        x ~ Sex + Age + Sex:Age, data = df_comput))$coefficients
    
    df_final_interactive <- tibble(
      pvalue.Int = fit_interactive["SexM:Age","Pr(>|t|)"],
      estimate.Int = fit_interactive["SexM:Age","Estimate"],
      se.Int = fit_interactive["SexM:Age","Std. Error"])
    
    return(df_final_interactive)}
  
  result_int <- apply(df_value, 2, int)
  
  result <- list(result_add,result_int)
  
  names(result) <- c("additive","interactive")
  
  return(result)
}

############################ FULL SET ##########################################

result <- anova_comput(dataset.protein)

# Generating data for supplemental table

result_age <- result$additive %>% 
  bind_rows(.id = "protein.id") %>% 
  select(protein.id, estimate.Age, se.Age, pvalue.Age) %>% 
  left_join(dataset.protein$annot.protein %>% 
              dplyr::select(protein.id, gene.id, symbol)) %>% 
  mutate(adjusted_p = p.adjust(pvalue.Age, method = "BH"),
         symbol = paste0("\"", symbol, "\"")) %>% 
  rename(Protein_ID = protein.id, 
         Gene_ID = gene.id,
         Gene_symbol = symbol,
         Effect = estimate.Age,
         Effect_SE = se.Age,
         p = pvalue.Age) %>% 
  select(Protein_ID, Gene_ID, Gene_symbol, Effect, Effect_SE, p, adjusted_p)

write.csv(result_age,
         file = "Downstream_Analysis/Results/DEA_age_FULL_proteins.csv",
          row.names = FALSE)

result_sex <- result$additive %>% 
  bind_rows(.id = "protein.id") %>% 
  select(protein.id, estimate.Sex, se.Sex, pvalue.Sex) %>% 
  left_join(dataset.protein$annot.protein %>% 
              dplyr::select(protein.id, gene.id, symbol)) %>% 
  mutate(adjusted_p = p.adjust(pvalue.Sex, method = "BH"),
         symbol = paste0("\"", symbol, "\"")) %>% 
  rename(Protein_ID = protein.id,
         Gene_ID = gene.id, 
         Gene_symbol = symbol,
         Effect = estimate.Sex,
         Effect_SE = se.Sex,
         p = pvalue.Sex) %>% 
  select(Protein_ID, Gene_ID, Gene_symbol, Effect, Effect_SE, p, adjusted_p)

write.csv(result_sex,
          file = "Downstream_Analysis/Results/DEA_sex_FULL_proteins.csv",
          row.names = FALSE)

result_age_sex <- result$interactive %>% 
  bind_rows(.id = "protein.id") %>% 
  left_join(dataset.protein$annot.protein %>% 
              dplyr::select(protein.id, gene.id, symbol)) %>% 
  mutate(adjusted_p = p.adjust(pvalue.Int, method = "BH"),
         symbol = paste0("\"", symbol, "\"")) %>% 
  rename(Protein_ID = protein.id, 
         Gene_ID = gene.id,
         Gene_symbol = symbol,
         Effect = estimate.Int,
         Effect_SE = se.Int,
         p = pvalue.Int) %>% 
  select(Protein_ID, Gene_ID, Gene_symbol, Effect, Effect_SE, p, adjusted_p)

write.csv(result_age_sex,
          file = "Downstream_Analysis/Results/DEA_int_FULL_proteins.csv",
          row.names = FALSE)

################################################################################

# Running DEA only for proteins that have transcripts

data <- t(dataset.protein$data$norm) %>% 
  data.frame() %>% 
  rownames_to_column("protein.id") %>% 
  left_join(dataset.protein$annot.protein %>% 
              select(gene.id, protein.id)) %>% 
  filter(gene.id %in% common_set) %>% 
  select(protein.id, contains("DO.")) %>% 
  column_to_rownames("protein.id") %>% 
  t() %>% 
  as.matrix()
# We have 4117 proteins that have transcripts

dataset.protein$data$norm <- data

result <- anova_comput(dataset.protein)

# Generating data for supplemental table

result_age <- result$additive %>% 
  bind_rows(.id = "protein.id") %>% 
  select(protein.id, estimate.Age, se.Age, pvalue.Age) %>% 
  left_join(dataset.protein$annot.protein %>% 
              dplyr::select(protein.id, gene.id, symbol)) %>% 
  mutate(adjusted_p = p.adjust(pvalue.Age, method = "BH"),
         symbol = paste0("\"", symbol, "\"")) %>% 
  rename(Protein_ID = protein.id, 
         Gene_ID = gene.id,
         Gene_symbol = symbol,
         Effect = estimate.Age,
         Effect_SE = se.Age,
         p = pvalue.Age) %>% 
  select(Protein_ID, Gene_ID, Gene_symbol, Effect, Effect_SE, p, adjusted_p)

write.csv(result_age, file = "Downstream_Analysis/Results/DEA_age_proteins.csv",
          row.names = FALSE)

result_age %>% filter(adjusted_p < 0.01) %>% dim()
#[1] 2084    7

result_sex <- result$additive %>% 
  bind_rows(.id = "protein.id") %>% 
  select(protein.id, estimate.Sex, se.Sex, pvalue.Sex) %>% 
  left_join(dataset.protein$annot.protein %>% 
              dplyr::select(protein.id, gene.id, symbol)) %>% 
  mutate(adjusted_p = p.adjust(pvalue.Sex, method = "BH"),
         symbol = paste0("\"", symbol, "\"")) %>% 
  rename(Protein_ID = protein.id,
         Gene_ID = gene.id, 
         Gene_symbol = symbol,
         Effect = estimate.Sex,
         Effect_SE = se.Sex,
         p = pvalue.Sex) %>% 
  select(Protein_ID, Gene_ID, Gene_symbol, Effect, Effect_SE, p, adjusted_p)

write.csv(result_sex, file = "Downstream_Analysis/Results/DEA_sex_proteins.csv",
          row.names = FALSE)

result_sex %>% filter(adjusted_p < 0.01) %>% dim()
#[1] 408   7

result_age_sex <- result$interactive %>% 
  bind_rows(.id = "protein.id") %>% 
  left_join(dataset.protein$annot.protein %>% 
              dplyr::select(protein.id, gene.id, symbol)) %>% 
  mutate(adjusted_p = p.adjust(pvalue.Int, method = "BH"),
         symbol = paste0("\"", symbol, "\"")) %>% 
  rename(Protein_ID = protein.id, 
         Gene_ID = gene.id,
         Gene_symbol = symbol,
         Effect = estimate.Int,
         Effect_SE = se.Int,
         p = pvalue.Int) %>% 
  select(Protein_ID, Gene_ID, Gene_symbol, Effect, Effect_SE, p, adjusted_p)

write.csv(result_age_sex, file = "Downstream_Analysis/Results/DEA_int_proteins.csv",
          row.names = FALSE)

result_age_sex %>% filter(adjusted_p < 0.01) %>% dim()
#[1] 144  7

results$age_protein <- result_age

results$sex_protein <- result_sex

results$int_protein <- result_age_sex


rm(result_age, result_age_sex, result_sex, result, anova_comput, data)


# Running enrichment analysis with fGSEA

# Loading new libraries

library(fgsea)
library(msigdbr)

set.seed(2021)

# # Get the mouse genesets

# all_gene_sets <- msigdbr(species = "Mus musculus")

all_gene_sets <- readRDS(file="Downstream_Analysis/Data/all_gene_sets.RDS")

# convert a gene_set to a fgsea Pathway format (named list)

tmp <- all_gene_sets %>% 
  filter(gs_subcat %in% c("GO:BP","GO:CC","GO:MF")) %>% 
  mutate(gs_name = factor(gs_name)) %>%
  dplyr::select(gs_name, ensembl_gene)

GO_Pathways <- split(tmp$ensembl_gene, tmp$gs_name)

rm(tmp)

enrich_function <- function(x) {
  
  # To run this analysis we standardized all the transcripts and use them all as input
  
  myRanks <- x$Effect/x$Effect_SE
  
  names(myRanks) <- x$Gene_ID
  
  myRanks <- jitter(myRanks, amount = 0.0001)
  
  myRanks <- sort(myRanks, decreasing = TRUE)
  
  fgseaRes <- fgsea(pathways = GO_Pathways, 
                    stats    = myRanks,
                    minSize  = 15,
                    maxSize  = 500)
  
  mainPathways <- fgseaRes[pathway %in% 
                             collapsePathways(fgseaRes[order(pval)][padj < 0.05], 
                                              GO_Pathways, myRanks)$mainPathways]
  list_result <- list(enrichment = fgseaRes, 
                      collapsed_enrichment = mainPathways,
                      rank = myRanks)
  
  return(list_result)
  
}

enrich_list <- lapply(results, enrich_function)

#save(enrich_list, file = "Downstream_Analysis/Results/enrichments_FULL.RData")

#save(enrich_list, file = "Downstream_Analysis/Results/enrichments.RData")

# For each test, add the genes symbols within the pathways to see what they are
# and plot some nice pathways.

# 1) Age

# ------------------------------ Transcripts

enrichment_age <- enrich_list$age_rna$collapsed_enrichment$leadingEdge

names(enrichment_age) <-  enrich_list$age_rna$collapsed_enrichment$pathway

temp <- plyr::ldply(enrichment_age, rbind) %>% 
  rename(pathway = .id) %>% 
  gather("test","ensembl_gene", -pathway) %>% 
  dplyr::select(-test) %>% 
  filter(!is.na(ensembl_gene)) %>% 
  left_join(all_gene_sets %>% select(ensembl_gene, gene_symbol)) %>% 
  unique() 

list_symbols <- list()

for(i in unique(temp$pathway)) {
  
  list <- temp %>% filter(pathway == i) %>% 
    select(gene_symbol) %>% 
    unlist()
  
  list <- knitr::combine_words(
    as.vector(list),
    sep = "/",
    and = "")
  
  list_symbols[[i]] <- list
  
}

rm(list, i, temp)

list_symbols <- bind_rows(list_symbols, .id = "pathway") %>% 
  gather("pathway",symbols_list)

enrichment_age <- enrich_list$age_rna$collapsed_enrichment %>% 
  left_join(list_symbols)

rm(list_symbols)

# Looking for the most significant and different pathways. 

enrichment_age %>% 
  arrange(padj) %>% 
  View()

# Saving enrichment results for supplemental data 

enrichment_age <- enrichment_age %>% 
  data.frame() %>% 
  dplyr::rename(Pathway_name = pathway,
                p = pval,
                adjusted_p = padj,
                Log2Err = log2err,
                Size = size,
                Gene_IDs = leadingEdge,
                Gene_symbols = symbols_list) %>% 
  arrange(adjusted_p) %>% 
  select(Pathway_name, Log2Err, ES, NES, Size, p, adjusted_p, Gene_IDs, Gene_symbols)


# data.table::fwrite(
#   enrichment_age,
#   file = "Downstream_Analysis/Results/enrichment_age_FULL_transcripts.csv",
#   sep = ",",
#   na = "NA")

data.table::fwrite(
  enrichment_age,
  file = "Downstream_Analysis/Results/enrichment_age_transcripts.csv",
  sep = ",",
  na = "NA")

# Now, merging enrichment with DEA results.

result_age <- read_csv("Downstream_Analysis/Results/DEA_age_transcripts.csv")

enrichment <- enrichment_age$Gene_IDs

names(enrichment) <-  enrichment_age$Pathway_name

temp <- plyr::ldply(enrichment, rbind) %>%
  rename(pathway = .id) %>%
  gather("test","ensembl_gene", -pathway) %>%
  dplyr::select(-test) %>%
  filter(!is.na(ensembl_gene)) %>%
  left_join(all_gene_sets %>% select(ensembl_gene, gene_symbol)) %>%
  unique()

list_symbols <- list()

for(i in unique(temp$pathway)) {
  
  list <- temp %>% filter(pathway == i) %>%
    select(gene_symbol) %>%
    unlist()
  
  names(list) <- NULL
  
  list_symbols[[i]] <- list
  
}

rm(list, i, temp)

list_symbols <- plyr::ldply(list_symbols, rbind) %>%
  gather("rm","Gene_symbol", -`.id`) %>%
  rename(pathway = .id) %>%
  select(-rm) %>%
  filter(!is.na(Gene_symbol))

merge <- list_symbols %>%
  inner_join(result_age %>%
               filter(adjusted_p < 0.01))

# What are the pathways with the largest number of significant genes?

merge %>% group_by(pathway) %>% count() %>% arrange(desc(n))
# # A tibble: 77 x 2
# # Groups:   pathway [77]
# pathway                                               n
# <chr>                                             <int>
#   1 GOBP_IMMUNE_EFFECTOR_PROCESS                         30
# 2 GOBP_CELL_ACTIVATION                                 26
# 3 GOBP_LEUKOCYTE_MEDIATED_IMMUNITY                     25
# 4 GOBP_DEFENSE_RESPONSE                                23
# 5 GOBP_PROTEIN_FOLDING                                 23
# 6 GOBP_EXOCYTOSIS                                      22
# 7 GOBP_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS    19
# 8 GOCC_SECRETORY_GRANULE                               19
# 9 GOBP_ENDOCYTOSIS                                     17
# 10 GOBP_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN     17
# # … with 67 more rows

# What are pathways with most negative mean effect? 

merge %>% group_by(pathway) %>% summarise(mean = mean(Effect)) %>% arrange(mean)
# # A tibble: 77 x 2
# pathway                                            mean
# <chr>                                             <dbl>
#   1 GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_HEAT     -0.446
# 2 GOBP_CHAPERONE_MEDIATED_PROTEIN_FOLDING          -0.311
# 3 GOMF_HEAT_SHOCK_PROTEIN_BINDING                  -0.276
# 4 GOBP_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN -0.243
# 5 GOMF_UNFOLDED_PROTEIN_BINDING                    -0.230
# 6 GOBP_MICROTUBULE_BUNDLE_FORMATION                -0.209
# 7 GOBP_PROTEIN_FOLDING                             -0.209
# 8 GOMF_PROTEIN_KINASE_ACTIVITY                     -0.195
# 9 GOMF_HSP90_PROTEIN_BINDING                       -0.178
# 10 GOCC_MITOCHONDRIAL_MATRIX                        -0.176


# # For each interesting pathway, do a volcano plot.

info <- enrichment_age %>%
  select(Pathway_name, adjusted_p) %>%
  rename(pathway = Pathway_name,
         adjusted_p_pathway = adjusted_p) %>%
  arrange(adjusted_p_pathway)

result_age %>%
  mutate(log10_p = -log10(p)) %>% 
  filter(adjusted_p < 0.01) %>%
  arrange(log10_p)

# # A tibble: 206 x 7
# Gene_ID            Gene_symbol  Effect Effect_SE        p adjusted_p log10_p
# <chr>              <chr>         <dbl>     <dbl>    <dbl>      <dbl>   <dbl>
#   1 ENSMUSG00000076609 Igkc         1.21      0.334  0.000504    0.00991    3.30

# The smallest log10_p under the padjusted threshold of 0.01 is 3.3

merge <- list_symbols %>%
  right_join(result_age %>%
               mutate(log10_p = -log10(p))) %>%
  mutate(`FDR < 0.01` = ifelse(log10_p > 3.30,
                               TRUE,
                               FALSE)) %>%
  left_join(info) %>%
  arrange(adjusted_p_pathway) %>%
  mutate(pathway = factor(pathway, levels = info$pathway))

merge <- merge[!duplicated(merge$Gene_symbol),] %>% #Make sure that each gene
  # symbol is in only one category
  mutate(pathway = as.character(pathway)) %>%
  mutate(pathway = ifelse(
    !pathway %in% c("GOBP_CELL_ACTIVATION",
                    "GOBP_EXOCYTOSIS",
                    "GOBP_HUMORAL_IMMUNE_RESPONSE",
                    "GOBP_IMMUNE_EFFECTOR_PROCESS",
                    "GOBP_LEUKOCYTE_MEDIATED_IMMUNITY",
                    "GOBP_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN",
                    "GOBP_PROTEIN_FOLDING",
                    "GOMF_UNFOLDED_PROTEIN_BINDING",
                    "GOMF_HSP90_PROTEIN_BINDING",
                    "GOMF_HEAT_SHOCK_PROTEIN_BINDING"),
    "Other or no pathway",
    pathway)) %>% 
  mutate(adjusted_p_pathway = ifelse(
    !pathway %in% c("GOBP_CELL_ACTIVATION",
                    "GOBP_EXOCYTOSIS",
                    "GOBP_HUMORAL_IMMUNE_RESPONSE",
                    "GOBP_IMMUNE_EFFECTOR_PROCESS",
                    "GOBP_LEUKOCYTE_MEDIATED_IMMUNITY",
                    "GOBP_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN",
                    "GOBP_PROTEIN_FOLDING",
                    "GOMF_UNFOLDED_PROTEIN_BINDING",
                    "GOMF_HSP90_PROTEIN_BINDING",
                    "GOMF_HEAT_SHOCK_PROTEIN_BINDING"),
    NA,
    adjusted_p_pathway))

ORDER <- merge %>% 
  select(pathway, adjusted_p_pathway) %>% 
  unique() %>% 
  arrange(adjusted_p_pathway)

merge <- merge %>% 
  mutate(pathway = factor(pathway, levels = ORDER$pathway))

age_rna_plot <- ggplot(data = merge %>% filter(pathway == "Other or no pathway"),
       aes(x = Effect, y = log10_p)) +
  geom_point(size = 3, alpha = 0.5, color = "grey") +
  geom_point(data = merge %>% filter(pathway %in%
                                       c("GOBP_CELL_ACTIVATION",
                                         "GOBP_EXOCYTOSIS",
                                         "GOBP_HUMORAL_IMMUNE_RESPONSE",
                                         "GOBP_IMMUNE_EFFECTOR_PROCESS",
                                         "GOBP_LEUKOCYTE_MEDIATED_IMMUNITY",
                                         "GOBP_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN",
                                         "GOBP_PROTEIN_FOLDING",
                                         "GOMF_UNFOLDED_PROTEIN_BINDING",
                                         "GOMF_HSP90_PROTEIN_BINDING",
                                         "GOMF_HEAT_SHOCK_PROTEIN_BINDING")),
             aes(color = pathway), alpha = 0.75, size = 3) +
  scale_color_manual(
    values = rev(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
                   '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a'))) +
  geom_vline(xintercept = 0, col="black", linetype = "dashed", size = 1) +
  geom_hline(yintercept = 3.3, col="red", linetype = "dashed", size = 1) +
  geom_text(data = merge %>% filter(pathway == "Other or no pathway"),
            aes(label = ifelse(
              log10_p > 7.5 |
                log10_p > 3.3 & Effect < -0.3 |
                log10_p > 5 & Effect > 0.5,
              Gene_symbol,'')), hjust = 0,vjust = 0, size = 4) +
  geom_text(data = merge %>% filter(pathway %in%
                                      c("GOBP_CELL_ACTIVATION",
                                        "GOBP_EXOCYTOSIS",
                                        "GOBP_HUMORAL_IMMUNE_RESPONSE",
                                        "GOBP_IMMUNE_EFFECTOR_PROCESS",
                                        "GOBP_LEUKOCYTE_MEDIATED_IMMUNITY",
                                        "GOBP_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN",
                                        "GOBP_PROTEIN_FOLDING",
                                        "GOMF_UNFOLDED_PROTEIN_BINDING",
                                        "GOMF_HSP90_PROTEIN_BINDING",
                                        "GOMF_HEAT_SHOCK_PROTEIN_BINDING")),
            aes(label = ifelse(
              log10_p > 7.5 |
                log10_p > 3.3 & Effect < -0.3 |
                log10_p > 5 & Effect > 0.5,
              Gene_symbol,'')), hjust = 0,vjust = 0, size = 5.5) +
  labs(x = "Age_effect",
       y = "-log10(p)") +
  theme_bw() +
  xlim(-3,3) +
  theme(legend.text = element_text(size = 2)) +
  ggtitle("Transcripts")

rm(enrichment_age, info, list_symbols, merge, ORDER, enrichment, result_age)

graphics.off()

# ------------------------------ Proteins

enrichment_age <- enrich_list$age_protein$collapsed_enrichment$leadingEdge

names(enrichment_age) <-  enrich_list$age_protein$collapsed_enrichment$pathway

temp <- plyr::ldply(enrichment_age, rbind) %>% 
  rename(pathway = .id) %>% 
  gather("test","ensembl_gene", -pathway) %>% 
  dplyr::select(-test) %>% 
  filter(!is.na(ensembl_gene)) %>% 
  left_join(all_gene_sets %>% select(ensembl_gene, gene_symbol)) %>% 
  unique() 

list_symbols <- list()

for(i in unique(temp$pathway)) {
  
  list <- temp %>% filter(pathway == i) %>% 
    select(gene_symbol) %>% 
    unlist()
  
  list <- knitr::combine_words(
    as.vector(list),
    sep = "/",
    and = "")
  
  list_symbols[[i]] <- list
  
}

rm(list, i, temp)

list_symbols <- bind_rows(list_symbols, .id = "pathway") %>% 
  gather("pathway",symbols_list)

enrichment_age <- enrich_list$age_protein$collapsed_enrichment %>% 
  left_join(list_symbols)

rm(list_symbols)

# Looking for the most significant and different pathways. 

enrichment_age %>% 
  arrange(padj) %>% 
  View()

# Saving enrichment results for supplemental data 

enrichment_age <- enrichment_age %>% 
  data.frame() %>% 
  dplyr::rename(Pathway_name = pathway,
                p = pval,
                adjusted_p = padj,
                Log2Err = log2err,
                Size = size,
                Gene_IDs = leadingEdge,
                Gene_symbols = symbols_list) %>% 
  arrange(adjusted_p) %>% 
  select(Pathway_name, Log2Err, ES, NES, Size, p, adjusted_p, Gene_IDs, Gene_symbols)

# data.table::fwrite(
#   enrichment_age,
#   file = "Downstream_Analysis/Results/enrichment_age_FULL_proteins.csv",
#   sep = ",",
#   na = "NA")

data.table::fwrite(
  enrichment_age,
  file = "Downstream_Analysis/Results/enrichment_age_proteins.csv",
  sep = ",",
  na = "NA")

# Now, merging enrichment with DEA results.

result_age <- read_csv("Downstream_Analysis/Results/DEA_age_proteins.csv")

enrichment <- enrichment_age$Gene_IDs

names(enrichment) <-  enrichment_age$Pathway_name

temp <- plyr::ldply(enrichment, rbind) %>%
  rename(pathway = .id) %>%
  gather("test","ensembl_gene", -pathway) %>%
  dplyr::select(-test) %>%
  filter(!is.na(ensembl_gene)) %>%
  left_join(all_gene_sets %>% select(ensembl_gene, gene_symbol)) %>%
  unique()

list_symbols <- list()

for(i in unique(temp$pathway)) {
  
  list <- temp %>% filter(pathway == i) %>%
    select(gene_symbol) %>%
    unlist()
  
  names(list) <- NULL
  
  list_symbols[[i]] <- list
  
}

rm(list, i, temp)

list_symbols <- plyr::ldply(list_symbols, rbind) %>%
  gather("rm","Gene_symbol", -`.id`) %>%
  rename(pathway = .id) %>%
  select(-rm) %>%
  filter(!is.na(Gene_symbol))

merge <- list_symbols %>%
  inner_join(result_age %>%
               filter(adjusted_p < 0.01))

# What are the pathways with the largest number of significant genes?

merge %>% group_by(pathway) %>% count() %>% arrange(desc(n))
# # A tibble: 26 x 2
# # Groups:   pathway [26]
# pathway                                     n
# <chr>                                   <int>
#   1 GOCC_ORGANELLE_SUBCOMPARTMENT             221
# 2 GOCC_MEMBRANE_PROTEIN_COMPLEX             167
# 3 GOBP_RESPONSE_TO_NITROGEN_COMPOUND        166
# 4 GOCC_GOLGI_APPARATUS                      160
# 5 GOCC_ENDOSOME                             156
# 6 GOBP_CELLULAR_PROTEIN_CATABOLIC_PROCESS   136
# 7 GOCC_VESICLE_MEMBRANE                     119
# 8 GOBP_ENDOMEMBRANE_SYSTEM_ORGANIZATION     108
# 9 GOCC_VACUOLAR_MEMBRANE                     80
# 10 GOMF_TRANSPORTER_ACTIVITY                  70
# # … with 16 more rows

# What are pathways with most negative mean effect? 

merge %>% group_by(pathway) %>% summarise(mean = mean(Effect)) %>% arrange(mean)
# A tibble: 26 x 2
# pathway                                                 mean
# <chr>                                                  <dbl>
#   1 GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY -0.407
# 2 GOMF_STRUCTURAL_CONSTITUENT_OF_MUSCLE                 -0.397
# 3 GOMF_OXIDOREDUCTASE_ACTIVITY_ACTING_ON_NAD_P_H        -0.328
# 4 GOMF_UNFOLDED_PROTEIN_BINDING                         -0.184
# 5 GOBP_FATTY_ACID_BETA_OXIDATION                        -0.163
# 6 GOCC_AP_TYPE_MEMBRANE_COAT_ADAPTOR_COMPLEX             0.222
# 7 GOCC_CLATHRIN_COAT                                     0.231
# 8 GOCC_MEMBRANE_COAT                                     0.248
# 9 GOBP_REGULATION_OF_AUTOPHAGY                           0.285
# 10 GOBP_ENDOMEMBRANE_SYSTEM_ORGANIZATION                  0.289
# # … with 16 more rows

# # For each interesting pathway, do a volcano plot.

info <- enrichment_age %>%
  select(Pathway_name, adjusted_p) %>%
  rename(pathway = Pathway_name,
         adjusted_p_pathway = adjusted_p) %>%
  arrange(adjusted_p_pathway)

result_age %>%
  mutate(log10_p = -log10(p)) %>% 
  filter(adjusted_p < 0.01) %>%
  arrange(log10_p)

# # A tibble: 2,084 x 8
# Protein_ID         Gene_ID            Gene_symbol  Effect Effect_SE       p adjusted_p log10_p
# <chr>              <chr>              <chr>         <dbl>     <dbl>   <dbl>      <dbl>   <dbl>
#   1 ENSMUSP00000025904 ENSMUSG00000024953 Prdx5       -0.110     0.0388 0.00504    0.00996    2.30
# The smallest log10_p under the padjusted threshold of 0.01 is 2.3

merge <- list_symbols %>%
  right_join(result_age %>%
               mutate(log10_p = -log10(p))) %>%
  mutate(`FDR < 0.01` = ifelse(log10_p > 2.30,
                               TRUE,
                               FALSE)) %>%
  left_join(info) %>%
  arrange(adjusted_p_pathway) %>%
  mutate(pathway = factor(pathway, levels = info$pathway))

merge <- merge[!duplicated(merge$Gene_symbol),] %>% #Make sure that each gene
  # symbol is in only one category
  mutate(pathway = as.character(pathway)) %>%
  mutate(pathway = ifelse(
    !pathway %in% c("GOCC_ORGANELLE_SUBCOMPARTMENT",
                    "GOCC_VACUOLAR_MEMBRANE",
                    "GOCC_GOLGI_APPARATUS",
                    "GOCC_VESICLE_MEMBRANE",
                    "GOBP_REGULATION_OF_AUTOPHAGY",
                    "GOBP_FATTY_ACID_BETA_OXIDATION",
                    "GOMF_UNFOLDED_PROTEIN_BINDING",
                    "GOMF_STRUCTURAL_CONSTITUENT_OF_MUSCLE",
                    "GOBP_DE_NOVO_PROTEIN_FOLDING",
                    "GOBP_GLUCOSE_CATABOLIC_PROCESS"),
    "Other or no pathway",
    pathway)) %>% 
  mutate(adjusted_p_pathway = ifelse(
    !pathway %in% c("GOCC_ORGANELLE_SUBCOMPARTMENT",
                    "GOCC_VACUOLAR_MEMBRANE",
                    "GOCC_GOLGI_APPARATUS",
                    "GOCC_VESICLE_MEMBRANE",
                    "GOBP_REGULATION_OF_AUTOPHAGY",
                    "GOBP_FATTY_ACID_BETA_OXIDATION",
                    "GOMF_UNFOLDED_PROTEIN_BINDING",
                    "GOMF_STRUCTURAL_CONSTITUENT_OF_MUSCLE",
                    "GOBP_DE_NOVO_PROTEIN_FOLDING",
                    "GOBP_GLUCOSE_CATABOLIC_PROCESS"),
    NA,
    adjusted_p_pathway))

# ORDER <- merge %>% 
#   select(pathway, adjusted_p_pathway) %>% 
#   unique() %>% 
#   arrange(adjusted_p_pathway)

merge <- merge %>% 
  mutate(pathway = factor(pathway, levels = c("GOCC_ORGANELLE_SUBCOMPARTMENT",
                                              "GOCC_VACUOLAR_MEMBRANE",
                                              "GOCC_GOLGI_APPARATUS",
                                              "GOCC_VESICLE_MEMBRANE",
                                              "GOBP_REGULATION_OF_AUTOPHAGY",
                                              "GOBP_FATTY_ACID_BETA_OXIDATION",
                                              "GOMF_UNFOLDED_PROTEIN_BINDING",
                                              "GOMF_STRUCTURAL_CONSTITUENT_OF_MUSCLE",
                                              "GOBP_DE_NOVO_PROTEIN_FOLDING",
                                              "GOBP_GLUCOSE_CATABOLIC_PROCESS",
                                              "Other or no pathway")))

age_protein_plot <- ggplot(data = merge %>% filter(pathway == "Other or no pathway"),
                       aes(x = Effect, y = log10_p)) +
  geom_point(size = 3, alpha = 0.5, color = "grey") +
  geom_point(data = merge %>% filter(pathway %in%
                                       c("GOCC_ORGANELLE_SUBCOMPARTMENT",
                                         "GOCC_VACUOLAR_MEMBRANE",
                                         "GOCC_GOLGI_APPARATUS",
                                         "GOCC_VESICLE_MEMBRANE",
                                         "GOBP_REGULATION_OF_AUTOPHAGY",
                                         "GOBP_FATTY_ACID_BETA_OXIDATION",
                                         "GOMF_UNFOLDED_PROTEIN_BINDING",
                                         "GOMF_STRUCTURAL_CONSTITUENT_OF_MUSCLE",
                                         "GOBP_DE_NOVO_PROTEIN_FOLDING",
                                         "GOBP_GLUCOSE_CATABOLIC_PROCESS")),
             aes(color = pathway), alpha = 0.75, size = 3) +
  scale_color_manual(
    values = rev(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
                   '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a'))) +
  geom_vline(xintercept = 0, col="black", linetype = "dashed", size = 1) +
  geom_hline(yintercept = 3.3, col="red", linetype = "dashed", size = 1) +
  geom_text(data = merge %>% filter(pathway == "Other or no pathway"),
            aes(label = ifelse(
              log10_p > 30 |
                log10_p > 2.3 & Effect > 1.5|
                log10_p > 2.3 & Effect < -0.8,
              Gene_symbol,'')), hjust = 0,vjust = 0, size = 4) +
  geom_text(data = merge %>% filter(pathway %in%
                                      c("GOCC_ORGANELLE_SUBCOMPARTMENT",
                                        "GOCC_VACUOLAR_MEMBRANE",
                                        "GOCC_GOLGI_APPARATUS",
                                        "GOCC_VESICLE_MEMBRANE",
                                        "GOBP_REGULATION_OF_AUTOPHAGY",
                                        "GOBP_FATTY_ACID_BETA_OXIDATION",
                                        "GOMF_UNFOLDED_PROTEIN_BINDING",
                                        "GOMF_STRUCTURAL_CONSTITUENT_OF_MUSCLE",
                                        "GOBP_DE_NOVO_PROTEIN_FOLDING",
                                        "GOBP_GLUCOSE_CATABOLIC_PROCESS")),
            aes(label = ifelse(
              log10_p > 30 |
                log10_p > 2.3 & Effect > 1.5|
                log10_p > 2.3 & Effect < -0.8,
              Gene_symbol,'')), hjust = 0,vjust = 0, size = 5.5) +
  labs(x = "Age_effect",
       y = "-log10(p)") +
  theme_bw() +
  theme(legend.text = element_text(size = 2)) +
  xlim(-3.5,3.5) +
  ggtitle("Proteins")

pdf("Downstream_Analysis/Results/VolcanoPlot_enrichment_selection.pdf",
    width = 15, height = 5.5)

cowplot::plot_grid(age_rna_plot, age_protein_plot, ncol = 2)

dev.off()

rm(age_protein_plot, age_rna_plot, annot.samples, enrichment, enrichment_age,
   info, list_symbols, merge, ORDER, result_age, common_set, enrich_function)


################################################################################

# 2) Sex

# ------------------------------ Transcripts

enrichment_sex <- enrich_list$sex_rna$collapsed_enrichment$leadingEdge

names(enrichment_sex) <-  enrich_list$sex_rna$collapsed_enrichment$pathway

temp <- plyr::ldply(enrichment_sex, rbind) %>% 
  rename(pathway = .id) %>% 
  gather("test","ensembl_gene", -pathway) %>% 
  dplyr::select(-test) %>% 
  filter(!is.na(ensembl_gene)) %>% 
  left_join(all_gene_sets %>% select(ensembl_gene, gene_symbol)) %>% 
  unique() 

list_symbols <- list()

for(i in unique(temp$pathway)) {
  
  list <- temp %>% filter(pathway == i) %>% 
    select(gene_symbol) %>% 
    unlist()
  
  list <- knitr::combine_words(
    as.vector(list),
    sep = "/",
    and = "")
  
  list_symbols[[i]] <- list
  
}

rm(list, i, temp)

list_symbols <- bind_rows(list_symbols, .id = "pathway") %>% 
  gather("pathway",symbols_list)

enrichment_sex <- enrich_list$sex_rna$collapsed_enrichment %>% 
  left_join(list_symbols)

rm(list_symbols)

# Looking for the most significant and different pathways. 

enrichment_sex %>% 
  arrange(padj) %>% 
  View()

# Saving enrichment results for supplemental data 

enrichment_sex <- enrichment_sex %>% 
  data.frame() %>% 
  dplyr::rename(Pathway_name = pathway,
                p = pval,
                adjusted_p = padj,
                Log2Err = log2err,
                Size = size,
                Gene_IDs = leadingEdge,
                Gene_symbols = symbols_list) %>% 
  arrange(adjusted_p) %>% 
  select(Pathway_name, Log2Err, ES, NES, Size, p, adjusted_p, Gene_IDs, Gene_symbols)

# data.table::fwrite(
#   enrichment_sex,
#   file = "Downstream_Analysis/Results/enrichment_sex_FULL_transcripts.csv",
#   sep = ",",
#   na = "NA")

# data.table::fwrite(
#   enrichment_sex,
#   file = "Downstream_Analysis/Results/enrichment_sex_transcripts.csv",
#   sep = ",",
#   na = "NA")

rm(enrichment_sex)

# ------------------------------ Proteins

enrichment_sex <- enrich_list$sex_protein$collapsed_enrichment$leadingEdge

names(enrichment_sex) <-  enrich_list$sex_protein$collapsed_enrichment$pathway

temp <- plyr::ldply(enrichment_sex, rbind) %>% 
  rename(pathway = .id) %>% 
  gather("test","ensembl_gene", -pathway) %>% 
  dplyr::select(-test) %>% 
  filter(!is.na(ensembl_gene)) %>% 
  left_join(all_gene_sets %>% select(ensembl_gene, gene_symbol)) %>% 
  unique() 

list_symbols <- list()

for(i in unique(temp$pathway)) {
  
  list <- temp %>% filter(pathway == i) %>% 
    select(gene_symbol) %>% 
    unlist()
  
  list <- knitr::combine_words(
    as.vector(list),
    sep = "/",
    and = "")
  
  list_symbols[[i]] <- list
  
}

rm(list, i, temp)

list_symbols <- bind_rows(list_symbols, .id = "pathway") %>% 
  gather("pathway",symbols_list)

enrichment_sex <- enrich_list$sex_protein$collapsed_enrichment %>% 
  left_join(list_symbols)

rm(list_symbols)

# Looking for the most significant and different pathways. 

enrichment_sex %>% 
  arrange(padj) %>% 
  View()

# Saving enrichment results for supplemental data 

enrichment_sex <- enrichment_sex %>% 
  data.frame() %>% 
  dplyr::rename(Pathway_name = pathway,
                p = pval,
                adjusted_p = padj,
                Log2Err = log2err,
                Size = size,
                Gene_IDs = leadingEdge,
                Gene_symbols = symbols_list) %>% 
  arrange(adjusted_p) %>% 
  select(Pathway_name, Log2Err, ES, NES, Size, p, adjusted_p, Gene_IDs, Gene_symbols)

# data.table::fwrite(
#   enrichment_sex,
#   file = "Downstream_Analysis/Results/enrichment_sex_FULL_proteins.csv",
#   sep = ",",
#   na = "NA")

# data.table::fwrite(
#   enrichment_sex,
#   file = "Downstream_Analysis/Results/enrichment_sex_proteins.csv",
#   sep = ",",
#   na = "NA")

rm(result_sex, enrichment_sex)


################################################################################

# 3) Sex by age interaction

# ------------------------------ Transcripts

enrichment_int <- enrich_list$int_rna$collapsed_enrichment$leadingEdge

names(enrichment_int) <-  enrich_list$int_rna$collapsed_enrichment$pathway

temp <- plyr::ldply(enrichment_int, rbind) %>% 
  rename(pathway = .id) %>% 
  gather("test","ensembl_gene", -pathway) %>% 
  dplyr::select(-test) %>% 
  filter(!is.na(ensembl_gene)) %>% 
  left_join(all_gene_sets %>% select(ensembl_gene, gene_symbol)) %>% 
  unique() 

list_symbols <- list()

for(i in unique(temp$pathway)) {
  
  list <- temp %>% filter(pathway == i) %>% 
    select(gene_symbol) %>% 
    unlist()
  
  list <- knitr::combine_words(
    as.vector(list),
    sep = "/",
    and = "")
  
  list_symbols[[i]] <- list
  
}

rm(list, i, temp)

list_symbols <- bind_rows(list_symbols, .id = "pathway") %>% 
  gather("pathway",symbols_list)

enrichment_int <- enrich_list$int_rna$collapsed_enrichment %>% 
  left_join(list_symbols)

rm(list_symbols)

# Looking for the most significant and different pathways. 

enrichment_int %>% 
  arrange(padj) %>% 
  View()

# Saving enrichment results for supplemental data 

enrichment_int <- enrichment_int %>% 
  data.frame() %>% 
  dplyr::rename(Pathway_name = pathway,
                p = pval,
                adjusted_p = padj,
                Log2Err = log2err,
                Size = size,
                Gene_IDs = leadingEdge,
                Gene_symbols = symbols_list) %>% 
  arrange(adjusted_p) %>% 
  select(Pathway_name, Log2Err, ES, NES, Size, p, adjusted_p, Gene_IDs, Gene_symbols)

# data.table::fwrite(
#   enrichment_int,
#   file = "Downstream_Analysis/Results/enrichment_int_FULL_transcripts.csv",
#   sep = ",",
#   na = "NA")

# data.table::fwrite(
#   enrichment_int,
#   file = "Downstream_Analysis/Results/enrichment_int_transcripts.csv",
#   sep = ",",
#   na = "NA")

rm(enrichment_int)

# ------------------------------ Proteins

enrichment_int <- enrich_list$int_protein$collapsed_enrichment$leadingEdge

names(enrichment_int) <-  enrich_list$int_protein$collapsed_enrichment$pathway

temp <- plyr::ldply(enrichment_int, rbind) %>% 
  rename(pathway = .id) %>% 
  gather("test","ensembl_gene", -pathway) %>% 
  dplyr::select(-test) %>% 
  filter(!is.na(ensembl_gene)) %>% 
  left_join(all_gene_sets %>% select(ensembl_gene, gene_symbol)) %>% 
  unique() 

list_symbols <- list()

for(i in unique(temp$pathway)) {
  
  list <- temp %>% filter(pathway == i) %>% 
    select(gene_symbol) %>% 
    unlist()
  
  list <- knitr::combine_words(
    as.vector(list),
    sep = "/",
    and = "")
  
  list_symbols[[i]] <- list
  
}

rm(list, i, temp)

list_symbols <- bind_rows(list_symbols, .id = "pathway") %>% 
  gather("pathway",symbols_list)

enrichment_int <- enrich_list$int_protein$collapsed_enrichment %>% 
  left_join(list_symbols)

rm(list_symbols)

# Looking for the most significant and different pathways. 

enrichment_int %>% 
  arrange(padj) %>% 
  View()

# Saving enrichment results for supplemental data 

enrichment_int <- enrichment_int %>% 
  data.frame() %>% 
  dplyr::rename(Pathway_name = pathway,
                p = pval,
                adjusted_p = padj,
                Log2Err = log2err,
                Size = size,
                Gene_IDs = leadingEdge,
                Gene_symbols = symbols_list) %>% 
  arrange(adjusted_p) %>% 
  select(Pathway_name, Log2Err, ES, NES, Size, p, adjusted_p, Gene_IDs, Gene_symbols)

# data.table::fwrite(
#   enrichment_int,
#   file = "Downstream_Analysis/Results/enrichment_int_FULL_proteins.csv",
#   sep = ",",
#   na = "NA")

# data.table::fwrite(
#   enrichment_int,
#   file = "Downstream_Analysis/Results/enrichment_int_proteins.csv",
#   sep = ",",
#   na = "NA")

rm(result_int, enrichment_int)


# Now, making p-values plots for DEA age, sex and interaction for both transcripts
# and proteins. 

dea <- bind_rows(results, .id = "Type") %>% 
  separate(Type, c("Test","Type"), "_") %>% 
  mutate(Type = ifelse(Type == "rna","Transcript","Protein"),
         Test = ifelse(Test == "age","Age effect",
                       ifelse(Test == "sex", "Sex effect",
                              "Age by sex interaction effect"))) %>% 
  select(-Protein_ID) %>% 
  mutate(Type = factor(Type, 
                       levels = c("Transcript","Protein")),
         Test = factor(Test, 
                       levels = c("Age effect","Sex effect","Age by sex interaction effect")))

pdf("Downstream_Analysis/Results/p_value_histograms.pdf", 
    width = 8, height = 4)
dea %>% 
  ggplot(aes(x = p)) +
  geom_histogram(bins = 25) +
  facet_grid(Type ~ Test) +
  theme_bw()
dev.off()

# Checking the number of significant genes according to different thresholds to 
# make a table

dea %>% 
  group_by(Type, Test) %>% 
  filter(adjusted_p < 0.01) %>% 
  summarise(Count = length(Gene_ID))

dea %>% 
  group_by(Type, Test) %>% 
  filter(adjusted_p < 0.1) %>% 
  summarise(Count = length(Gene_ID))
  

# Now, gathering the common enrichments between protein and transcript for age

transcript <- read.csv("Downstream_Analysis/Results/enrichment_age_transcripts.csv") %>% 
  mutate(Gene_IDs = as.list(Gene_IDs))

protein <- read.csv("Downstream_Analysis/Results/enrichment_age_proteins.csv") %>% 
  mutate(Gene_IDs = as.list(Gene_IDs))

intersect(transcript$Pathway_name, protein$Pathway_name)
#[1] "GOMF_UNFOLDED_PROTEIN_BINDING"

# What about common genes but different pathway names?


# --------------------- For proteins

DEA_proteins <- read_csv("Downstream_Analysis/Results/DEA_age_proteins.csv")

enrichment <- protein$Gene_IDs

names(enrichment) <-  protein$Pathway_name

# Transforming these into character vector

for (i in names(enrichment)) {
  
  tmp <- gsub("[[:punct:]]"," ",enrichment[[i]])
  
  tmp <- str_split(tmp," ", simplify = TRUE) %>% as.vector()
  
  enrichment[[i]] <- tmp
}

rm(i, tmp)

temp <- plyr::ldply(enrichment, rbind) %>%
  rename(pathway = .id) %>%
  gather("test","ensembl_gene", -pathway) %>%
  dplyr::select(-test) %>%
  filter(!is.na(ensembl_gene)) %>%
  left_join(all_gene_sets %>% select(ensembl_gene, gene_symbol)) %>%
  unique()

list_symbols <- list()

for(i in unique(temp$pathway)) {
  
  list <- temp %>% filter(pathway == i) %>%
    select(gene_symbol) %>%
    unlist()
  
  names(list) <- NULL
  
  list_symbols[[i]] <- list
  
}

rm(list, i, temp)

list_symbols_protein <- plyr::ldply(list_symbols, rbind) %>%
  gather("rm","Gene_symbol", -`.id`) %>%
  rename(Pathway_name = .id) %>%
  select(-rm) %>%
  filter(!is.na(Gene_symbol)) %>% 
  left_join(protein %>% select(Pathway_name, adjusted_p)) %>% 
  arrange(adjusted_p)

rm(list_symbols, enrichment)

# For transcripts

DEA_transcript <- read_csv("Downstream_Analysis/Results/DEA_age_transcripts.csv")

enrichment <- transcript$Gene_IDs

names(enrichment) <-  transcript$Pathway_name

# Transforming these into character vector

for (i in names(enrichment)) {
  
  tmp <- gsub("[[:punct:]]"," ",enrichment[[i]])
  
  tmp <- str_split(tmp," ", simplify = TRUE) %>% as.vector()
  
  enrichment[[i]] <- tmp
}

rm(i, tmp)

temp <- plyr::ldply(enrichment, rbind) %>%
  rename(pathway = .id) %>%
  gather("test","ensembl_gene", -pathway) %>%
  dplyr::select(-test) %>%
  filter(!is.na(ensembl_gene)) %>%
  left_join(all_gene_sets %>% select(ensembl_gene, gene_symbol)) %>%
  unique()

list_symbols <- list()

for(i in unique(temp$pathway)) {
  
  list <- temp %>% filter(pathway == i) %>%
    select(gene_symbol) %>%
    unlist()
  
  names(list) <- NULL
  
  list_symbols[[i]] <- list
  
}

rm(list, i, temp)

list_symbols_transcript <- plyr::ldply(list_symbols, rbind) %>%
  gather("rm","Gene_symbol", -`.id`) %>%
  rename(Pathway_name = .id) %>%
  select(-rm) %>%
  filter(!is.na(Gene_symbol)) %>% 
  left_join(transcript %>% select(Pathway_name, adjusted_p)) %>% 
  arrange(adjusted_p)

rm(list_symbols, enrichment)

# Joining the datasets

list_symbols_transcript <- list_symbols_transcript %>% 
  mutate(Transcript_enrichment = TRUE) %>% 
  select(-adjusted_p)


list_symbols_protein <- list_symbols_protein %>% 
  mutate(Protein_enrichment = TRUE) %>% 
  select(-adjusted_p)


merge <- list_symbols_transcript %>%
  full_join(list_symbols_protein, by = "Gene_symbol") %>% 
  rename(Pathway_name_transcript = Pathway_name.x,
         Pathway_name_protein = Pathway_name.y) %>% 
  select(Pathway_name_transcript, Pathway_name_protein, Gene_symbol,
         Transcript_enrichment, Protein_enrichment) %>% 
  mutate(Transcript_enrichment = ifelse(
    is.na(Transcript_enrichment), FALSE,Transcript_enrichment),
    Protein_enrichment = ifelse(
      is.na(Protein_enrichment), FALSE,Protein_enrichment),
  ) %>% 
  left_join(DEA_transcript %>%
              mutate(log10_p_transcript = -log10(p))) %>%
  mutate(Transcript_significant = ifelse(log10_p_transcript > 3.3,
                                         TRUE,
                                         FALSE)) %>% 
  rename(Effect_transcript = Effect,
         Effect_SE_transcript = Effect_SE) %>% 
  mutate(STD_Effect_transcript = Effect_transcript/Effect_SE_transcript) %>% 
  select(Pathway_name_transcript, Pathway_name_protein,Transcript_enrichment, 
         Protein_enrichment, Gene_ID, Gene_symbol, Effect_transcript, 
         Effect_SE_transcript, STD_Effect_transcript, log10_p_transcript, 
         Transcript_significant) %>% 
  left_join(DEA_proteins %>%
              mutate(log10_p_protein = -log10(p)), by = c("Gene_ID","Gene_symbol")) %>%
  mutate(Protein_significant = ifelse(log10_p_protein > 2.3,
                                      TRUE,
                                      FALSE)) %>% 
  rename(Effect_protein = Effect,
         Effect_SE_protein = Effect_SE) %>% 
  mutate(STD_Effect_protein = Effect_protein/Effect_SE_protein) %>% 
  select(Pathway_name_transcript, Pathway_name_protein, Transcript_enrichment, 
         Protein_enrichment, Gene_ID, Gene_symbol, Effect_transcript, 
         Effect_SE_transcript,STD_Effect_transcript, log10_p_transcript,
         Transcript_significant, Effect_protein, Effect_SE_protein,
         STD_Effect_protein, log10_p_protein,Protein_significant) %>% 
  unique()

merge %>% filter(Transcript_enrichment & Protein_enrichment) %>% 
  select(Gene_symbol) %>% unlist() %>% unique() %>% length()
# There are 250 that appear both on transcript and protein enrichments

merge %>% filter(Transcript_enrichment & Protein_enrichment) %>% 
  select(Pathway_name_transcript) %>% unlist() %>% unique() %>% length()
# They are spread over 79 categories for transcripts

merge %>% filter(Transcript_enrichment & Protein_enrichment) %>% 
  select(Pathway_name_protein) %>% unlist() %>% unique() %>% length()
# and 24 categories for proteins
  
# Checking the different directions of change and their pathways

merge %>% filter(Transcript_enrichment & Protein_enrichment, 
                 Effect_transcript > 0, Effect_protein > 0) %>% View()
# Transcripts: GOBP_CELL_ACTIVATION, GOBP_EXOCYTOSIS, GOCC_SECRETORY_GRANULE and
# GOBP_IMMUNE_EFFECTOR_PROCESS.
# Proteins: GOBP_ENDOMEMBRANE_SYSTEM_ORGANIZATION, GOCC_ORGANELLE_SUBCOMPARTMENT,
# GOCC_GOLGI_APPARATUS, GOBP_CELLULAR_PROTEIN_CATABOLIC_PROCESS.

merge %>% filter(Transcript_enrichment & Protein_enrichment, 
                 Effect_transcript > 0, Effect_protein < 0) %>% View()
# Transcripts: GOBP_CELL_ACTIVATION, GOBP_EXOCYTOSIS, GOCC_SECRETORY_GRANULE and
# GOBP_IMMUNE_EFFECTOR_PROCESS.
# Proteins: GOMF_OXIDOREDUCTASE_ACTIVITY_ACTING_ON_NAD_P_H, 
# GOMF_UNFOLDED_PROTEIN_BINDING

merge %>% filter(Transcript_enrichment & Protein_enrichment, 
                 Effect_transcript < 0, Effect_protein < 0) %>% View()
# Transcripts: GOBP_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN, GOCC_SARCOPLASM, 
# GOMF_UNFOLDED_PROTEIN_BINDING and GOBP_FATTY_ACID_CATABOLIC_PROCESS
# Proteins: GOMF_UNFOLDED_PROTEIN_BINDING, 
# GOMF_STRUCTURAL_CONSTITUENT_OF_MUSCLE, GOBP_FATTY_ACID_BETA_OXIDATION

merge %>% filter(Transcript_enrichment & Protein_enrichment, 
                 Effect_transcript < 0, Effect_protein > 0) %>% View()
# Transcripts: GOBP_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN, GOMF_PROTEIN_KINASE_ACTIVITY, 
# GOBP_CELLULAR_RESPONSE_TO_ALCOHOL 
# Proteins: GOCC_ORGANELLE_SUBCOMPARTMENT, 
# GOCC_GOLGI_APPARATUS, GOCC_ENDOSOME

pdf("Downstream_Analysis/Results/enrichment_overlap_scatterplot.pdf", 
    width = 10, height = 6)

merge %>% 
  select(Transcript_enrichment, Protein_enrichment, Gene_symbol, 
         STD_Effect_transcript, STD_Effect_protein, Transcript_significant, 
         Protein_significant) %>% 
  unique() %>% 
  mutate(Significance_level = ifelse(
    Transcript_significant & Protein_significant, 
    "Both",
    ifelse(!Transcript_significant & !Protein_significant,
           "Neither",
           ifelse(Transcript_significant & !Protein_significant, 
           "Transcript only",
                     "Protein only")))) %>% 
  filter(Transcript_enrichment & Protein_enrichment) %>% 
  ggplot(aes(STD_Effect_transcript, STD_Effect_protein)) +
  geom_point(size = 4, alpha = 0.5, aes(color = Significance_level)) +
  geom_vline(xintercept = 0, col="black", linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, col="black", linetype = "dashed", size = 1) +
  scale_colour_manual(breaks = c('Transcript only', 'Protein only', 'Both', 'Neither'),
                        values = c("firebrick","royalblue","mediumorchid4","gray")) +  
  geom_text(aes(label = ifelse(
    abs(STD_Effect_transcript) > 4 |
    STD_Effect_protein > 7.5 |
      STD_Effect_protein < -2,
    Gene_symbol,"")), hjust = 0, vjust = 0, size = 5) +
  scale_x_continuous(limits = c(-12,12)) +
  scale_y_continuous(limits = c(-12,12)) +
  theme_bw()

dev.off()


# Now, trying to plot two different volcano plots for transcripts and proteins

# merge_edit <- merge %>% 
#   filter(Transcript_enrichment & Protein_enrichment) %>%
#   mutate(Pathway_name_transcript = ifelse(
#     !Pathway_name_transcript %in% c("GOBP_CELL_ACTIVATION",
#                                     "GOBP_EXOCYTOSIS",
#                                     "GOBP_IMMUNE_EFFECTOR_PROCESS",
#                                     "GOBP_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN",
#                                     "GOCC_SARCOPLASM",
#                                     "GOBP_FATTY_ACID_CATABOLIC_PROCESS"),
#     "Other or no pathway",
#     Pathway_name_transcript)) %>% 
#   mutate(Pathway_name_protein = ifelse(
#     !Pathway_name_protein %in% c("GOBP_ENDOMEMBRANE_SYSTEM_ORGANIZATION",
#                                  "GOCC_GOLGI_APPARATUS",
#                                  "GOBP_CELLULAR_PROTEIN_CATABOLIC_PROCESS",
#                                  "GOMF_UNFOLDED_PROTEIN_BINDING",
#                                  "GOMF_STRUCTURAL_CONSTITUENT_OF_MUSCLE",
#                                  "GOBP_FATTY_ACID_BETA_OXIDATION"),
#     "Other or no pathway",
#     Pathway_name_protein)) %>% 
#   gather("Type_pathway", "pathway", 
#          Pathway_name_transcript, Pathway_name_protein) %>% 
#   mutate(Type_pathway = ifelse(
#     Type_pathway == "Pathway_name_transcript",
#     "Transcript",
#     "Protein")) %>% 
#   gather("Type_effect","Effect",Effect_transcript, Effect_protein) %>% 
#     mutate(Type_effect = ifelse(
#       Type_effect == "Effect_transcript",
#       "Transcript",
#       "Protein")) %>% 
#   gather("Type_log10_p","log10_p",log10_p_transcript, log10_p_protein) %>% 
#   mutate(Type_log10_p = ifelse(
#       Type_log10_p == "log10_p_transcript",
#       "Transcript",
#       "Protein")) %>% 
#   filter(Type_pathway == Type_effect &
#            Type_effect == Type_log10_p) %>% 
#   select(Type_pathway, pathway, Gene_symbol, Effect, log10_p) %>% 
#   rename(Type = Type_pathway) %>% 
#   unique()
# 
# 
# pdf("Downstream_Analysis/Results/enrichment_overlap_volcanoplot.pdf", 
#     width = 15, height = 5.5)
# 
# subset <- merge_edit %>% 
#   filter(Type == "Transcript") %>% 
#   mutate(pathway = factor(pathway, levels = c("GOBP_CELL_ACTIVATION",
#                                               "GOBP_EXOCYTOSIS",
#                                               "GOBP_IMMUNE_EFFECTOR_PROCESS",
#                                               "GOBP_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN",
#                                               "GOCC_SARCOPLASM",
#                                               "GOBP_FATTY_ACID_CATABOLIC_PROCESS",
#                                               "Other or no pathway"))) %>% 
#   arrange(pathway)
# 
# subset <- subset[!duplicated(subset$Gene_symbol),]
# 
# Transcript_plot <- ggplot(data = subset %>%  
#                             filter(pathway == "Other or no pathway"),
#                           aes(x = Effect, y = log10_p)) +
#   geom_point(size = 3, alpha = 0.5, color = "grey") +
#   geom_point(data = subset %>% 
#                        filter(pathway %in% c("GOBP_CELL_ACTIVATION",
#                                              "GOBP_EXOCYTOSIS",
#                                              "GOBP_IMMUNE_EFFECTOR_PROCESS",
#                                              "GOBP_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN",
#                                              "GOCC_SARCOPLASM",
#                                              "GOBP_FATTY_ACID_CATABOLIC_PROCESS")) %>% 
#                         mutate(pathway = factor(pathway, 
#                                                 levels = c("GOBP_CELL_ACTIVATION",
#                                                            "GOBP_EXOCYTOSIS",
#                                                            "GOBP_IMMUNE_EFFECTOR_PROCESS",
#                                                            "GOBP_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN",
#                                                            "GOCC_SARCOPLASM",
#                                                            "GOBP_FATTY_ACID_CATABOLIC_PROCESS"))),
#              aes(color = pathway), alpha = 0.75, size = 3) +
#   scale_color_manual(
#     values = c('#ae017e','#f768a1','#fbb4b9','#225ea8','#41b6c4','#a1dab4')) +
#   geom_vline(xintercept = 0, col="black", linetype = "dashed", size = 1) +
#   geom_hline(yintercept = 3.3, col="red", linetype = "dashed", size = 1) +
#   geom_text(data = subset %>% 
#                     filter(pathway == "Other or no pathway"),
#             aes(label = ifelse(
#               log10_p > 4.5 | 
#                 log10_p > 3.3 &
#                 abs(Effect) > 1,
#               Gene_symbol,'')), hjust = 0,vjust = 0, size = 4) +
#   geom_text(data = subset %>% 
#               filter(pathway %in% c("GOBP_CELL_ACTIVATION",
#                                      "GOBP_EXOCYTOSIS",
#                                      "GOBP_IMMUNE_EFFECTOR_PROCESS",
#                                      "GOBP_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN",
#                                      "GOCC_SARCOPLASM",
#                                      "GOBP_FATTY_ACID_CATABOLIC_PROCESS")),
#             aes(label = ifelse(
#               log10_p > 4.5 | 
#                 log10_p > 3.3 &
#                 abs(Effect) > 1,
#               Gene_symbol,'')), hjust = 0,vjust = 0, size = 5.5) +
#   labs(x = "Age_effect",
#        y = "-log10(p)") +
#   theme_bw() +
#   theme(legend.text = element_text(size = 2)) +
#   xlim(-2,2) +
#   ggtitle("Transcripts")
# 
# 
# subset <- merge_edit %>% 
#   filter(Type == "Protein") %>% 
#   mutate(pathway = factor(pathway, levels = c("GOBP_ENDOMEMBRANE_SYSTEM_ORGANIZATION",
#                                               "GOCC_GOLGI_APPARATUS",
#                                               "GOBP_CELLULAR_PROTEIN_CATABOLIC_PROCESS",
#                                               "GOMF_UNFOLDED_PROTEIN_BINDING",
#                                               "GOMF_STRUCTURAL_CONSTITUENT_OF_MUSCLE",
#                                               "GOBP_FATTY_ACID_BETA_OXIDATION",
#                                               "Other or no pathway"))) %>% 
#   arrange(pathway)
# 
# subset <- subset[!duplicated(subset$Gene_symbol),]
# 
# 
# Protein_plot <- ggplot(data = subset %>%  
#                          filter(Type == "Protein",
#                                 pathway == "Other or no pathway"),
#                        aes(x = Effect, y = log10_p)) +
#   geom_point(size = 3, alpha = 0.5, color = "grey") +
#   geom_point(data = subset %>% 
#                filter(Type == "Protein",
#                       pathway %in% c("GOBP_ENDOMEMBRANE_SYSTEM_ORGANIZATION",
#                                      "GOCC_GOLGI_APPARATUS",
#                                      "GOBP_CELLULAR_PROTEIN_CATABOLIC_PROCESS",
#                                      "GOMF_UNFOLDED_PROTEIN_BINDING",
#                                      "GOMF_STRUCTURAL_CONSTITUENT_OF_MUSCLE",
#                                      "GOBP_FATTY_ACID_BETA_OXIDATION")) %>% 
#                mutate(pathway = factor(pathway, 
#                                        levels = c("GOBP_ENDOMEMBRANE_SYSTEM_ORGANIZATION",
#                                                   "GOCC_GOLGI_APPARATUS",
#                                                   "GOBP_CELLULAR_PROTEIN_CATABOLIC_PROCESS",
#                                                   "GOMF_UNFOLDED_PROTEIN_BINDING",
#                                                   "GOMF_STRUCTURAL_CONSTITUENT_OF_MUSCLE",
#                                                   "GOBP_FATTY_ACID_BETA_OXIDATION"))),
#              aes(color = pathway), alpha = 0.75, size = 3) +
#   scale_color_manual(
#     values = c('#ae017e','#f768a1','#fbb4b9','#225ea8','#41b6c4','#a1dab4')) +
#   geom_vline(xintercept = 0, col="black", linetype = "dashed", size = 1) +
#   geom_hline(yintercept = 2.3, col="red", linetype = "dashed", size = 1) +
#   geom_text(data = subset %>% 
#               filter(pathway == "Other or no pathway"),
#             aes(label = ifelse(
#             log10_p > 15 | 
#               log10_p > 2.3 &
#               Effect > 0.7 |
#               log10_p > 2.3 &
#               Effect < -0.3,
#               Gene_symbol,'')), hjust = 0,vjust = 0, size = 4) +
#   geom_text(data = subset %>% 
#               filter(pathway %in% c("GOBP_ENDOMEMBRANE_SYSTEM_ORGANIZATION",
#                                     "GOCC_GOLGI_APPARATUS",
#                                     "GOBP_CELLULAR_PROTEIN_CATABOLIC_PROCESS",
#                                     "GOMF_UNFOLDED_PROTEIN_BINDING",
#                                     "GOMF_STRUCTURAL_CONSTITUENT_OF_MUSCLE",
#                                     "GOBP_FATTY_ACID_BETA_OXIDATION")),
#             aes(label = ifelse(
#               log10_p > 15 | 
#                 log10_p > 2.3 &
#                 Effect > 0.7 |
#                 log10_p > 2.3 &
#                 Effect < -0.3,
#               Gene_symbol,'')), hjust = 0,vjust = 0, size = 5.5) +
#   labs(x = "Age_effect",
#        y = "-log10(p)") +
#   theme_bw() +
#   theme(legend.text = element_text(size = 2)) +
#   xlim(-1.5,1.5) +
#   ggtitle("Proteins")
# 
# 
# cowplot::plot_grid(Transcript_plot, Protein_plot, ncol = 2)
# 
# 
# dev.off()

# Using something similar to dotplot from clusterProfiler to plot a summary
# of the enrichments. ince we have a lot of categories, selecting the top 5 that
# are up-regulated and top 5 down regulated.

transcript <- enrich_list$age_rna$collapsed_enrichment %>% 
                filter(ES > 0) %>% 
                 arrange(padj)

transcript_up <- transcript[1:5,]

transcript <- enrich_list$age_rna$collapsed_enrichment %>% 
                filter(ES < 0) %>% 
                arrange(padj)

transcript_down <- transcript[c(1:4,6),] # The 5th category has the same exact genes 
# as the GOMF_HSP90_PROTEIN_BINDING.

transcript <- rbind(
  transcript_up, transcript_down
)

ORDER <- transcript$pathway

plot1 <- transcript %>% 
  mutate(pathway = factor(pathway, levels = rev(ORDER))) %>% 
  ggplot(aes(x = ES, y = pathway)) +
  geom_point(aes(color = padj, size = size)) +
  scale_color_gradientn(colours = rev(c('#bdc9e1','#67a9cf','#1c9099','#016c59'))) +
  scale_x_continuous(limits = c(-0.7,0.7), breaks = c(-0.7, -0.35, 0, 0.35, 0.7)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  theme_bw() +
  theme(axis.text = element_text(size = 1))


protein <- enrich_list$age_protein$collapsed_enrichment %>% 
  filter(ES > 0) %>% 
  arrange(padj)

protein_up <- protein[1:5,]

protein <- enrich_list$age_protein$collapsed_enrichment %>% 
  filter(ES < 0) %>% 
  arrange(padj)

protein_down <- protein[1:5,]

protein <- rbind(
  protein_up, protein_down
)

ORDER <- protein$pathway

plot2 <- protein %>% 
  mutate(pathway = factor(pathway, levels = rev(ORDER))) %>% 
  ggplot(aes(x = ES, y = pathway)) +
  geom_point(aes(color = padj, size = size)) +
  scale_color_gradientn(colours = rev(c('#bdc9e1','#67a9cf','#1c9099','#016c59'))) +
  scale_x_continuous(limits = c(-0.6,0.6), breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  theme_bw() +
  theme(axis.text = element_text(size = 1))


pdf("Downstream_Analysis/Results/enrichment_dotplot.pdf", width = 12, height = 4)

cowplot::plot_grid(plot1,plot2, ncol = 2)

dev.off()

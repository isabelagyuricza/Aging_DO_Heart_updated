####################################################################################################################
#
#   This script contains functions to debatch protein and peptide data using LMER
#
#  Author: Greg Keele
#  Date:   March 4, 2020
#  Modified: 
#  E-mail: greg.keele@jax.org
#
####################################################################################################################
#### Function to regress out batch (and possibly sex)
debatch_protein_lmer <- function(long_dat, 
                                 filter = 0.5,
                                 filter_value = c("NA", "0"),
                                 regress_out_sex = FALSE) {
  filter_value <- filter_value[1]
  
  wide_dat <- long_dat %>%
    dplyr::select(protein.id, Intensity, Sex, Age, Generation, mouse.id, Batch, tag) %>%
    spread(key = protein.id, value = Intensity)
  
  ## Filter out proteins with lots of NAs
  if (filter_value == "NA") {
    na_prop <- wide_dat %>%
      dplyr::select(contains("ENSMUS")) %>%
      as.matrix %>%
      apply(2, function(x) mean(is.na(x)))
    proteins <- names(na_prop)[na_prop <= filter]
  }
  ## Filter out proteins with lots of 0s
  if (filter_value == "0") {
    zero_prop <- wide_dat %>%
      dplyr::select(contains("ENSMUS")) %>%
      as.matrix %>%
      apply(2, function(x) mean(x == 0))
    proteins <- names(zero_prop)[zero_prop <= filter]
  }

  resid_dat <- matrix(NA, nrow = nrow(wide_dat), ncol = length(proteins))
  rownames(resid_dat) <- wide_dat$mouse.id
  colnames(resid_dat) <- proteins
  for (i in 1:length(proteins)) {
    print(paste(i, "out of", length(proteins)))
      
    lmer_fit <- lme4::lmer(formula(paste(proteins[i], "~ Sex + Age + (1 | Batch)")), data = wide_dat)

    resid_dat[,i] <- wide_dat[,proteins[i], drop = TRUE] - lme4::ranef(lmer_fit)$Batch[wide_dat$Batch, 1]
    if (regress_out_sex) {
      resid_dat[,i] <- resid_dat[,proteins[i]] - lme4::fixef(lmer_fit)["SexM"] * model.matrix(~Sex, data = wide_dat)[, -1, drop = FALSE]
    }
  }

  resid_dat <- resid_dat %>% 
    as.data.frame %>%
    rownames_to_column("mouse.id") %>%
    left_join(wide_dat %>%
                dplyr::select(mouse.id, Sex, Age, Generation, Batch, tag))

  resid_dat
}


#################################################################################################################
#################################################################################################################
######
######                Title: protein_rollup.R
######                Description: Functions for roll up protein abundance estimates from
######                             their component peptides
######                
######                Manuscript: Keele & Zhang et al. 
######
######                Author: Greg Keele
######
#################################################################################################################
#################################################################################################################

## Timmed mean function
## Default does not include the minimum or maximum value
trimmed_mean <- function(x, 
                         trim = 1) {
  sort_x <- sort(x)
  if (trim > 0) {
    use_values <- sort_x[-c(1:trim, (length(x) - trim + 1):length(x))]
  }
  else {
    use_values <- sort_x
  }  
  trimmed_mean <- mean(use_values)
  trimmed_mean
}

## Expects batch and tag to named: Batch, tag
protein_rollup <- function(long_peptide_dat,
                           bridge_sample = NULL,
                           scale_to_100 = FALSE,
                           scale_trim_mean = FALSE,
                           trim = 1,
                           filter_out_peptides = NULL,
                           na_as_zero = FALSE,
                           zero_as_na = FALSE,
                           batch_scaling = TRUE) {
  
  ## Derive normalization factors within batch
  ## Ratio of max cumulative intensities for a tag to each cumulative tag intensity 
  if (batch_scaling) {
    tmt_max <- long_peptide_dat %>%
      group_by(Batch, tag) %>%
      summarise(sum_inten = sum(Intensity, na.rm = TRUE)) %>%
      ungroup %>%
      group_by(Batch) %>%
      summarise(max_sum_inten = max(sum_inten))
    tmt_norm <- long_peptide_dat %>%
      group_by(sample.id) %>%
      summarise(sum_inten = sum(Intensity, na.rm = TRUE)) %>%
      ungroup %>%
      mutate(Batch = gsub(x = sample.id, pattern = "~.*$", replacement = "", perl = TRUE)) %>%
      left_join(tmt_max) %>%
      mutate(norm_factor = max_sum_inten/sum_inten) %>%
      dplyr::select(-c(sum_inten, max_sum_inten, Batch))
  }

  ## Filtering out peptides
  if (!is.null(filter_out_peptides)) {
    long_peptide_dat <- long_peptide_dat %>%
      filter(!peptide.id %in% filter_out_peptides)
  }

  ## Sum peptides to quantify a protein
  long_protein_dat <- long_peptide_dat %>%
    ## Summing peptides for a protein
    group_by(protein.id, sample.id) %>%
    summarize(Intensity = sum(Intensity)) %>%
    ungroup %>%
    # Remove to add NAs
    complete(protein.id, sample.id) %>%
    # Add back Batch
    mutate(Batch = gsub(x = sample.id, pattern = "~.*$", replacement = "", perl = TRUE),
           tag = gsub(x = sample.id, pattern = "^[A-Z0-9]*~", replacement = "", perl = TRUE))
  
  ## Turn NAs into 0s
  if (na_as_zero) {
    long_protein_dat <- long_protein_dat %>%
      mutate(Intensity = ifelse(is.na(Intensity), 0, Intensity))
  }
  ## Turn 0s into NAs
  if (zero_as_na) {
    long_protein_dat <- long_protein_dat %>%
      mutate(Intensity = ifelse(Intensity == 0, NA, Intensity))
  }
  
  if (batch_scaling) {
    long_protein_dat <- long_protein_dat %>%
      ## Merge in tag-specific normalizations
      left_join(tmt_norm) %>%
      # Standardize across samples in a Batch
      mutate(Intensity = Intensity * norm_factor) %>%
      # Remove norm_factor
      dplyr::select(-norm_factor)
  }
  
  ## Scale to sum to 100 option
  if (scale_to_100) {
    scale100_dat <- long_protein_dat %>%
      group_by(protein.id, Batch) %>%
      summarize(scale = 100/sum(Intensity)) %>%
      ungroup %>%
      mutate(scale = ifelse(is.infinite(scale), 1, scale))
    long_protein_dat <- long_protein_dat %>%
      left_join(scale100_dat) %>%
      mutate(Intensity = Intensity * scale) %>%
      dplyr::select(protein.id, sample.id, Intensity, Batch, tag, -scale)
  }
  
  ## Scale by bridge sample option
  if (!is.null(bridge_sample)) {
    ## Shift 0s to 0.01 to avoid issues in log scale
    long_protein_dat <- long_protein_dat %>%
      #mutate(Intensity = ifelse(Intensity == 0, 0.01, Intensity))
      mutate(Intensity = log2(Intensity + 1))
    ## Pull out bridge samples
    bridge_dat <- long_protein_dat %>%
      filter(tag == bridge_sample) %>%
      dplyr::select(-c(sample.id, tag))
    long_protein_dat <- long_protein_dat %>%
      filter(tag != bridge_sample) %>%
      # Merging in bridge intensity
      left_join(bridge_dat %>%
                  dplyr::rename(bridge_intensity = Intensity)) %>%
      mutate(Intensity = Intensity - bridge_intensity) %>%
      dplyr::select(-bridge_intensity)
  }
  else {
    long_protein_dat <- long_protein_dat %>%
      mutate(Intensity = log2(Intensity + 1))
  }
  
  ## Subtract trimmed mean - represents a harsh batch correction
  if (scale_trim_mean) {
    trimmed_mean_dat <- long_protein_dat %>%
      group_by(protein.id, Batch) %>%
      summarize(mean_intensity = trimmed_mean(Intensity, trim = trim))
    long_protein_dat <- long_protein_dat %>%
      left_join(trimmed_mean_dat) %>%
      mutate(Intensity = Intensity - mean_intensity) %>%
      dplyr::select(-mean_intensity)
  }
  long_protein_dat
}

## Function to put intensities within a batch
## on the same scale without rolling up to proteins
peptide_batch_scaling <- function(long_peptide_dat,
                                  bridge_sample = NULL,
                                  na_as_zero = TRUE,
                                  zero_as_na = FALSE,
                                  batch_scaling = TRUE) {
  
  ## Derive normalization factors within batch
  ## Ratio of max cumulative intensities for a tag to each cumulative tag intensity 
  if (batch_scaling) {
    tmt_max <- long_peptide_dat %>%
      group_by(Batch, tag) %>%
      summarise(sum_inten = sum(Intensity, na.rm = TRUE)) %>%
      ungroup %>%
      group_by(Batch) %>%
      summarise(max_sum_inten = max(sum_inten))
    tmt_norm <- long_peptide_dat %>%
      group_by(sample.id) %>%
      summarise(sum_inten = sum(Intensity, na.rm = TRUE)) %>%
      ungroup %>%
      mutate(Batch = gsub(x = sample.id, pattern = "~.*$", replacement = "", perl = TRUE)) %>%
      left_join(tmt_max) %>%
      mutate(norm_factor = max_sum_inten/sum_inten) %>%
      dplyr::select(-c(sum_inten, max_sum_inten, Batch))
  }

  ## Sum peptides to quantify a protein
  long_peptide_dat <- long_peptide_dat %>%
    # Remove to add NAs
    complete(sample.id, nesting(peptide.id)) %>%
    # Add back protein ensembl id and gene symbol
    dplyr::select(peptide.id, sample.id, Intensity) %>%
    left_join(long_peptide_dat %>%
                dplyr::select(protein.id, peptide.id, symbol, balance, polymorphic) %>%
                distinct) %>% 
    # Add back Batch
    mutate(Batch = gsub(x = sample.id, pattern = "~.*$", replacement = "", perl = TRUE),
           tag = gsub(x = sample.id, pattern = "^[A-Z0-9]*~", replacement = "", perl = TRUE))
  
  ## Convert NAs to 0s
  if (na_as_zero) {
    long_peptide_dat <- long_peptide_dat %>%
      mutate(Intensity = ifelse(is.na(Intensity), 0, Intensity))
  }
  ## Turn 0s into NAs
  if (zero_as_na) {
    long_peptide_dat <- long_peptide_dat %>%
      mutate(Intensity = ifelse(Intensity == 0, NA, Intensity))
  }
  
  if (batch_scaling) {
    ## Merge in tag-specific normalizations
    long_peptide_dat <- long_peptide_dat %>%
      left_join(tmt_norm) %>%
      # Standardize across samples in a Batch
      mutate(Intensity = Intensity * norm_factor) %>%
      # Remove norm_factor
      dplyr::select(-norm_factor) %>%
      arrange(protein.id, peptide.id, sample.id)
  }
      
  ## Scaling by bridge sample
  if (!is.null(bridge_sample)) {
    bridge_dat <- long_peptide_dat %>%
      filter(tag == bridge_sample) %>%
      #mutate(Intensity = ifelse(Intensity == 0, 0.01, Intensity)) %>%
      dplyr::select(-c(sample.id, tag))
    long_peptide_dat <- long_peptide_dat %>%
      filter(tag != bridge_sample) %>%
      # Merging in bridge intensity
      left_join(bridge_dat %>%
                  dplyr::rename(bridge_intensity = Intensity)) %>%
      #mutate(Intensity = ifelse(Intensity == 0, 0.01, Intensity)) %>%
      mutate(Intensity = log2((Intensity + 1)/(bridge_intensity + 1))) %>%
      dplyr::select(-bridge_intensity)
  }
  else {
    long_peptide_dat <- long_peptide_dat %>%
      mutate(Intensity = log2(Intensity + 1))
  }
  long_peptide_dat
}

## Function to perform put intensities within a batch
## on the same scale
take_protein_sample_means <- function(long_peptide_dat) {
  
  ## Derive normalization factors within batch
  ## Ratio of max cumulative intensities for a tag to each cumulative tag intensity 
  tmt_max <- long_peptide_dat %>%
    group_by(Batch, tag) %>%
    summarise(sum_inten = sum(Intensity, na.rm = TRUE)) %>%
    ungroup %>%
    group_by(Batch) %>%
    summarise(max_sum_inten = max(sum_inten))
  tmt_norm <- long_peptide_dat %>%
    group_by(sample.id) %>%
    summarise(sum_inten = sum(Intensity, na.rm = TRUE)) %>%
    ungroup %>%
    mutate(Batch = gsub(x = sample.id, pattern = "~.*$", replacement = "", perl = TRUE)) %>%
    left_join(tmt_max) %>%
    mutate(norm_factor = max_sum_inten/sum_inten) %>%
    dplyr::select(-c(sum_inten, max_sum_inten, Batch))

  long_peptide_dat <- long_peptide_dat %>%
    left_join(tmt_norm) %>%
    mutate(Intensity = Intensity * norm_factor) %>%
    dplyr::select(-norm_factor)

  ## Sum peptides to quantify a protein
  long_protein_dat_mean <- long_peptide_dat %>%
    ## Summing peptides for a protein
    group_by(protein.id, sample.id) %>%
    summarize(Intensity = mean(log2(Intensity + 0.01))) %>%
    ungroup
  long_protein_dat_se <- long_peptide_dat %>%
    ## Summing peptides for a protein
    group_by(protein.id, sample.id) %>%
    summarise_at("Intensity", ~sd(log2(. + 0.01))/sqrt(n())) %>%
    mutate(se = Intensity) %>%
    ungroup %>%
    dplyr::select(-Intensity)
    
    # Add back Batch
  long_protein_dat <- inner_join(long_protein_dat_mean, long_protein_dat_se) %>%
    mutate(Batch = gsub(x = sample.id, pattern = "~.*$", replacement = "", perl = TRUE),
           tag = gsub(x = sample.id, pattern = "^TMT[0-9]*~", replacement = "", perl = TRUE))
    
  long_protein_dat
}



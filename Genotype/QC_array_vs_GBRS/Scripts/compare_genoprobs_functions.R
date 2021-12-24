# Ronnegard & Valdar 2011
## Meant for the CC
calc_SIM_genoprobs <- function(genoprobs,
                               map_list,
                               use_subjects = NULL, 
                               p_null = 0.125, 
                               just_these_loci = NULL) {
  
  ## Convert map list to map data.frame
  map_df <- qtl2convert::map_list_to_df(map_list)
  
  ## Grab loci
  loci <- map_df$marker
  chr <- map_df$chr
  pos <- map_df$pos
  if (!is.null(just_these_loci)) { 
    loci <- just_these_loci 
    chr <- map_df[loci,]$chr
    pos <- map_df[loci,]$pos
  }
  
  ## Grab subjects
  subjects <- dimnames(genoprobs[[1]])[[1]]
  if(!is.null(use_subjects)){
    subjects <- use_subjects
  }
  
  SIM_mat <- matrix(NA, nrow = length(subjects), ncol = length(loci))
  for (i in 1:length(loci)) {
    probs <- genoprobs[[chr[i]]][subjects,,loci[i]]
    
    SIM_vec <- rowSums(get_SIM(p_post = probs, p_null = p_null))
    SIM_mat[,i] <- SIM_vec
  }
  rownames(SIM_mat) <- subjects
  colnames(SIM_mat) <- loci
  
  ## Results
  SIM_list <- list(SIM_mat = SIM_mat, 
                   pos = pos,
                   chr = chr)
  SIM_list
}

get_SIM <- function(p_post, 
                    p_null = NULL) {
  
  if (is.null(p_null)) { p_null <- 1/length(p_post) }
  p_post[p_post == 0] <- .Machine$double.eps
  SIM_ind <- p_post*log(p_post/p_null)
  SIM_ind/log(1/p_null)
}

plot_SIM <- function(SIM_object, chr = 1, main = "") {
  
  use_chr <- SIM_object$chr == chr
  pos <- SIM_object$pos[use_chr]
  SIM <- SIM_object$SIM_mat[,use_chr]
  o <- order(pos)
  plot(x = pos[o], y = SIM[1,o], 
       type = "l", col = "grey", 
       ylim = c(0, 1), ylab = "Selective information content (SIC)", xlab = paste("Chr", chr, "position (Mb)"), main = main, las = 1, frame.plot = FALSE)
  for (i in 2:nrow(SIM)) {
    lines(x = pos[o], y = SIM[i,o], col="grey")
  }
  lines(x = pos[o], y = colMeans(SIM)[o], col = "black", lwd = 3)
}


compare_dist_genoprobs <- function(genoprobs1, genoprobs2, 
                                   map_list,
                                   use_dosages = FALSE) {
  ## Convert map list to map data.frame
  map_df <- qtl2convert::map_list_to_df(map_list)
  
  ## Grab loci
  loci <- map_df$marker
  chr <- map_df$chr
  pos <- map_df$pos
  
  ## Grab overlap in subjects
  subjects <- intersect(dimnames(genoprobs1[[1]])[[1]], dimnames(genoprobs2[[1]])[[1]])
  
  dist_mat <- matrix(NA, nrow = length(subjects), ncol = nrow(map_df))
  for (i in 1:length(loci)) {
    if (loci[i] %in% dimnames(genoprobs1[[chr[i]]])[[3]] & loci[i] %in% dimnames(genoprobs2[[chr[i]]])[[3]]) {
      probs1 <- genoprobs1[[chr[i]]][subjects,,loci[i]]
      probs2 <- genoprobs2[[chr[i]]][subjects,,loci[i]]
      
      dist_mat[,i] <- sapply(1:nrow(probs1), function(j) sqrt(sum((probs1[j,] - probs2[j,])^2)))
    }
    else {
      dist_mat[,i] <- rep(NA, nrow(probs1))
    }
  }
  rownames(dist_mat) <- subjects
  colnames(dist_mat) <- loci
  
  ## Remove NAs
  remove_loci_index <- which(is.na(dist_mat[1,]))
  if (length(remove_loci_index) > 0) {
    dist_mat <- dist_mat[,-remove_loci_index]
    pos <- pos[-remove_loci_index]
    chr <- chr[-remove_loci_index]
  }
  
  ## Results
  dim_list <- list(dist_mat = dist_mat, 
                   pos = pos,
                   chr = chr)
  dim_list
}

compare_dist_genoprobs_sex <- function(genoprobs_minimuga, 
                                       map_list,
                                       use_dosages = FALSE) {
  
  ## Convert map list to map data.frame
  map_df <- qtl2convert::map_list_to_df(map_list)
  
  ## Grab loci
  loci <- map_df$marker
  chr <- map_df$chr
  pos <- map_df$pos
  
  ## Grab overlap in subjects
  subjects <- dimnames(genoprobs_minimuga[[1]])[[1]]
  #browser()
  
  strains <- unique(gsub(x = subjects, pattern = "_.$", replacement = "", perl = TRUE))
  ## Check for male
  strains <- strains[paste(strains, "M", sep = "_") %in% subjects]
  ## Check for male
  strains <- strains[paste(strains, "F", sep = "_") %in% subjects]
  
  dist_mat <- matrix(NA, nrow = length(strains), ncol = nrow(map_df))
  for (i in 1:length(loci)) {
    full_probs <- genoprobs_minimuga[[chr[i]]][,,loci[i]]
    
    male_probs <- full_probs[paste(strains, "M", sep = "_"),]
    female_probs <- full_probs[paste(strains, "F", sep = "_"),]
    
    dist_mat[,i] <- sapply(1:length(strains), function(j) sqrt(sum((male_probs[j,] - female_probs[j,])^2)))
  }
  rownames(dist_mat) <- strains
  colnames(dist_mat) <- loci
  
  ## Remove NAs
  remove_loci_index <- which(is.na(dist_mat[1,]))
  if (length(remove_loci_index) > 0) {
    dist_mat <- dist_mat[,-remove_loci_index]
    pos <- pos[-remove_loci_index]
    chr <- chr[-remove_loci_index]
  }
  ## Results
  dim_list <- list(dist_mat = dist_mat, 
                   pos = pos,
                   chr = chr)
  dim_list
}


plot_dist <- function(dist_object,
                      chr = 1,
                      main=""){
  
  use_chr <- dist_object$chr == chr
  pos <- dist_object$pos[use_chr]
  dist <- dist_object$dist_mat[,use_chr]
  max <- which.max(rowMeans(dist_object$dist_mat))
  min <- which.min(rowMeans(dist_object$dist_mat))
  
  o <- order(pos)
  plot(x = pos[o], y = dist[1, o], type = "l", pch=20, col="grey", 
       ylim = c(0, max(dist_object$dist_mat)), 
       ylab = "Euclidian distance", xlab = paste("Chr", chr, "position (Mb)"), 
       frame.plot = FALSE,
       main = main, 
       las = 1)
  for (i in 2:nrow(dist)){
    lines(x = pos[o], y = dist[i,o], type="l", pch=20, col="grey")
  }
  lines(x = pos[o], y = colMeans(dist)[o], col = "black", lwd = 3)
  lines(x = pos[o], y = dist[max,o], type="l", pch=20, col="red")
  lines(x = pos[o], y = dist[min,o], type="l", pch=20, col="blue")
  abline(h = dist(rbind(c(1,0,0,0,0,0,0,0), c(0,1,0,0,0,0,0,0))), lty = 2, col = "cornflowerblue")
}

# compare_genoprobs_from <- function(genoprobs1, mouse_id, genoprobs2, map_list) {
#   
#   ## Convert map list to map data.frame
#   map_df <- qtl2convert::map_list_to_df(map_list)
#   
#   ## Grab loci
#   loci <- map_df$marker
#   chr <- map_df$chr
#   pos <- map_df$pos
#   
#   ## Grab overlap in subjects
#   subjects <- dimnames(genoprobs2[[1]])[[1]]
#   
#   dist_mat <- matrix(NA, nrow = length(subjects), ncol = nrow(map_df))
#   for (i in 1:length(loci)) {
#     if (loci[i] %in% dimnames(genoprobs1[[chr[i]]])[[3]] & loci[i] %in% dimnames(genoprobs2[[chr[i]]])[[3]]) {
#       probs1 <- genoprobs1[[chr[i]]][mouse_id,,loci[i]]
#       probs2 <- genoprobs2[[chr[i]]][,,loci[i]]
#       
#       dist_mat[,i] <- sapply(1:nrow(probs2), function(j) sqrt(sum((probs1 - probs2[j,])^2)))
#     }
#     else {
#       dist_mat[,i] <- rep(NA, nrow(probs1))
#     }
#   }
#   rownames(dist_mat) <- subjects
#   colnames(dist_mat) <- loci
#   
#   ## Remove NAs
#   remove_loci_index <- which(is.na(dist_mat[1,]))
#   if (length(remove_loci_index) > 0) {
#     dist_mat <- dist_mat[,-remove_loci_index]
#     pos <- pos[-remove_loci_index]
#     chr <- chr[-remove_loci_index]
#   }
#   
#   ## Results
#   dim_list <- list(dist_mat = dist_mat, 
#                    pos = pos,
#                    chr = chr)
#   dim_list
#   
# }

check_duplicate_genoprobs <- function(genoprobs, mouse_id, map_list) {
  
  ## Convert map list to map data.frame
  map_df <- qtl2convert::map_list_to_df(map_list)
  
  ## Grab loci
  loci <- map_df$marker
  chr <- map_df$chr
  pos <- map_df$pos
  
  ## Grab overlap in subjects
  subjects <- dimnames(genoprobs[[1]])[[1]]
  
  dist_mat <- matrix(NA, nrow = length(subjects), ncol = nrow(map_df))
  for (i in 1:length(loci)) {
    if (loci[i] %in% dimnames(genoprobs[[chr[i]]])[[3]]) {
      probs1 <- genoprobs[[chr[i]]][rep(mouse_id, length(subjects)),,loci[i]]
      probs2 <- genoprobs[[chr[i]]][,,loci[i]]
      
      dist_mat[,i] <- sapply(1:nrow(probs2), function(j) sqrt(sum((probs1[j,] - probs2[j,])^2)))
    }
    else {
      dist_mat[,i] <- rep(NA, nrow(probs1))
    }
  }
  rownames(dist_mat) <- subjects
  colnames(dist_mat) <- loci
  
  ## Remove NAs
  remove_loci_index <- which(is.na(dist_mat[1,]))
  if (length(remove_loci_index) > 0) {
    dist_mat <- dist_mat[,-remove_loci_index]
    pos <- pos[-remove_loci_index]
    chr <- chr[-remove_loci_index]
  }
  
  ## Results
  dim_list <- list(dist_mat = dist_mat, 
                   pos = pos,
                   chr = chr)
  dim_list
}

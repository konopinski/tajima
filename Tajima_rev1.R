library(parallel)
library(Biostrings)
library(vcfR)
library(dartR)
library(wesanderson)
library(ggplot2)

#########################################

folder <- "/mnt/Dane/Szop/RADy/do_plotor/RAD/" 
resFolder <- "/mnt/Dane/Szop/RADy/do_plotor/RAD/tajima/"
popmap <- paste0(folder,"pop_mapy/popmap_RAD_B_S.txt")
popsFile <- read.table(popmap, header = FALSE, sep = "\t", col.names = c("ID", "pop"))

vcf_path <- paste0(folder,"vcf_maf001.recode.vcf")

#########################################
get_strict_categories <- function(vcf_path) {
  
  vcf <- read.vcfR(vcf_path, verbose = FALSE)
  info_data <- vcf@fix[, "INFO"]
  
  # extracting ANN=...
  ann_list <- sub(".*ANN=([^;]+).*", "\\1", info_data)
  
  # List of strong effects (blocking "Strict 4-fold" status)
  strong_effects <- c(
    "missense_variant", "stop_gained", "stop_lost", "start_lost",
    "splice_region_variant&intron_variant", 
    "splice_region_variant&synonymous_variant",
    "missense_variant&splice_region_variant", 
    "splice_donor_variant&intron_variant",
    "splice_acceptor_variant&intron_variant",
    "splice_region_variant&stop_retained_variant", 
    "stop_gained&splice_region_variant",
    "stop_lost&splice_region_variant", 
    "inframe_insertion", "inframe_deletion", "frameshift_variant"
  )
  
  # all coding effects 
  coding_effects <- c(strong_effects, "synonymous_variant", "stop_retained_variant", "start_retained_variant")
  
  # 4-fold degenerated aminoacids
  four_fold_aas <- c("p.Ala", "p.Gly", "p.Pro", "p.Thr", "p.Val", "p.Arg", "p.Leu", "p.Ser")
  
  # Helper function - codon position check
  is_third_pos <- function(c_string) {
    if (is.na(c_string) || c_string == "") return(FALSE)
    # 1. must start with 'c.' (coding)
    if (!grepl("^c\\.", c_string)) return(FALSE)
    # 2. rejecting UTRs ("*", "-) and Introns ("+")
    if (grepl("[\\*\\-\\+]", c_string)) return(FALSE)
    # 3. esctraction of position in CDS
    num_match <- regmatches(c_string, regexpr("[0-9]+", c_string))
    if (length(num_match) == 0) return(FALSE)
    # 4. checking if the SNP is located at 3rd codon position
    return(as.numeric(num_match) %% 3 == 0)
  }
  
  # processing snps
  categories <- t(sapply(strsplit(ann_list, ","), function(transcripts) {
    details <- strsplit(transcripts, "\\|")
    
    # extracting data with safety length check
    effects <- sapply(details, "[", 2)
    hgvs_c  <- sapply(details, function(x) if(length(x) >= 10) x[10] else "")
    hgvs_p  <- sapply(details, function(x) if(length(x) >= 11) x[11] else "")
    
    # 1. Is it coding snp?
    is_coding_snp <- any(effects %in% coding_effects)
    
    # 2. Is it non-coding?
    is_noncoding_snp <- !is_coding_snp
    
    # 3. Is it strict 4-fold:
    is_4fold_transcript <- (effects == "synonymous_variant") & 
      grepl(paste(four_fold_aas, collapse="|"), hgvs_p) &
      sapply(hgvs_c, is_third_pos)
    
    has_4fold_hit  <- any(is_4fold_transcript)
    has_deleterious <- any(effects %in% strong_effects)
    
    is_strict_4fold <- has_4fold_hit && !has_deleterious
    
    return(c(coding = is_coding_snp, 
             noncoding = is_noncoding_snp, 
             fourfold = is_strict_4fold))
  }))
  
  rm(vcf)
  gc()
  
  return(as.data.frame(categories))
}


######### FUNCTION: Prepare Genlight (Fixes Compliance & Metadata) #########
prepare_genlight <- function(vcf_path, pops_df=popsFile) {
  vcf <- read.vcfR(vcf_path)
  gl  <- vcfR2genlight(vcf)
  
  # Initial compliance check to initialize slots
  gl <- gl.compliance.check(gl)
  
  # Manual metadata injection
  gl$loc.all <- apply(vcf@fix[,4:5], 1, function(x) paste0(x, collapse = "/"))
  gl@ploidy  <- rep(2L, length(gl@ind.names))
  
  # Secure population mapping
  if (pops_df$ID[1] == "ind") pops_df <- pops_df[-1,] 
  pop(gl) <- pops_df$pop[match(indNames(gl), pops_df$ID)]
  
  # Filter monomorphs and final compliance check
  gl <- gl.filter.monomorphs(gl)
  gl <- gl.compliance.check(gl)
  
  rm(vcf); gc()
  return(gl)
}

######### FUNCTION: Calculate SNP Metrics #########
calculate_snp_metrics <- function(gl_obj, label,pops) {
  sapply(pops, simplify = FALSE, USE.NAMES = TRUE, FUN = function(p) {
    # Keep population and recalculate loc.metrics (maf) for compliance
    popGl <- gl.keep.pop(gl_obj, p, recalc = TRUE)
    
    segregating <- popGl@other$loc.metrics$maf != 0
    gl_mat      <- as.matrix(popGl)
    n           <- 2L * as.integer(colSums(!is.na(gl_mat)))
    
    alt_counts  <- as.integer(colSums(gl_mat, na.rm=TRUE))
    pi_snp      <- ifelse(n > 1, (n^2 - (alt_counts^2 + (n - alt_counts)^2)) / (n * (n - 1)), NA)
    
    a1_snp <- ifelse(n > 1, sapply(n, function(i) sum(1 / seq_len(i - 1))), 0)
    a2_snp <- ifelse(n > 1, sapply(n, function(i) sum(1 / seq_len(i - 1)^2)), 0)
    b1_snp <- ifelse(n > 1, (n + 1) / (3 * (n - 1)), 0)
    b2_snp <- ifelse(n > 1, 2 * (n^2 + n + 3) / (9 * n * (n - 1)), 0)
    
    c1_snp <- ifelse(n > 1, b1_snp - (1 / a1_snp), 0)
    c2_snp <- ifelse(n > 1, b2_snp - ((n + 2) / (a1_snp * n)) + (a2_snp / a1_snp^2), 0)
    e1_snp <- ifelse(n > 1, c1_snp / a1_snp, 0)
    e2_snp <- ifelse(n > 1, c2_snp / (a1_snp^2 + a2_snp), 0)
    
    thetaW_snp <- ifelse(n > 1, as.integer(segregating) * (1/a1_snp), 0)
    tajD_snp   <- ifelse(sqrt(e1_snp) > 0, (pi_snp - thetaW_snp) / sqrt(e1_snp), NA)
    
    df_out <- data.frame(
      "CHROM" = popGl@chromosome, "SNP" = popGl@position, "locus" = locNames(popGl), 
      "pi" = round(pi_snp, 5), "thetaW" = round(thetaW_snp, 5), "D" = round(tajD_snp, 5)
    )
    # write.table(df_out, file = paste0(resFolder, p, "_", label, "_SNP_metrics.xls"), 
    #             sep = "\t", quote = FALSE, row.names = FALSE)
    
    return(df_out)
  })
}

############# FUNCTION preparing files with all pops indexes ##############


indexWise <- function(PopEstim, pops, functions){
  info <- PopEstim[[1]][,1:3]
  pi <- cbind(info,sapply(PopEstim, function(x)x[,4]),functions)
  theta <- cbind(info,sapply(PopEstim, function(x)x[,5]),functions)
  D <- cbind(info,sapply(PopEstim, function(x)x[,6]),functions)
  write.table(pi, file = paste0(resFolder, "pi.xls"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(theta, file = paste0(resFolder, "theta.xls"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(D, file = paste0(resFolder, "D.xls"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  return(list(pi,D,theta))
}

######### EXECUTION #########
# For vcf_maf001.recode.vcf
gl_standard  <- prepare_genlight(vcf_path, popsFile)
functions <- get_strict_categories(vcf_path)
PopEstim     <- calculate_snp_metrics(gl_standard, "standard", pops = levels(pop(gl_standard)))
piDtheta <- indexWise(PopEstim,pops,functions)



############# 4-fold codons  ##############
fourFolds <- c("p.Val", "p.Thr", "p.Ala", "p.Gly", "p.Pro", "p.Arg", "p.Leu", "p.Ser")
# 

# 1. Check if it is synonymous variant in coding region
is_cds_synonymous <- grepl("synonymous_variant", ann_field) &
  grepl("c\\.[0-9]+[A-Z]>[A-Z]", ann_field)

# 2. Extract positions only if condition 1 is fulfilled 
sapply(is_cds_synonymous,function(x)
if(x) {
  # Extract position number 
  pos_match <- regmatches(ann_field, regexpr("c\\.([0-9]+)", ann_info))
  pos_num <- as.numeric(gsub("c\\.", "", pos_match))

  # Check if it 3rd position
  is_4fold <- (pos_num %% 3 == 0)
} else {
  is_4fold <- FALSE
})


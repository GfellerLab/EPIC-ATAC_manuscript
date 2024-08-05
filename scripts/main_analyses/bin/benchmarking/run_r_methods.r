# Libraries ---------------------------------------------------------------
library(Metrics)
library(EPICATAC)
library(GenomicRanges)
library(readxl)

load("PBMC_ATAC_Ref_Major.rda")
load("Tumor_ATAC_Ref_Major.rda")
load("PBMC_ATAC_Ref_Subtypes.rda")
load("Tumor_ATAC_Ref_Subtypes.rda")

# Parameters --------------------------------------------------------------
bulk_path = as.character(commandArgs(TRUE)[1])
bulk_data_name = as.character(commandArgs(TRUE)[2])
with_subtypes = as.logical(commandArgs(TRUE)[3])
deconv_method = as.character(commandArgs(TRUE)[4])


# Load markers and profiles  -------------------------------------------------
if(grepl(pattern = 'PBMC', bulk_data_name)){
  ref_major = PBMC_ATAC_Ref_Major
  if(with_subtypes){
    ref_withSubtypes = PBMC_ATAC_Ref_Subtypes
  } 
}else{
  ref_major = Tumor_ATAC_Ref_Major
  if(with_subtypes){
    ref_withSubtypes = Tumor_ATAC_Ref_Subtypes
  } 
}  

if(deconv_method == "EPIC"){
  # Load bulk  ---------------------------------------------------------------
  load(paste0(bulk_data_name, "_bulk.Rdata"))
  bulks_data_raw = bulks_data$bulk[, colnames(bulks_data$obs)[which(colnames(bulks_data$obs) != "cell_type")], drop = F]
  
}else{
  matched_bulks_data = as.data.frame(read_excel(bulk_path, sheet = 1, col_names = T))
  rownames(matched_bulks_data) = matched_bulks_data[, 1]
  matched_bulks_data = matched_bulks_data[, 2:ncol(matched_bulks_data), drop = F]
  
}
  




# Run EPIC, QuantiSeq, MPC_counter and ABIS -------------------------------
if(deconv_method == "EPIC"){
  
  if(bulk_data_name %in% c("Gynecological_cancers","PBMC_multiome", "PBMC_multiome_simulated", "PBMC_experiment", "HNSCC","BRCA","CRC","CESC","SKCM","PDAC", "OV", "CEAD", "GBM", "UCEC")){
    genome_version = "hg38" 
  }else{
    genome_version = "hg19" 
  }   
  
  ref = get(ifelse(with_subtypes, "ref_withSubtypes", "ref_major"))
  epic_new_res = EPIC_ATAC(
    bulk = as.matrix(bulks_data_raw),
    reference = ref,
    withOtherCells = ifelse(grepl("PBMC", bulk_data_name), F, T),
    ATAC = T, genome_version = genome_version
  )
  
  epic_output = as.data.frame(epic_new_res$cellFractions)
  epic_output$samples = rownames(epic_output)
  write.table(epic_output, file = paste0(bulk_data_name, "_EPIC_predictions.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
 
}

if(deconv_method == "quantiSeq"){
  ref = get(ifelse(with_subtypes, "ref_withSubtypes", "ref_major"))
  quantiseq_output <- quantiseqr:::quanTIseq(currsig = ref$refProfiles[ref$sigPeaks,],
                                             currmix = matched_bulks_data[ref$sigPeaks, , drop=F] ,
                                             scaling = rep(1, ncol(ref$refProfiles)),
                                             method = "lsei")
  quantiseq_output <- as.data.frame(quantiseq_output)
  quantiseq_output$samples <- rownames(quantiseq_output)
  write.table(quantiseq_output, file = paste0(bulk_data_name, "_quantiseq_predictions.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
  
}

if(deconv_method == "MPC_counter"){
  ref = get(ifelse(with_subtypes, "ref_withSubtypes", "ref_major"))
  MCPcounter_res <- MCPcounter:::appendSignatures(matched_bulks_data, ref$markers)
  MCPcounter_output <- MCPcounter_res
  MCPcounter_output$samples <- rownames(MCPcounter_output)
  write.table(MCPcounter_output, file = paste0(bulk_data_name, "_MCPcounter_predictions.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
  
}

if(deconv_method == "ABIS"){
  ref = get(ifelse(with_subtypes, "ref_withSubtypes", "ref_major"))
  features <- intersect(rownames(matched_bulks_data), ref$sigPeaks)
  Dec <- (apply(matched_bulks_data[features, , drop = F] , 2, function(x) stats::coef(MASS::rlm(as.matrix(ref$refProfiles[features, ]), x, maxit = 100)))) * 100
  ABIS_res <- t(signif(Dec, 3))
  ABIS_output <- (ABIS_res - apply(ABIS_res, 1, min, na.rm = TRUE)) / (apply(ABIS_res, 1, max, na.rm = TRUE) - apply(ABIS_res, 1, min, na.rm = TRUE))
  ABIS_output <- as.data.frame(ABIS_output / rowSums(ABIS_output))
  ABIS_output$samples <- rownames(ABIS_output)
  write.table(ABIS_output, file = paste0(bulk_data_name, "_ABIS_predictions.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
  
}
  

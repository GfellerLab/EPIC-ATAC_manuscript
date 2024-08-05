# Libraries ---------------------------------------------------------------
library(GenomicRanges)
library(openxlsx)
library(liftOver)
library(EPICATAC)

# Functions ---------------------------------------------------------------
load("PBMC_ATAC_Ref_Major.rda")
load("Tumor_ATAC_Ref_Major.rda")
load("PBMC_ATAC_Ref_Subtypes.rda")
load("Tumor_ATAC_Ref_Subtypes.rda")

# Parameters --------------------------------------------------------------
bulk_path=as.character(commandArgs(TRUE)[1])  
bulk_data_name = as.character(commandArgs(TRUE)[2])
deconPeaker_sign_path = as.character(commandArgs(TRUE)[3])
with_subtypes = as.logical(commandArgs(TRUE)[4])

if(bulk_data_name %in% c("Gynecological_cancers","PBMC_multiome", "PBMC_experiment", "PBMC_multiome_simulated", "HNSCC","BRCA","CRC","CESC","SKCM","PDAC", "OV", "CEAD", "GBM", "UCEC")){
  genome_version = "hg38" 
}else{
  genome_version = "hg19" 
}   
  
# Load bulk  ---------------------------------------------------------------
load(bulk_path)
bulks_data$bulk = bulks_data$bulk[, colnames(bulks_data$obs)[which(colnames(bulks_data$obs) != "cell_type")], drop = F]

# Load markers and profiles  -------------------------------------------------

if(grepl(pattern = 'PBMC',bulk_data_name)){
  if(with_subtypes){
    ref = PBMC_ATAC_Ref_Subtypes
  } else{
    ref = PBMC_ATAC_Ref_Major
  }  
}else{
  if(with_subtypes){
    ref = Tumor_ATAC_Ref_Subtypes
  } else{
    ref = Tumor_ATAC_Ref_Major
  }  
}  

# Match peaks from bulk and our reference -------------------------------------
if (genome_version == "hg19") {
  lifted_matrix = EPICATAC:::run_liftOver(bulk_matrix = bulks_data$bulk, from = genome_version)
} else {
  lifted_matrix = bulks_data$bulk
} 
matched_bulk = EPICATAC:::match_peaks(bulk_matrix = lifted_matrix, profile_features = rownames(ref$refProfiles)) 


# # Save EPIC input files ---------------------------------------------------
write.xlsx(list(bulk_data=matched_bulk,
                bulk_prop=bulks_data$obs),
           file = paste0(bulk_data_name,"_epic_bulk_input.xlsx"),rowNames=T)

# Save bulk true proportions ----------------------------------------------
bulk_proportions = bulks_data$obs
bulk_proportions$cell_types = rownames(bulk_proportions)
write.table(bulk_proportions,file=paste0(bulk_data_name,"_bulk_proportions.txt"),quote = F,row.names = F,col.names = T, sep = "\t")

# Save CIBERSORT and deconPeaker mixture --------------------------------------------------
write.table(cbind(data.frame(sample = rownames(matched_bulk)),matched_bulk),
            file = paste0(bulk_data_name,"_cibersort_mixture.txt"),row.names = F,col.names = T,sep="\t",quote = F)

peaks_ids <- as.data.frame(do.call(rbind, strsplit(rownames(matched_bulk), "[:-]")))
colnames(peaks_ids) <- c("chr","start","stop")
write.table(cbind(peaks_ids, matched_bulk),
            file = paste0(bulk_data_name,"_deconPeaker_bulk_input.txt"),row.names = F,col.names = T,sep="\t",quote = F)


# Save our reference with our markers as deconpeaker and cibersort --------

# Deconpeaker ref input 
peaks_ids <- as.data.frame(do.call(rbind, strsplit(ref$sigPeaks, "[:-]")))
colnames(peaks_ids) <- c("chr","start","stop")
write.table(cbind(peaks_ids, ref$refProfiles[ref$sigPeaks,]),
            file = paste0(bulk_data_name,"_deconPeaker_ref_input.txt"),row.names = F,col.names = T,sep="\t",quote = F)

# CIBERSORTX ref input 
write.table(cbind(data.frame(NAME = ref$sigPeaks), ref$refProfiles[ref$sigPeaks,]),
            file = paste0(bulk_data_name,"_CIBERSORTx_ref_input.txt"),row.names = F,col.names = T,sep="\t",quote = F)

# Match peaks to the DeconPeaker reference signature matrix ---------------
deconPeaker_signature = read.table(deconPeaker_sign_path, header = T)
deconPeaker_peaks = sapply(1:nrow(deconPeaker_signature),function(i) paste(c(deconPeaker_signature$chrom[i],
                                                                           deconPeaker_signature$start[i],
                                                                           deconPeaker_signature$end[i]), collapse = "-"))

# Deconpeaker original markers are in hg19, so we need to liftover peaks coordinates of datasets processed with hg38 
run_liftOver <- function(bulk_matrix, from = "hg19"){
  
  peaks_gr <- EPICATAC:::StringToGRanges(rownames(bulk_matrix), sep = c("[:-]","-")) 
  peaks_gr$id <- rownames(bulk_matrix)
  
  # Lift over coordinates 
  if(from == "hg19"){
    path <- system.file(package = "liftOver", "extdata", "hg19ToHg38.over.chain")
  }else if(from == "hg38"){
    path <- system.file(package = "liftOver", "extdata", "hg38ToHg19.over.chain")
  } 
  ch <- rtracklayer::import.chain(path)
  
  GenomeInfoDb::seqlevelsStyle(ch) <- "UCSC"
  regions_gr_converted <- unlist(rtracklayer::liftOver(peaks_gr, ch))
  
  new_regions <- data.frame(
    chr = regions_gr_converted@seqnames,
    start = rtracklayer::start(regions_gr_converted),
    end = rtracklayer::end(regions_gr_converted),
    id = regions_gr_converted$id
  )
  
  # remove duplicated regions 
  if(sum(duplicated(new_regions$id)) != 0){
    new_regions <- new_regions[-which(duplicated(new_regions$id)), ]
  } 
  
  new_regions$new_names <- paste0(new_regions$chr, "-", new_regions$start, "-", new_regions$end)
  rownames(new_regions) <- new_regions$id
  
  # remove duplicated regions 
  if(sum(duplicated(new_regions$new_names)) != 0){
    new_regions <- new_regions[-which(duplicated(new_regions$new_names)), ]
  } 
  
  lifted_matrix <- bulk_matrix[which(rownames(bulk_matrix) %in% new_regions$id),,drop = F]
  rownames(lifted_matrix) <- new_regions[rownames(lifted_matrix), "new_names"]
  
  return(lifted_matrix)
} 
if (genome_version == "hg19") {
  lifted_matrix = bulks_data$bulk
} else {
  lifted_matrix = run_liftOver(bulk_matrix = bulks_data$bulk, from = "hg38")
}
matched_bulk <- EPICATAC:::match_peaks(bulk_matrix = lifted_matrix, profile_features = deconPeaker_peaks) 

peaks_ids <- data.frame(chrom=sapply(rownames(matched_bulk),function(i) unlist(strsplit(i,"-"))[1]),
                        start=sapply(rownames(matched_bulk),function(i) unlist(strsplit(i,"-"))[2]),
                        stop=sapply(rownames(matched_bulk),function(i) unlist(strsplit(i,"-"))[3]) )
write.table(cbind(peaks_ids,matched_bulk),
            file = paste0(bulk_data_name,"_OM_deconPeaker_bulk_input.txt"),row.names = F,col.names = T,sep="\t",quote = F)


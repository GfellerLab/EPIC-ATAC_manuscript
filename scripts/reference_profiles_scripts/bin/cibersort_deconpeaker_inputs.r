library(openxlsx)
library(GenomicRanges)
library(rtracklayer)

# Parameters --------------------------------------------------------------
counts_path = as.character(commandArgs(TRUE)[1])
metadata_path = as.character(commandArgs(TRUE)[2]) 
major_groups = as.logical(commandArgs(TRUE)[3])

# Load counts and metadata --------------------------------------------------
counts = read.table(counts_path, header = T, sep = "\t", row.names = 1)
metadata = read.table(metadata_path, header = T, sep = "\t")

if(major_groups){
  metadata$groups = plyr::revalue(metadata$groups,replace = c("Naive_CD4_Tcells"="CD4_Tcells","Non_Naive_CD4_Tcells"="CD4_Tcells", "Tregs"="CD4_Tcells" ,
                                                               "Naive_CD8_Tcells"="CD8_Tcells" ,"Non_Naive_CD8_Tcells"="CD8_Tcells" ))
}else{
  if(length(which(metadata$groups %in% c("CD4_Tcells","CD8_Tcells"))!=0)){
    counts = counts[, -which(metadata$groups %in% c("CD4_Tcells","CD8_Tcells"))] 
    metadata = metadata[-which(metadata$groups %in% c("CD4_Tcells","CD8_Tcells")), ] 
  } 

}   

peaks_gr = Signac::StringToGRanges(rownames(counts),sep = c(":","-"))
region.length = peaks_gr@ranges@width
names(region.length) = rownames(counts)
counts <-  counts / region.length[rownames(counts)]
tpm.mat <- t(t(counts) * 1e6 / colSums(counts))

# Save Cibersortx input files to generate profiles ------------------------
summary(colnames(tpm.mat) == metadata$sample)
tpm_subset = tpm.mat[, which(! metadata$groups %in% c("Monocytes"))] 
metadata_subset = metadata[which(! metadata$groups %in% c("Monocytes")), ] 

cibersort_sorted_samples = tpm_subset
colnames(cibersort_sorted_samples) = metadata_subset$groups
genes = gsub(":","-",rownames(tpm_subset))
cibersort_sorted_samples = cbind(genes, cibersort_sorted_samples)
write.table(cibersort_sorted_samples, file = "cibersortx_input_ref.txt", sep = "\t", row.names = F, col.names = F, quote = F)

pheno_data = vector()
for(cell_type in unique(metadata_subset$groups)){
  pheno_status = rep(2, nrow(metadata_subset))
  pheno_status[which(metadata_subset$groups == cell_type)] = 1
  pheno_data = rbind(pheno_data, pheno_status)
} 
colnames(pheno_data) = metadata_subset$sample
pheno_data = cbind(data.frame(cell_type = unique(metadata_subset$groups)), pheno_data)
write.table(pheno_data, file = "cibersortx_input_pheno.txt", sep = "\t", row.names = F, col.names = F, quote = F)


# Save Deconpeaker input files to generate profiles -----------------------
peaks_gr = Signac::StringToGRanges(rownames(counts), sep = c(":", "-"))
peaks_coordinates = data.frame(chrom = rtracklayer::chrom(peaks_gr), start = start(peaks_gr), end = end(peaks_gr))
deconpeaker_sorted_samples = cbind(peaks_coordinates, tpm_subset)

write.table(deconpeaker_sorted_samples, file = "deconpeaker_input_ref.txt", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(pheno_data, file = "deconpeaker_input_pheno.txt", sep = "\t", row.names = F, col.names = F, quote = F)



# Save a PBMC reference: without the fibroblasts, the endothelial and macrophages 
tpm_subset = tpm.mat[, which(! metadata$groups %in% c("Macrophages", "Fibroblasts", "Endothelial"))] 
cibersort_sorted_samples = tpm_subset
metadata_subset = metadata[which(! metadata$groups %in% c("Macrophages", "Fibroblasts", "Endothelial")), ] 

colnames(cibersort_sorted_samples) = metadata_subset$groups
genes = gsub(":","-",rownames(tpm_subset))
cibersort_sorted_samples = cbind(genes, cibersort_sorted_samples)
write.table(cibersort_sorted_samples, file = "cibersortx_input_PBMCref.txt", sep = "\t", row.names = F, col.names = F, quote = F)

pheno_data = vector()
for(cell_type in unique(metadata_subset$groups)){
  pheno_status = rep(2, nrow(metadata_subset))
  pheno_status[which(metadata_subset$groups == cell_type)] = 1
  pheno_data = rbind(pheno_data, pheno_status)
} 
colnames(pheno_data) = metadata_subset$sample
pheno_data = cbind(data.frame(cell_type = unique(metadata_subset$groups)), pheno_data)
write.table(pheno_data, file = "cibersortx_input_PBMCpheno.txt", sep = "\t", row.names = F, col.names = F, quote = F)


peaks_gr = Signac::StringToGRanges(rownames(counts), sep = c(":", "-"))
peaks_coordinates = data.frame(chrom = rtracklayer::chrom(peaks_gr), start = start(peaks_gr), end = end(peaks_gr))
deconpeaker_sorted_samples = cbind(peaks_coordinates, tpm_subset)

write.table(deconpeaker_sorted_samples, file = "deconpeaker_input_PBMCref.txt", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(pheno_data, file = "deconpeaker_input_PBMCpheno.txt", sep = "\t", row.names = F, col.names = F, quote = F)

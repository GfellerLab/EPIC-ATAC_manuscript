library(openxlsx)
library(GenomicRanges)
library(rtracklayer)

# Parameters --------------------------------------------------------------
counts_path = as.character(commandArgs(TRUE)[1])
metadata_path = as.character(commandArgs(TRUE)[2]) 
markers_path = as.character(commandArgs(TRUE)[3])
major_groups = as.logical(commandArgs(TRUE)[4])


# Load counts and metadata --------------------------------------------------
counts = read.table(counts_path, header = T, sep = "\t", row.names = 1)
metadata = read.table(metadata_path, header = T, sep = "\t")
markers_df = read.table(markers_path, header = T, sep = "\t")

if(major_groups){
  metadata$groups = plyr::revalue(metadata$groups,replace = c("Naive_CD4_Tcells"="CD4_Tcells","Non_Naive_CD4_Tcells"="CD4_Tcells", "Tregs"="CD4_Tcells" ,
                                                                "Naive_CD8_Tcells"="CD8_Tcells" ,"Non_Naive_CD8_Tcells"="CD8_Tcells" ))
}else{
  if(length(which(metadata$groups %in% c("CD4_Tcells","CD8_Tcells"))!=0)){
    counts = counts[, -which(metadata$groups %in% c("CD4_Tcells","CD8_Tcells"))] 
    metadata = metadata[-which(metadata$groups %in% c("CD4_Tcells","CD8_Tcells")), ] 
  } 
 
}  


# Save all features profiles ----------------------------------------------
peaks_gr = Signac::StringToGRanges(rownames(counts),sep = c(":","-"))
region.length = peaks_gr@ranges@width
names(region.length) = rownames(counts)
counts <-  counts / region.length[rownames(counts)]

tpm.mat <- t(t(counts) * 1e6 / colSums(counts))
matrix_to_use = tpm.mat
rownames(matrix_to_use) = gsub(":","-",rownames(matrix_to_use))
markers_df$peak_id = gsub(":","-",markers_df$peak_id)

profiles = data.frame(row.names = rownames(matrix_to_use))
profiles.var = data.frame(row.names = rownames(matrix_to_use))

for(i in 1:length(unique(metadata$groups))) {
  c_type = as.vector(unique(metadata$groups))[i]
  print(c_type)
  cell_index = which(metadata$groups %in% c_type)
  profiles = cbind(profiles, matrixStats::rowMedians(matrix_to_use[, cell_index]))
  colnames(profiles)[ncol(profiles)] = c_type
  
  profiles.var = cbind(profiles.var, matrixStats::rowIQRs(matrix_to_use[, cell_index]))
  colnames(profiles.var)[ncol(profiles.var)] = c_type
}
profiles <- t(t(profiles) * 1e6 / colSums(profiles))
profiles_list = list(
  refProfiles = profiles,
  refProfiles.var = profiles.var
)


save(profiles_list, markers_df, file = paste0("profile.Rdata"))




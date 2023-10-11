
# Libraries ---------------------------------------------------------------
library(rtracklayer)


# Functions --------------------------------------------------------------
check_correlation <- function(markers_list, ATAC_ref, min_markers = 2, use_caret = T, cor_thr = 0.75){
  markers_list_updated = markers_list
  for(cell_type in names(markers_list) ){
    markers_subset = gsub(":","-",markers_list[[cell_type]])

    atac_mat_subset = log(ATAC_ref[which(rownames(ATAC_ref) %in% markers_subset),]+1)
    markers_subset = rownames(atac_mat_subset)
    cor_matrix_pval = psych::corr.test(t(atac_mat_subset),method = "pearson", ci = F, adjust = "none")$p  
    cor_matrix_estimate = cor(t(atac_mat_subset),method = "pearson")  
    cor_matrix_pval[which(cor_matrix_estimate<0)] = 1
    nb_cor = sapply(rownames(cor_matrix_pval), function(i) sum(cor_matrix_pval[i,colnames(cor_matrix_pval) != i] < 0.05/nrow(cor_matrix_pval)) )
    prop_cor = nb_cor/nrow(cor_matrix_pval)
    
    prop_min = quantile(prop_cor,cor_thr)
    print("after correlation")
    if(use_caret){
      markers_subset = markers_subset[caret::findCorrelation(cor_matrix_estimate, cutoff = quantile(cor_matrix_estimate, cor_thr))]
      
    } else{
      markers_subset = markers_subset[which(prop_cor>=prop_min)]
      
    } 
    
    if(length(markers_subset) >= min_markers){
      markers_list_updated[[cell_type]] = markers_subset
    } 
    
  }
  return(markers_list_updated)
} 

# Parameters --------------------------------------------------------------

profile_path = as.character(commandArgs(TRUE)[1])
TCGA_data = as.character(commandArgs(TRUE)[2])
human_atlas_peaks = as.character(commandArgs(TRUE)[3])


# Load reference profiles -------------------------------------------------
load(profile_path)

#  Human Atlas filtering --------------------------------------------------
print("Human atlas filtering")

overlapping_markers <- vector()
cell_types_modules <- list("immune" = c(8:25), "fibroblasts" = c(41:49, 139:150), "endothelial" = c(26:35))
cell_types_names <- list("immune" = c("Macrophages","CD8_Tcells","NK","Monocytes","CD4_Tcells","Neutrophils","DCs","Bcells",
                                      "Naive_CD8_Tcells","Non_Naive_CD8_Tcells","Naive_CD4_Tcells","Non_Naive_CD4_Tcells","Tregs"),
                         "fibroblasts" = c("Fibroblasts"), "endothelial" = c("Endothelial"))

for(cell_type in c("immune", "fibroblasts", "endothelial")){
  if(length(which(markers_df$cell_type %in% cell_types_names[[cell_type]])) != 0){
    print("Human atlas filtering")
    cCRE_38 <- read.table(human_atlas_peaks, header = T, sep = "\t")
    cCRE_38_subset = cCRE_38[-which(cCRE_38$cCDRE_module %in% cell_types_modules[[cell_type]]), ]
    
    regions.1 = GenomicRanges::makeGRangesFromDataFrame(cCRE_38_subset, keep.extra.columns = T)
    regions.2 = Signac::StringToGRanges(markers_df$peak_id[which(markers_df$cell_type %in% cell_types_names[[cell_type]])], sep = c("-", "-"))
    
    overlap = IRanges::subsetByOverlaps(regions.2, regions.1)
    overlapping_markers = c(overlapping_markers, paste0(chrom(overlap), "-", start(overlap), "-", end(overlap)))
  } 
}
markers_df$in_atlas = F
markers_df$in_atlas[which(markers_df$peak_id %in% overlapping_markers)] = T

markers_df$remaining_afterFiltering = !markers_df$in_atlas



# COR filtering ----------------------------------------------------------

print("COR correlation filtering")

markers_list = list()
sapply(unique(markers_df$cell_type), function(i){ markers_list[[i]] <<- markers_df$peak_id[which(markers_df$cell_type == i & markers_df$in_atlas == F)]  })

# Load TCGA reference data
ATAC_counts_norm <- readRDS(TCGA_data)
markers_list_updated_PBMC = markers_list 
markers_list_updated_TCGA = check_correlation(markers_list = markers_list, ATAC_ref = ATAC_counts_norm, min_markers = 2, cor_thr = 0.9, use_caret = T)

markers_df$remaining_afterFiltering_PBMC = markers_df$peak_id %in% unique(unlist(markers_list_updated_PBMC))
markers_df$remaining_afterFiltering_TCGA = markers_df$peak_id %in% unique(unlist(markers_list_updated_TCGA))


# Save PBMC profiles after filtering -------------------------------------------

profiles_list_save <- profiles_list
profiles_list$refProfiles = profiles_list_save$refProfiles[,grep('Macrophages|Fibroblasts|Endothelial',colnames(profiles_list_save$refProfiles),invert = T)]
profiles_list$refProfiles.var = profiles_list_save$refProfiles.var[,grep('Macrophages|Fibroblasts|Endothelial',colnames(profiles_list_save$refProfiles.var),invert = T)]

profiles_list$marker_peaks = markers_df$peak_id[which(!(markers_df$cell_type %in% c("Macrophages", "Fibroblasts", "Endothelial")) &
                                                        markers_df$remaining_afterFiltering == T)]

save(profiles_list, markers_df, file = paste0("PBMC_profile.Rdata"))


# Save TME profile --------------------------------------------------------

profiles_list$refProfiles = profiles_list_save$refProfiles[,grep('Monocytes',colnames(profiles_list_save$refProfiles),invert = T)]
profiles_list$refProfiles.var = profiles_list_save$refProfiles.var[,grep('Monocytes',colnames(profiles_list_save$refProfiles.var),invert = T)]
profiles_list$marker_peaks = markers_df$peak_id[which(!(markers_df$cell_type %in% c("Monocytes")) &
                                                        markers_df$remaining_afterFiltering_TCGA == T)]

save(profiles_list, markers_df, file = paste0("Tumor_profile.Rdata"))



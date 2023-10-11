
# Libraries ---------------------------------------------------------------

library(data.table)
library(Signac)
library(openxlsx)
library(GenomicRanges)
library(scales)
library(rtracklayer)

# Functions ---------------------------------------------------------------
StringToGRanges <- function(regions, sep = c("-", "-")) {
  ranges.df <- data.frame(ranges = regions)
  ranges.df <- tidyr::separate(
    data = ranges.df,
    col = "ranges",
    sep = paste0(sep[[1]], "|", sep[[2]]),
    into = c("chr", "start", "end")
  )
  granges <- makeGRangesFromDataFrame(df = ranges.df)
  return(granges)
}

# Parameters --------------------------------------------------------------

DA_res_file = as.character(commandArgs(TRUE)[1]) 
da_method = unlist(strsplit(fs::path_ext_remove(DA_res_file),"Peaks_pairwise_"))[[2]]

# Load DEG results  --------------------------------------------------------
DEG = read.table(DA_res_file, header = T)

cells_to_consider = unique(DEG$cell_type1)
cells_to_consider = cells_to_consider[which(cells_to_consider %in% c("Macrophages","CD8_Tcells","NK","Endothelial","Monocytes","CD4_Tcells","Neutrophils","DCs","Bcells","Fibroblasts",
                                                                     "Naive_CD4_Tcells","Non_Naive_CD4_Tcells","Tregs","Naive_CD8_Tcells","Non_Naive_CD8_Tcells"))] 

# Process DEG to get markers ----------------------------------------------
logFC_thr = 0.5
nb_markers = 200

DEG_pos = DEG[which(DEG$avg_log2FC > logFC_thr), ]
dim(DEG_pos)

peaks = unique(DEG_pos$gene)
markers_list=list()
for(i in cells_to_consider){
  print(i)
  DA_summary = as.data.frame(matrix(0, nrow = length(peaks), ncol = length(unique(DEG$cell_type1))))
  
  rownames(DA_summary) = peaks
  colnames(DA_summary) = unique(DEG$cell_type1)
  for(j in colnames(DA_summary)[which(colnames(DA_summary) != i)]){
    DA_summary[peaks[which(peaks %in% DEG_pos$gene[which( DEG_pos$cell_type1==i & DEG_pos$cell_type2==j)])],j] = 1
  }
  DA_summary$nb_diff = rowSums(DA_summary)
  
  nb_diff_var = table(DA_summary$nb_diff)
  print(nb_diff_var)
  if(nb_diff_var[length(nb_diff_var)]<10){
    nb_diff_keep = as.numeric(names(nb_diff_var)[c(length(nb_diff_var), (length(nb_diff_var)-1))])
  }else{
    nb_diff_keep = as.numeric(names(nb_diff_var)[c(length(nb_diff_var))])
  }
  
  DEG_subset = DEG_pos[which(DEG_pos$gene %in% rownames(DA_summary)[which(DA_summary$nb_diff %in% nb_diff_keep)] &
                               DEG_pos$cell_type1==i),]
  peaks_subset = unique(DEG_subset$gene)
  
  max_p = aggregate(p_val_adj ~ gene, DEG_subset, max)
  rownames(max_p) = max_p$gene
  
  markers_list[[i]] = max_p$gene[order(max_p$p_val_adj, decreasing = F)][1:min(nrow(max_p), nb_markers)]
}


# save reference signatures
library(plyr)
markers_df <- ldply(markers_list, data.frame)
colnames(markers_df)=c("cell_type","peak_id")

write.table(markers_df, file = paste0("markers_",da_method,".txt"), quote = F, row.names = F, col.names = T, sep = "\t")


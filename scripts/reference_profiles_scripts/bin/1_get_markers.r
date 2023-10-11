
# Libraries ---------------------------------------------------------------
library(limma)
library(edgeR)


# Parameters --------------------------------------------------------------
counts_path = as.character(commandArgs(TRUE)[1])
metadata_path = as.character(commandArgs(TRUE)[2]) 
major_groups = as.logical(commandArgs(TRUE)[3])
ENCODE_counts_path = as.character(commandArgs(TRUE)[4])
ENCODE_metadata_path = as.character(commandArgs(TRUE)[5])
cross_validation_file = as.character(commandArgs(TRUE)[6])
samples_selection_id = as.numeric(commandArgs(TRUE)[7])

# Load reference data to learn the markers --------------------------------

counts = read.table(counts_path, header = T, sep = "\t", row.names = 1) 
metadata = read.table(metadata_path, header = T, sep = "\t") 

if(major_groups){
  metadata$groups = plyr::revalue(metadata$groups,replace = c("Naive_CD4_Tcells"="CD4_Tcells","Non_Naive_CD4_Tcells"="CD4_Tcells", "Tregs"="CD4_Tcells" ,
                                                               "Naive_CD8_Tcells"="CD8_Tcells" ,"Non_Naive_CD8_Tcells"="CD8_Tcells" ))
}else{
  if(length(which(metadata$groups %in% c("CD4_Tcells", "CD8_Tcells")) != 0)){
    counts = counts[,-which(metadata$groups %in% c("CD4_Tcells", "CD8_Tcells"))]
    metadata = metadata[-which(metadata$groups %in% c("CD4_Tcells", "CD8_Tcells")),] 
  } 
} 

if(samples_selection_id != 0){
  samples_subset = read.table(cross_validation_file, header = T, sep = "\t")
  counts = counts[, which(metadata$sample %in% samples_subset$SRA_ID )] 
  metadata = metadata[which(metadata$sample %in% samples_subset$SRA_ID), ] 
  samples_selection_id = paste0("-cv",samples_selection_id)
} else{
  samples_selection_id = "-all"
} 



# Load ENTEX counts -------------------------------------------------------
entex_counts = read.table(ENCODE_counts_path, header = T, sep = "\t", row.names = 1) 
# Load samples metadata
metadata_entex = read.table(ENCODE_metadata_path, header = T, sep = "\t")
rownames(metadata_entex) = metadata_entex$sample
metadata_entex = metadata_entex[colnames(entex_counts),] 

entex_counts = entex_counts[rownames(counts),] 

counts = cbind(counts, entex_counts)
metadata = rbind(metadata, metadata_entex)


# Run differential analysis -----------------------------------------------
cell_types <- metadata$groups

NormFactor <- calcNormFactors(object = counts, method = "TMM")
d0 <- DGEList(counts, norm.factors = NormFactor)
drop <- which(apply(edgeR::cpm(d0), 1, max) < 5)
if(length(drop)!=0){
  d <- d0[-drop,] 
  dim(d)
}else{
  d <- d0
}
mm <- model.matrix( ~ 0 + cell_types)
pdf("voom_plot.pdf")
y <- voom(d, design = mm, plot = T) 
dev.off()
limma_fit <- lmFit(y, mm)

print(head(coef(limma_fit)))
cells_to_consider <- unique(cell_types)
deg_tests_limma <- vector()
for(i in 1:length(cells_to_consider)){
  c_type <- as.character(cells_to_consider[i])
  print(c_type)
  markers_list_eachComp = list()
  for(j in cells_to_consider[which(cells_to_consider != c_type)]){
    contrast <- c(rep(0, length(cells_to_consider)))
    names(contrast) <- c(sort(cells_to_consider)) 
    
    contrast[c_type] <- 1
    contrast[j] <- -1
    contr <- makeContrasts(contrast, levels = colnames(coef(limma_fit)))
    tmp <- contrasts.fit(limma_fit, contr)
    tmp <- eBayes(tmp)
    
    top.table <- topTable(tmp, sort.by = "P", n = Inf)
    print(length(which(top.table$adj.P.Val < 0.01)))
    da_tests <- top.table[which(top.table$adj.P.Val < 0.3 & top.table$logFC >= 0.2), ]
    da_tests$cell_type1 <- rep(c_type, nrow(da_tests))
    da_tests$cell_type2 <- rep(j, nrow(da_tests))
    da_tests$gene <- rownames(da_tests)
    deg_tests_limma <- rbind(deg_tests_limma, da_tests)
  }
}
head(deg_tests_limma)

colnames(deg_tests_limma)[c(1, 5)] = c("avg_log2FC", "p_val_adj")

data.table::fwrite(deg_tests_limma,
                   file = paste0("Peaks_pairwise_limma", samples_selection_id, ".txt"),
                   quote = F, row.names = F, col.names = T, sep = "\t")





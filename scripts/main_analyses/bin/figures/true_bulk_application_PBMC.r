library(EPICATAC)
library(ggplot2)
library(tidyr)
library(rtracklayer)
library(dplyr)
library(ComplexHeatmap)
library(gridExtra)

# Function ----------------------------------------------------------------

get_gene_length = function(gtf){
  GTF <- import.gff(paste0(gtf), format="gtf", feature.type="exon")
  gene_matching=data.frame(ensembl_id=GTF$gene_id,gene_name=GTF$gene_name)
  gene_matching=unique(gene_matching)
  grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
  reducedGTF <- unlist(grl, use.names=T)
  elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementNROWS(grl))
  elementMetadata(reducedGTF)$gene_name <- rep(names(grl), elementNROWS(grl))
  elementMetadata(reducedGTF)$widths <- width(reducedGTF)
  
  length_df = as.data.frame(elementMetadata(reducedGTF))  %>% group_by(gene_id) %>% summarise(sum(widths))
  length_df = as.data.frame(length_df)
  rownames(length_df) = length_df$gene_id
  gene_matching$gene_length = length_df[gene_matching$ensembl_id,2]
  colnames(gene_matching)=c("ensembl_gene_id_1","hgnc_symbol","gene_length")
  gene_matching$ensembl_gene_id=sapply(gene_matching$ensembl_gene_id_1,function(i) unlist(strsplit(i,"[.]"))[1])
  summary(duplicated(gene_matching$ensembl_gene_id))
  
  gene_matching_unique = aggregate(gene_length ~ hgnc_symbol, data = gene_matching, FUN =  mean)
  gene_matching_subset = gene_matching[-which(duplicated(gene_matching$hgnc_symbol)), ]
  rownames(gene_matching_subset) = gene_matching_subset$hgnc_symbol
  gene_matching_unique$ensembl_gene_id = gene_matching_subset[gene_matching_unique$hgnc_symbol, "ensembl_gene_id"] 
  
  return(gene_matching_subset)
  
}


# Parameters --------------------------------------------------------------
data_path = as.character(commandArgs(TRUE)[1])  
gtf_file = as.character(commandArgs(TRUE)[2])  
fig_path = as.character(commandArgs(TRUE)[3])

# Load ATAC, RNA and cytometry data ---------------------------------------

atac_counts <- read.delim(paste0(data_path, "GSE193140_atac_counts_qval_0.05.txt.gz"), sep = " ")
dim(atac_counts)
rownames(atac_counts) <- paste0(atac_counts$chr, ":", atac_counts$start, "-", atac_counts$end)
atac_counts <- atac_counts[, -c(1:3)]

atac_samples <- sapply(colnames(atac_counts), function(i) paste(strsplit(i, "_")[[1]][1:2], collapse = "_"))
atac_counts <- atac_counts[, -which(duplicated(atac_samples))]
colnames(atac_counts) <- atac_samples[-which(duplicated(atac_samples))]

rna_counts <- read.delim(paste0(data_path, "GSE193141_rna_counts_011622.csv.gz"), sep = ",")
dim(rna_counts)
rna_counts <- rna_counts[which(!is.na(rna_counts$GeneSymbol)), ]
rna_counts <- rna_counts[-which(duplicated(rna_counts$GeneSymbol)), ]
rownames(rna_counts) <- rna_counts$GeneSymbol
rna_counts <- rna_counts[, -c(1:2)]
rna_samples <- sapply(colnames(rna_counts), function(i) strsplit(i, "RNA_")[[1]][2])
colnames(rna_counts) <- rna_samples

summary(colnames(atac_counts) %in% colnames(rna_counts))
summary(colnames(rna_counts) %in% colnames(atac_counts))
rna_counts <- rna_counts[, colnames(atac_counts)]

metadata <- read.table(paste0(data_path, "11357_2023_986_MOESM7_ESM.tsv"), sep = "\t")
summary(metadata$PassesQC)
metadata <- metadata[which(metadata$PassesQC == T), ]
rownames(metadata) <- metadata$Subject
colnames(metadata) <- plyr::revalue(colnames(metadata), c("Granulocytes" = "Neutrophils",
                                                          "B_Cells" = "Bcells",
                                                          "NK_Cells" = "NKcells",
                                                          "CD4_T_Cells" = "CD4_Tcells",
                                                          "CD8_T_Cells" = "CD8_Tcells"))
metadata <- metadata[metadata$Subject %in% colnames(atac_counts), ]

# Normalize proportions so that they sum to 1 ------------------------------

cytometry <- metadata[, c("Monocytes","Neutrophils","Bcells","NKcells","CD4_Tcells","CD8_Tcells", "T_Cells", "Lymphocytes")]
cytometry <- cytometry[-which(rownames(cytometry) == "CR_118"),] # CR_118 has NA values

cytometry[,c("CD4_Tcells","CD8_Tcells")] <-  cytometry[,c("CD4_Tcells","CD8_Tcells")] / rowSums(cytometry[,c("CD4_Tcells","CD8_Tcells")]) * cytometry$T_Cells
cytometry[,c("Bcells","NKcells", "CD4_Tcells","CD8_Tcells")] <- cytometry[,c("Bcells","NKcells", "CD4_Tcells","CD8_Tcells")] / rowSums(cytometry[,c("Bcells","NKcells", "CD4_Tcells","CD8_Tcells")]) * cytometry$Lymphocytes
cytometry <- cytometry[, c("Monocytes","Neutrophils","Bcells","NKcells","CD4_Tcells","CD8_Tcells")] 

cytometry = as.data.frame(t(cytometry))/100
cytometry = t(as.data.frame(t(t(cytometry)/colSums(cytometry,na.rm = T))))

summary(colnames(atac_counts) %in% metadata$Subject)
summary(colnames(rna_counts) %in% metadata$Subject)
atac_counts <- atac_counts[, metadata$Subject]
rna_counts <- rna_counts[, metadata$Subject]


#  run EPIC-ATAC on atacseq samples ---------------------------------------

peaks_gr = Signac::StringToGRanges(rownames(atac_counts),sep = c(":","-"))
region.length = peaks_gr@ranges@width
names(region.length) = rownames(atac_counts)
atac_counts <-  atac_counts / region.length[rownames(atac_counts)]
atac_counts <- t(t(atac_counts) * 1e6 / colSums(atac_counts))

epicatac_res <- EPICATAC::EPIC_ATAC(atac_counts, reference = EPICATAC::atacRef_PBMC, withOtherCells = F,
                                    ATAC = T, genome_version = "hg38", nb_iter = 10000)
epicatac_res <- as.data.frame(epicatac_res$cellFractions)
colnames(epicatac_res)[2] <- "NKcells"
epicatac_res <- epicatac_res[-which(rownames(epicatac_res) %in% c("CR_118")),]


# run EPIC-RNA on RNAseq samples -----------------------------------------

gft_file_10x = paste0(gtf_file)
gene_matching <- get_gene_length(gft_file_10x)
summary(rownames(rna_counts) %in% rownames(gene_matching))
rna_counts <- rna_counts[which(rownames(rna_counts) %in% rownames(gene_matching)),] 
summary(EPICATAC::BRef$sigGenes %in% rownames(rna_counts))

gene_lengths <- gene_matching[rownames(rna_counts), "gene_length"]
rna_counts <- rna_counts / gene_lengths
rna_counts <- t(t(rna_counts) * 1e6 / colSums(rna_counts))

epicrna_res <- EPICATAC::EPIC_ATAC(rna_counts, reference = EPICATAC::BRef, withOtherCells = F, ATAC = F)
epicrna_res <- as.data.frame(epicrna_res$cellFractions)
epicrna_res <- epicrna_res[-which(rownames(epicrna_res) == "CR_118"),]


# Compare predictions -----------------------------------------------------

rmse_matrix = cor_matrix = cor_text = matrix(nrow = 3, ncol = 7,
                                             dimnames = list(c("EPIC-ATAC", "EPIC-RNA","EPIC-ATAC/RNA"),
                                                             c("Bcells", "CD4_Tcells", "CD8_Tcells",  "Neutrophils", "NKcells", "Monocytes", "All")))

for(cell_type in c("Bcells", "CD4_Tcells", "CD8_Tcells",  "Neutrophils", "NKcells", "Monocytes")){ 
  
  method = "EPIC-ATAC"
  rmse_matrix[method, cell_type] = round(Metrics::rmse(epicatac_res[, cell_type], cytometry[rownames(epicatac_res), cell_type]),3)
  
  correlation_test = cor.test(epicatac_res[, cell_type], cytometry[rownames(epicatac_res), cell_type], alternative = "two.sided", method = "pearson")
  if(correlation_test$p.value < 0.001){pval = "***"} else if (correlation_test$p.value < 0.01){pval = "**"} else if (correlation_test$p.value <= 0.05){ pval = "*"} else{pval = ""}
  cor_text[method, cell_type] = paste0(round(correlation_test$estimate,3), "\n", pval)
  cor_matrix[method, cell_type] = round(correlation_test$estimate,3)
  
  method = "EPIC-RNA"
  rmse_matrix[method, cell_type] = round(Metrics::rmse(epicrna_res[, cell_type], cytometry[rownames(epicrna_res), cell_type]),3)
  
  correlation_test = cor.test(epicrna_res[, cell_type], cytometry[rownames(epicrna_res), cell_type], alternative = "two.sided", method = "pearson")
  if(correlation_test$p.value < 0.001){pval = "***"} else if (correlation_test$p.value < 0.01){pval = "**"} else if (correlation_test$p.value <= 0.05){ pval = "*"} else{pval = ""}
  cor_text[method, cell_type] = paste0(round(correlation_test$estimate,3), "\n", pval)
  cor_matrix[method, cell_type] = round(correlation_test$estimate,3)
  
  method = "EPIC-ATAC/RNA"
  rmse_matrix[method, cell_type] = round(Metrics::rmse( rowMeans(cbind(epicrna_res[, cell_type], epicatac_res[rownames(epicrna_res), cell_type])), cytometry[rownames(epicrna_res), cell_type]),3)
  
  correlation_test = cor.test(rowMeans(cbind(epicrna_res[, cell_type], epicatac_res[rownames(epicrna_res), cell_type])), cytometry[rownames(epicrna_res), cell_type], 
                              alternative = "two.sided", method = "pearson")
  if(correlation_test$p.value < 0.001){pval = "***"} else if (correlation_test$p.value < 0.01){pval = "**"} else if (correlation_test$p.value <= 0.05){ pval = "*"} else{pval = ""}
  cor_text[method, cell_type] = paste0(round(correlation_test$estimate,3), "\n", pval)
  cor_matrix[method, cell_type] = round(correlation_test$estimate,3)
  
} 


all_pred_rna <- unlist(lapply(c("Bcells", "CD4_Tcells", "CD8_Tcells",  "Neutrophils", "NKcells", "Monocytes"), function(i) epicrna_res[, i]))
all_pred_atac <- unlist(lapply(c("Bcells", "CD4_Tcells", "CD8_Tcells",  "Neutrophils", "NKcells", "Monocytes"), function(i) epicatac_res[rownames(epicatac_res), i]))
all_truth <- unlist(lapply(c("Bcells", "CD4_Tcells", "CD8_Tcells",  "Neutrophils", "NKcells", "Monocytes"), function(i) cytometry[rownames(epicatac_res), i]))
cell_type = "All"
method = "EPIC-ATAC"
rmse_matrix[method, cell_type] = round(Metrics::rmse(all_pred_atac, all_truth),3)
correlation_test = cor.test(all_pred_atac, all_truth, alternative = "two.sided", method = "pearson")
if(correlation_test$p.value < 0.001){pval = "***"} else if (correlation_test$p.value < 0.01){pval = "**"} else if (correlation_test$p.value <= 0.05){ pval = "*"} else{pval = ""}
cor_text[method, cell_type] = paste0(round(correlation_test$estimate,3), "\n", pval)
cor_matrix[method, cell_type] = round(correlation_test$estimate,3)

method = "EPIC-RNA"
rmse_matrix[method, cell_type] = round(Metrics::rmse(all_pred_rna, all_truth),3)
correlation_test = cor.test(all_pred_rna, all_truth, alternative = "two.sided", method = "pearson")
if(correlation_test$p.value < 0.001){pval = "***"} else if (correlation_test$p.value < 0.01){pval = "**"} else if (correlation_test$p.value <= 0.05){ pval = "*"} else{pval = ""}
cor_text[method, cell_type] = paste0(round(correlation_test$estimate,3), "\n", pval)
cor_matrix[method, cell_type] = round(correlation_test$estimate,3)

method = "EPIC-ATAC/RNA"
rmse_matrix[method, cell_type] = round(Metrics::rmse(rowMeans(cbind(all_pred_rna, all_pred_atac)), all_truth),3)
correlation_test = cor.test(rowMeans(cbind(all_pred_rna, all_pred_atac)), all_truth, alternative = "two.sided", method = "pearson")
if(correlation_test$p.value < 0.001){pval = "***"} else if (correlation_test$p.value < 0.01){pval = "**"} else if (correlation_test$p.value <= 0.05){ pval = "*"} else{pval = ""}
cor_text[method, cell_type] = paste0(round(correlation_test$estimate,3), "\n", pval)
cor_matrix[method, cell_type] = round(correlation_test$estimate,3)


axis_text_size = 11; text_size = 9
col_fun <- colorRampPalette(c("#FC2232","#FFD85B","white"))(256)
cn = colnames(rmse_matrix)
ht_rmse <- Heatmap(rmse_matrix, col = col_fun, show_row_names = T, row_names_gp = gpar(fontsize = axis_text_size), column_names_gp = gpar(fontsize = axis_text_size),
                   show_column_names = F, row_names_side = "left",
                   bottom_annotation = HeatmapAnnotation(
                     text = anno_text(cn, rot = 45)
                   ),
                   row_title= NULL, column_title = "RMSE", cluster_rows  = F, cluster_columns = F,
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(rmse_matrix[i, j], x, y, gp = gpar(fontsize = text_size))
                   },
                   na_col = "snow3", 
                   column_split = 1:ncol(rmse_matrix), border = TRUE ,heatmap_legend_param = list(title = "RMSE"))

col_fun <- colorRampPalette(c("white","#FFD85B","#FC2232"))(256)
cn = colnames(cor_matrix)
ht_cor <- Heatmap(cor_matrix, col = col_fun, show_row_names = T, row_names_gp = gpar(fontsize = axis_text_size), column_names_gp = gpar(fontsize = axis_text_size),
                  show_column_names = F, row_names_side = "left",
                  bottom_annotation = HeatmapAnnotation(
                    text = anno_text(cn, rot = 45)
                  ),
                  row_title= NULL, column_title = "Correlation", cluster_rows  = F, cluster_columns = F,
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(cor_text[i, j], x, y, gp = gpar(fontsize = text_size))
                  },
                  na_col = "snow3", column_split = 1:ncol(cor_matrix), border = TRUE ,heatmap_legend_param = list(title = "Correlation"))

combined_heatmap <- ht_cor + ht_rmse

pdf(paste0(fig_path, "Fig7_panelC.pdf"), w = 9, h = 3)
draw(combined_heatmap, ht_gap = unit(2,"cm"))
dev.off()


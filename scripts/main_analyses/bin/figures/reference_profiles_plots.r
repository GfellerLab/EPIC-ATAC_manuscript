
# Libraries ---------------------------------------------------------------
library(ggplot2)
library(plyr)
library(dplyr)
require(pheatmap)
require(RColorBrewer)
library(gridExtra)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(ComplexHeatmap)
library(readxl)
library(mclust)
library(circlize)
library(RColorBrewer)
library(plot3D)

# Functions ---------------------------------------------------------------

StringToGRanges <- function(regions, sep = c("-", "-"), ...){
  ranges.df <- data.frame(ranges = regions)
  ranges.df <- tidyr::separate(data = ranges.df, col = "ranges", sep = paste0(sep[[1]], "|", sep[[2]]), into = c("chr", "start", "end"))
  rownames(ranges.df) <- regions
  granges <- GenomicRanges::makeGRangesFromDataFrame(df = ranges.df, ...)
  return(granges)
}

# Parameters --------------------------------------------------------------
descriminative_colors = c("Bcells"="darkgreen","CD4_Tcells"="#3288bdff","CD8_Tcells"="#4daf4aff","Dendritic"="#984ea3ff", "DCs"="#984ea3ff",  
                          "Macrophages"="#ff7f00ff","Monocytes"="#ffd92fff","Myeloid"="#a0451fff", "DCs + Macrophages"="#984ea3ff", "NKcells"="#70c4b4ff","NK"="#70c4b4ff",
                          "NK + T cells"="#bea0ccff", "Tcells"="#660033","Uncharacterized"="#d90017ff", "Neutrophils"="black","Endothelial"="#a0451fff","Fibroblasts"="gray",
                          "Naive_CD4_Tcells"="#bea0ccff","Non_Naive_CD4_Tcells"="#99FFFF", "Memory_and_Helper_CD4" = "#99FFFF","Tregs"="#3288bdff", "Naive_CD8_Tcells"="#660033","Non_Naive_CD8_Tcells"="#4daf4aff",
                          "bladder" = '#660033', "breast" = "#CC0000", "colon" = "#99FFFF", "liver" = "#66CC00", "lung" = "#FFFF66", "ovary" = "#FFB266", "pancreas" = "#333300", "thyroid" = "#6666FF")

withSuptypes = as.logical(commandArgs(TRUE)[1])
data_path = as.character(commandArgs(TRUE)[2])
fig_path = as.character(commandArgs(TRUE)[3])

fasta_file <- paste0(data_path, "annotation_files/hg38.fa")
ega_counts_path <- paste0(data_path, "/markers_validation/validation_bulk_counts.txt")
all_sra_counts_path <- paste0(data_path, "/markers_validation/sra_counts.txt")
ENTEX_path <- paste0(data_path, "/markers_validation/encode_counts.txt") 
sup_tables <- paste0(data_path, "/Supplementary_tables.xlsx")
human_atlas_counts <- paste0(data_path, "/markers_validation/human_atlas_our_peaks_counts.rds")
human_atlas_metadata <- paste0(data_path,"/markers_validation/GSE184462_metadata.tsv")

if(withSuptypes){
  fig_path <- paste0(fig_path,"/with_Subtypes/")
}
dir.create(fig_path,recursive = T)


# Load profiles and markers -----------------------------------------------
if(!withSuptypes){
  # PBMC_profile 
  load(paste0(data_path, "/reference_profiles/PBMC_profile_noSubtypes.Rdata"))
  
  markers <- profiles_list$marker_peaks
  markers_df <- markers_df[which(markers_df$peak_id %in% markers & markers_df$cell_type %in% colnames(profiles_list$refProfiles)), ]
  PBMC_profile = profiles_list
  PBMC_profile$sigGenes <- markers
  PBMC_markers = markers_df[which(!markers_df$cell_type %in% c("Endothelial","Macrophages","Fibroblasts")),] 
  
  
  # Tumor profile
  load(paste0(data_path, "/reference_profiles/TME_profile_noSubtypes.Rdata"))
  
  markers <- profiles_list$marker_peaks
  markers_df <- markers_df[which(markers_df$peak_id %in% markers & markers_df$cell_type %in% colnames(profiles_list$refProfiles)), ]
  
  Tumor_profile = profiles_list
  Tumor_profile$sigGenes <- markers
  Tumor_markers = markers_df[which(!markers_df$cell_type %in% c("Monocytes")),] 
}else{
  # PBMC_profile 
  load(paste0(data_path, "/reference_profiles/PBMC_profile_withSubtypes.Rdata"))
  markers <- profiles_list$marker_peaks
  markers_df <- markers_df[which(markers_df$peak_id %in% markers & markers_df$cell_type %in% colnames(profiles_list$refProfiles)), ]
  PBMC_profile = profiles_list
  PBMC_profile$sigGenes <- markers
  PBMC_markers = markers_df[which(!markers_df$cell_type %in% c("Endothelial","Macrophages","Fibroblasts")),] 
  
  # Tumor profile
  load(paste0(data_path, "/reference_profiles/TME_profile_withSubtypes.Rdata"))
  markers <- profiles_list$marker_peaks
  markers_df <- markers_df[which(markers_df$peak_id %in% markers & markers_df$cell_type %in% colnames(profiles_list$refProfiles)), ]
  Tumor_profile = profiles_list
  Tumor_profile$sigGenes <- markers
  Tumor_markers = markers_df[which(!markers_df$cell_type %in% c("Monocytes")),] 
} 

markers = unique(c(PBMC_profile$sigGenes, Tumor_profile$sigGenes))
markers_df = unique(rbind(PBMC_markers, Tumor_markers))

# Load SRA counts ---------------------------------------------------------
sra_counts = read.table(all_sra_counts_path, skip = 1, header = T)
rownames_save = sra_counts$Geneid
sra_counts = sra_counts[, c(7:ncol(sra_counts))]
rownames(sra_counts) = rownames_save
colnames(sra_counts) = sapply(colnames(sra_counts), function(i) unlist(strsplit(i, '_sort_dedup.bam'))[[1]])


# Load EGA counts ---------------------------------------------------------
ega_counts = read.table(ega_counts_path, header = T, sep = "\t", row.names = 1) 
ega_counts = ega_counts[, c(6:ncol(ega_counts))]
colnames(ega_counts) = sapply(colnames(ega_counts), function(i) unlist(strsplit(i, '_sort_dedup.bam'))[[1]])
colnames(ega_counts) = sapply(colnames(ega_counts), function(i) unlist(strsplit(i, '[.]bam_files[.]'))[[2]])


# Load ENTEX counts -------------------------------------------------------
entex_counts = read.table(ENTEX_path, header = T, sep = "\t")
dim(entex_counts)


# Combine all samples together --------------------------------------------
combined_samples = as.matrix(cbind(sra_counts, ega_counts[rownames(sra_counts), ], entex_counts[rownames(sra_counts),]))
rownames(combined_samples) = gsub(":", "-", rownames(combined_samples))

combined_samples = combined_samples[grep("chrY|chrUn|random",rownames(combined_samples),invert = T),]
dim(combined_samples)


# Load metadata -----------------------------------------------------------
combined_metadata <- as.data.frame(read_excel(sup_tables, sheet = 2))
rownames(combined_metadata) <- combined_metadata$Sample_id
combined_metadata$annotation = combined_metadata$Cell_type

if(!withSuptypes){
  combined_metadata$Cell_type = plyr::revalue(combined_metadata$Cell_type, c("Non_Naive_CD8_Tcells" = "CD8_Tcells",
                                                                             "Naive_CD8_Tcells" = "CD8_Tcells",
                                                                             "Memory_and_Helper_CD4" = "CD4_Tcells",
                                                                             "Non_Naive_CD4_Tcells" = "CD4_Tcells","Tregs" = "CD4_Tcells", "Naive_CD4_Tcells" = "CD4_Tcells"))
} else{
  combined_samples = combined_samples[, -which(colnames(combined_samples) %in% rownames(combined_metadata)[which(combined_metadata$Cell_type %in% c("CD4_Tcells", "CD8_Tcells"))] )] 
} 


combined_samples = combined_samples[, which(colnames(combined_samples) %in% rownames(combined_metadata))]
combined_metadata = combined_metadata[colnames(combined_samples), ]


# Keep only reference samples ----------------------------------------------
subset_samples = combined_samples[, which(colnames(combined_samples) %in% combined_metadata$Sample_id[combined_metadata$Usage == "Reference sample"])]
subset_metadata = combined_metadata[colnames(subset_samples), ]

# Figure 2 A, nb samples per study and cell type 
plot_df = subset_metadata %>% group_by(Study, Cell_type) %>% tally()

axis_text_size = 25
text_size = 25
legend_text_size = 18
my_theme <- theme(axis.text = element_text(size = axis_text_size, face = "bold"),
                  legend.title = element_text(size = text_size, face = "bold"),
                  legend.text = element_text(size = legend_text_size),
                  axis.title = element_text(size = text_size, face = "bold"),
                  strip.text = element_text(size = text_size, face = "bold"),
                  axis.title.y = element_text(vjust = +1.5),
                  axis.title.x = element_text(vjust = -0.5),
                  legend.key.size = unit(1.5, 'lines')) 

pdf(paste0(fig_path, "fig2_panel_A.pdf"))
study_barplot <- ggplot(data=plot_df,aes(fill=Study, y=n, x=Cell_type)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("#F0A0FF","#0075DC","#993F00","#4C005C","#191919","#005C31","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00"))+ 
  ggtitle("") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + my_theme +
  ylab("Number of samples")+
  guides(fill=guide_legend("Study")) + theme(legend.position="top") +
  xlab("")
print(study_barplot)
dev.off()


# Get TPM like counts
peaks_gr <- StringToGRanges(rownames(subset_samples))
region.length <- peaks_gr@ranges@width
names(region.length) <- names(peaks_gr@ranges)
x <-  subset_samples / region.length[rownames(subset_samples)]
tpm.mat <- t(t(x) * 1e6 / colSums(x))
rm(x)

# generate UMAP/PCA for the reference samples ---------------------------------

gc_norm <- function(input_matrix, grouping_var, fa_file, by_group = F){
  ff <- Rsamtools::FaFile(fa_file)
  gr <- Signac::StringToGRanges(rownames(input_matrix))
  gr = GRanges(seqnames=seqnames(gr), ranges=IRanges(start(gr), end(gr)), strand="*", mcols=data.frame(peakID=rownames(input_matrix)))
  
  peakSeqs <- getSeq(x=ff, gr)
  gcContentPeaks = letterFrequency(peakSeqs, "GC",as.prob=TRUE)[,1]
  
  if(by_group){
    dataNorm <- vector()
    for(group in unique(grouping_var)){
      print(group)
      input_matrix_subset <- input_matrix[, which(grouping_var == group)] 
      dataWithin <- EDASeq::withinLaneNormalization(input_matrix_subset, y = gcContentPeaks, num.bins = 20, which = "full")#adjusting for gc content  
      dataNorm_subset <- EDASeq::betweenLaneNormalization(dataWithin, which = "full")# adjusting for sequencing depth 
      
      dataNorm <- cbind(dataNorm, dataNorm_subset)
    } 
  } else{
    dataWithin = EDASeq::withinLaneNormalization(input_matrix, y = gcContentPeaks, num.bins = 20, which = "full")#adjusting for gc content
    dataNorm = EDASeq::betweenLaneNormalization(dataWithin, which = "full")# adjusting for sequencing depth
  } 
  
  return(dataNorm)
} 

get_silhouette_values <- function(pcaData, meta_data, nb_pc = 10){
  dd <- as.matrix(dist(pcaData[,1:nb_pc]))
  num_factors <- 1:length(unique(meta_data$Study))
  names(num_factors) <- unique(meta_data$Study)
  sil_study <- summary(cluster::silhouette(num_factors[meta_data$Study], dd))$avg.width
  
  num_factors <- 1:length(unique(meta_data$Cell_type))
  names(num_factors) <- unique(meta_data$Cell_type)
  sil_cell_types <- summary(cluster::silhouette(num_factors[meta_data$Cell_type], dd))$avg.width
  
  ari <- mclust::adjustedRandIndex(meta_data$cluster, meta_data$Cell_type)
  sil_text <- paste0("Silhouette_study: ", round(sil_study, 2), "\n Silhouette_cellType: ", round(sil_cell_types, 2),
                     "\n ARI: ", round(ari, 2))
  return(sil_text)
} 
dataNorm <- gc_norm(input_matrix = subset_samples, 
                    grouping_var = subset_metadata$Cell_type, 
                    fa_file = fasta_file, 
                    by_group = F)
dataNorm <- dataNorm[,subset_metadata$Sample_id] 


# PBMC_profiles markers only
pca = prcomp(t(log(dataNorm[PBMC_profile$sigGenes,]+1 )))
print(summary(pca))
pcaData = as.data.frame(pca$x)
rownames(pcaData) = rownames(subset_metadata)

subset_metadata$cluster <- paste0("C", mclust::Mclust(pcaData[, 1:10])$classification )

umap_res = umap::umap(t(log(dataNorm[PBMC_profile$sigGenes,]+1 )), verbose = T, min_dist = 0.6)

umap_dr=umap_res$layout
colnames(umap_dr)=c("UMAP1","UMAP2")

my_theme <- theme_bw() + 
  theme(axis.text=element_text(size=15,face="bold"), 
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        axis.title=element_text(size=15,face="bold")) 

to_plot=merge(umap_dr,subset_metadata[,c("Sample_id","Cell_type","Study", "Sample_type","annotation","cluster")],by.x="row.names",by.y="Sample_id")
umap_refSamples_PBMCmarkers <- ggplot(to_plot, aes(x=UMAP1, y=UMAP2, color=Cell_type)) +
  geom_point(size=1, stroke = 0.3) + 
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(subset_metadata$Cell_type))]) +
  my_theme + 
  guides(col=guide_legend("Cell types", ncol = 1, override.aes = list(size=4)))+ 
  annotation_custom(grobTree(textGrob(get_silhouette_values(pcaData, subset_metadata, nb_pc = 10), x=0.1,  y=0.95, hjust=0,
                                      gp=gpar(col="black", fontsize=15))))


# pca 
to_plot=merge(pcaData,subset_metadata[,c("Sample_id","Cell_type","Study", "Sample_type","annotation","cluster")],by.x="row.names",by.y="Sample_id")
pca_refSamples_PBMCmarkers <- ggplot(to_plot, aes(x=PC1, y=PC2, color=Cell_type)) +
  geom_point(size=2, stroke = 0.3) + 
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(subset_metadata$Cell_type))]) +
  my_theme + 
  guides(col=guide_legend("Cell types", ncol = 1, override.aes = list(size=4)))+ 
  annotation_custom(grobTree(textGrob(get_silhouette_values(pcaData, subset_metadata, nb_pc = 10), x=0.1,  y=0.95, hjust=0,
                                      gp=gpar(col="black", fontsize=15))))

pdf(paste0(fig_path,"fig2_panel_B_sup1.pdf"), w = 12, h = 12)
par(cex.axis = 1.7)
pairs(pcaData[,1:3], col = descriminative_colors[subset_metadata$Cell_type], pch = 20, cex = 2,cex.labels = 5)
dev.off()


pdf(paste0(fig_path, "PCA_3D_PBMCmarkers.pdf"))
col_int <- 1:length(unique(subset_metadata$Cell_type))
names(col_int) <- unique(subset_metadata$Cell_type)
colors <- descriminative_colors[unique(subset_metadata$Cell_type)] #col_int[subset_metadata$Cell_type]


theta <- 40
phi <- 30
# for(theta in seq(0,180,10)){
#   for(phi in seq(0,180,10)){
    
    plot3D::scatter3D(x = pcaData$PC1, y = pcaData$PC2, z = pcaData$PC3, bty="b", #type ="h", 
                            pch = 18, theta = theta , phi = phi,
                            colvar = as.vector(col_int[subset_metadata$Cell_type]), 
                            col = as.vector(colors), alpha = 1,
                            xlab = "PC1", ylab="PC2", zlab = "PC3", main = paste0("theta:", theta, " Phi:", phi))
    print(segments3D(x0 = pcaData$PC1, y0 = pcaData$PC2, z0 = rep(min(pcaData$PC3), length(pcaData$PC3)),
               x1 = pcaData$PC1, y1 = pcaData$PC2, z1 = pcaData$PC3, col = "gray", alpha = 0.2, add = TRUE))
    
    
    # segments3D(x0 = pcaData$PC1, y0 = pcaData$PC2, z0 = rep(min(pcaData$PC3), length(pcaData$PC3)), xlab = "PC1", ylab="PC2", zlab = "PC3", 
    #            x1 = pcaData$PC1, y1 = pcaData$PC2, z1 = pcaData$PC3, col = rgb(0.5, 0.5, 0.5, alpha = 0.5), lwd = 1, theta = theta , phi = phi, main = paste0("theta:", theta, " Phi:", phi))
    # 
    # print(plot3D::scatter3D(x = pcaData$PC1, y = pcaData$PC2, z = pcaData$PC3, bty="b", #type ="h", 
    #                         pch = 18, theta = theta , phi = phi,
    #                         colvar = as.vector(col_int[subset_metadata$Cell_type]), 
    #                         col = as.vector(colors), 
    #                         alpha = 1, add = TRUE))
#   }
# }
dev.off()


# Tumor_profiles markers only
pca = prcomp(t(log(dataNorm[Tumor_profile$sigGenes,]+1 )))
print(summary(pca))
pcaData = as.data.frame(pca$x)
rownames(pcaData) = rownames(subset_metadata)
subset_metadata$cluster <- paste0("C", mclust::Mclust(pcaData[, 1:10])$classification )#paste0("C",kmeans(pcaData[, 1:10], centers = 10)$cluster)

umap_res = umap::umap(t(log(dataNorm[Tumor_profile$sigGenes,]+1 )), verbose = T, min_dist = 0.6)

umap_dr=umap_res$layout
colnames(umap_dr)=c("UMAP1","UMAP2")

to_plot=merge(umap_dr,subset_metadata[,c("Sample_id","Cell_type","Study", "Sample_type","annotation","cluster")],by.x="row.names",by.y="Sample_id")
umap_refSamples_Tumormarkers <- ggplot(to_plot, aes(x=UMAP1, y=UMAP2, color=Cell_type)) +
  geom_point(size=1, stroke = 0.3) + 
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(subset_metadata$Cell_type))]) +
  my_theme + 
  guides(col=guide_legend("Cell types", ncol = 1, override.aes = list(size=4)))+ 
  annotation_custom(grobTree(textGrob(get_silhouette_values(pcaData, subset_metadata, nb_pc = 10), x=0.1,  y=0.95, hjust=0,
                                      gp=gpar(col="black", fontsize=15))))

# pca 
to_plot=merge(pcaData,subset_metadata[,c("Sample_id","Cell_type","Study", "Sample_type","annotation","cluster")],by.x="row.names",by.y="Sample_id")
pca_refSamples_Tumormarkers <- ggplot(to_plot, aes(x=PC1, y=PC2, color=Cell_type)) +
  geom_point(size=2, stroke = 0.3) + 
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% unique(subset_metadata$Cell_type))]) +
  my_theme + 
  guides(col=guide_legend("Cell types", ncol = 1, override.aes = list(size=4)))+ 
  annotation_custom(grobTree(textGrob(get_silhouette_values(pcaData, subset_metadata, nb_pc = 10), x=0.1,  y=0.95, hjust=0,
                                      gp=gpar(col="black", fontsize=15))))

pdf(paste0(fig_path,"fig2_panel_B_sup2.pdf"), w = 12, h = 12)
par(cex.axis = 1.7)
pairs(pcaData[,1:3], col = descriminative_colors[subset_metadata$Cell_type], pch = 20, cex = 2,cex.labels = 5)
dev.off()

pdf(paste0(fig_path,"fig2_panel_B_UMAP.pdf"), w = 12, h = 4)
cowplot::plot_grid(
  umap_refSamples_PBMCmarkers + theme(legend.position = "none"), NULL,
  umap_refSamples_Tumormarkers + theme(legend.position = "none"), NULL, cowplot::get_legend(umap_refSamples_Tumormarkers),
  labels = c("", "", ""), ncol = 5,  label_size = 20, rel_widths = c(0.8, 0.05, 0.8, 0.05,0.5)
)
dev.off()

pdf(paste0(fig_path,"fig2_panel_B_PCA.pdf"), w = 12, h = 4)
cowplot::plot_grid(
  pca_refSamples_PBMCmarkers + theme(legend.position = "none"), NULL,
  pca_refSamples_Tumormarkers + theme(legend.position = "none"), NULL, cowplot::get_legend(pca_refSamples_Tumormarkers),
  labels = c("", "", ""), ncol = 5,  label_size = 20, rel_widths = c(0.8, 0.05, 0.8, 0.05,0.5)
)
dev.off()

pdf(paste0(fig_path, "PCA_3D_TMEmarkers.pdf"))
col_int <- 1:length(unique(subset_metadata$Cell_type))
names(col_int) <- unique(subset_metadata$Cell_type)
colors <- descriminative_colors[unique(subset_metadata$Cell_type)] 

theta <- 130
phi <- 30
# for(theta in seq(0,180,10)){
#   for(phi in seq(0,180,10)){

    plot3D::scatter3D(x = pcaData$PC1, y = pcaData$PC2, z = pcaData$PC3, bty="b", 
                      pch = 18, theta = theta , phi = phi,
                      colvar = as.vector(col_int[subset_metadata$Cell_type]), 
                      col = as.vector(colors), alpha = 1,
                      xlab = "PC1", ylab="PC2", zlab = "PC3", main = paste0("theta:", theta, " Phi:", phi))
    print(segments3D(x0 = pcaData$PC1, y0 = pcaData$PC2, z0 = rep(min(pcaData$PC3), length(pcaData$PC3)),
                     x1 = pcaData$PC1, y1 = pcaData$PC2, z1 = pcaData$PC3, col = "gray", alpha = 0.2, add = TRUE))
    
    # segments3D(x0 = pcaData$PC1, y0 = pcaData$PC2, z0 = rep(min(pcaData$PC3), length(pcaData$PC3)), xlab = "PC1", ylab="PC2", zlab = "PC3",
    #            x1 = pcaData$PC1, y1 = pcaData$PC2, z1 = pcaData$PC3, col = rgb(0.5, 0.5, 0.5, alpha = 0.5), lwd = 1, theta = theta , phi = phi, main = paste0("theta:", theta, " Phi:", phi))
    # 
    # print(plot3D::scatter3D(x = pcaData$PC1, y = pcaData$PC2, z = pcaData$PC3, bty="b", #type ="h", 
    #                         pch = 18, #theta = theta , phi = phi,
    #                         colvar = as.vector(col_int[subset_metadata$Cell_type]), 
    #                         col = as.vector(colors), 
    #                         alpha = 1, add = TRUE))
#   }
# }

dev.off()


# generate averaged profile heatmap + sample heatmap 
get_markers_heatmaps <- function(counts_data, norm_data, meta_data, peaks, markers_data, profile_heatmap_title, min_max_scaling = F, sample_heatmap = F){
  
  peaks_gr <- StringToGRanges(gsub(":", "-",rownames(counts_data)))
  region.length <- peaks_gr@ranges@width
  names(region.length) <- names(peaks_gr@ranges)
  x <-  counts_data[peaks, ] / region.length[rownames(counts_data[peaks, ])]
  tpm.mat <- t(t(x) * 1e6 / colSums(x))
  rm(x)
  
  averaged_profiles = data.frame(row.names = rownames(tpm.mat))
  for(i in 1:length(unique(meta_data$Cell_type))) {
    c_type = as.vector(unique(meta_data$Cell_type))[i]
    print(c_type)
    cell_index = which(meta_data$Cell_type %in% c_type)
    if(length(cell_index) != 1){
      averaged_profiles = cbind(averaged_profiles, matrixStats::rowMedians(tpm.mat[, cell_index]))
    }else{
      averaged_profiles = cbind(averaged_profiles, tpm.mat[, cell_index])
    }
    colnames(averaged_profiles)[ncol(averaged_profiles)] = c_type
  }
  head(averaged_profiles)
  
  if(withSuptypes){
    cell_types = c("Bcells","Naive_CD4_Tcells", "Memory_and_Helper_CD4", "Tregs","Naive_CD8_Tcells", "Non_Naive_CD8_Tcells","NK","Monocytes","Macrophages","DCs" ,"Neutrophils","Endothelial","Fibroblasts",
                   "bladder", "breast", "colon", "liver", "lung", "ovary", "pancreas", "thyroid")
    deconv_cell_types = c("Bcells","Naive_CD4_Tcells", "Memory_and_Helper_CD4", "Tregs","Naive_CD8_Tcells", "Non_Naive_CD8_Tcells", "NK","Monocytes","Macrophages","DCs" ,"Neutrophils","Endothelial","Fibroblasts") 
    markers_data$cell_type <- plyr::revalue(markers_data$cell_type, c("Non_Naive_CD4_Tcells" = "Memory_and_Helper_CD4"))
  } else{
    cell_types = c("Bcells","CD4_Tcells","CD8_Tcells","NK","Monocytes","Macrophages","DCs" ,"Neutrophils","Endothelial","Fibroblasts",
                   "bladder", "breast", "colon", "liver", "lung", "ovary", "pancreas", "thyroid")
    deconv_cell_types = c("Bcells","CD4_Tcells","CD8_Tcells", "NK","Monocytes","Macrophages","DCs" ,"Neutrophils","Endothelial","Fibroblasts") 
  } 
  
  markers_data <- markers_data[which(markers_data$cell_type %in% deconv_cell_types[which(deconv_cell_types %in% colnames(averaged_profiles))]),] 
  markers_data_subset = markers_data %>% arrange(factor(cell_type, levels = deconv_cell_types[which(deconv_cell_types %in% colnames(averaged_profiles))] ))
  profile_heatmap <- pheatmap(
    as.matrix(averaged_profiles[unique(markers_data_subset$peak_id[which(markers_data_subset$peak_id %in% peaks )]), cell_types[which(cell_types %in% colnames(averaged_profiles))] ]),
    scale = "row",
    show_rownames = F,# color = hcl.colors(50, "YlOrRd",rev = T),
    cluster_rows = F,
    cluster_cols = F,fontsize = 10, main = profile_heatmap_title, heatmap_legend_param = list(title = "")
  )
  
  if(sample_heatmap){

    meta_data <- meta_data %>%dplyr::arrange(factor(Cell_type, levels = cell_types), Study)
    ordered_mat = norm_data[markers_data_subset$peak_id[which(markers_data_subset$peak_id %in% peaks)], meta_data$Sample_id]
    print(ordered_mat[1:5,1:5])
    if(min_max_scaling){
      ordered_mat = (ordered_mat - rowMeans(ordered_mat, na.rm = TRUE)) / (apply(ordered_mat, 1, max, na.rm = TRUE) - apply(ordered_mat, 1, min, na.rm = TRUE))
    } 
    
    studyCol <-  c("#F0A0FF","#0075DC","#993F00","#4C005C","#191919","#005C31","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088")
    names(studyCol) <- unique(meta_data$Study)
    
    cell_typeCol <- descriminative_colors[which(names(descriminative_colors) %in% unique(meta_data$Cell_type))]
    
    ha = HeatmapAnnotation(
      "Cell type" = meta_data$Cell_type, 
      Study = meta_data$Study,
      col = list("Cell type" = cell_typeCol, Study = studyCol)
    )
    rownames(markers_data) = markers_data$peak_id
    ha_row = rowAnnotation("Cell type" = markers_data[rownames(ordered_mat), "cell_type"] , col = list("Cell type" = cell_typeCol))
    heatmap_ref <- Heatmap(ordered_mat, name = "Norm. counts", 
                           top_annotation = ha,
                           split = factor(markers_data[rownames(ordered_mat), "cell_type"], levels = deconv_cell_types), 
                           show_row_names = F, show_column_names = F, cluster_rows = F, cluster_columns = F, row_title_rot = 0, use_raster = T)
    
    return(list(profile_heatmap, heatmap_ref))
  } else{
    return(profile_heatmap)
    
  } 
  
  
}

ref_profile_heatmap <- get_markers_heatmaps(counts_data = subset_samples, norm_data = dataNorm[markers,] , meta_data = subset_metadata, 
                                            peaks = markers, markers_data = markers_df[which(markers_df$remaining_afterFiltering == T),], profile_heatmap_title = "Reference samples", min_max_scaling = T, sample_heatmap = F)
if(withSuptypes){
  pdf(paste0(fig_path,"fig5_sup3A.pdf"), w = 4, h = 6)
  print(ref_profile_heatmap)
  dev.off()
  
}else{
  pdf(paste0(fig_path,"fig2_panel_C.pdf"), w = 4, h = 6)
  print(ref_profile_heatmap)
  dev.off()
  
}


# Keep only ref samples + entex samples -----------------------------------
subset_samples = combined_samples[, which(colnames(combined_samples) %in% combined_metadata$Sample_id[combined_metadata$Usage == "Reference sample"] | combined_metadata$Study == "ENCODE")]
subset_metadata = combined_metadata[colnames(subset_samples), ]

dataNorm <- gc_norm(input_matrix = subset_samples, 
                    grouping_var = subset_metadata$Cell_type, 
                    fa_file = fasta_file, 
                    by_group = F)
dataNorm <- dataNorm[,subset_metadata$Sample_id] 

heatmap_plots <- get_markers_heatmaps(counts_data = subset_samples, norm_data = dataNorm[markers,] , meta_data = subset_metadata, 
                                      peaks = markers, markers_data = markers_df[which(markers_df$remaining_afterFiltering == T),], profile_heatmap_title = "Reference and ENCODE tissues samples", min_max_scaling = T, sample_heatmap = T)

pdf(paste0(fig_path,"fig1_sup1.pdf"))
print(heatmap_plots[[2]])
dev.off()



# Keep validation samples -------------------------------------------------
subset_samples = combined_samples[, which(!colnames(combined_samples) %in% combined_metadata$Sample_id[combined_metadata$Usage == "Reference sample"] & !combined_metadata$Study %in% c("ENCODE"))]
subset_metadata = combined_metadata[colnames(subset_samples), ]

dataNorm <- gc_norm(input_matrix = subset_samples, 
                    grouping_var = subset_metadata$Cell_type, 
                    fa_file = fasta_file, 
                    by_group = F)
dataNorm <- dataNorm[,subset_metadata$Sample_id] 

validation_heatmap <- get_markers_heatmaps(counts_data = subset_samples, norm_data = dataNorm[markers,], meta_data = subset_metadata, 
                                           peaks = markers, markers_data = markers_df[which(markers_df$remaining_afterFiltering == T),], profile_heatmap_title = "Validation in bulk samples", min_max_scaling = T, sample_heatmap = F)

if(withSuptypes){
  pdf(paste0(fig_path, "fig5_sup3B.pdf"), w = 4, h = 6)
  print(validation_heatmap)
  dev.off()
}else{
  pdf(paste0(fig_path, "fig2_panel_D.pdf"), w = 4, h = 6)
  print(validation_heatmap)
  dev.off()
}




# Heatmap for the single-cell data validation -----------------------------
if(!withSuptypes){
  atac_data <- readRDS(human_atlas_counts)
  
  head(atac_data@meta.data)
  atac_data@meta.data$cell_ID = paste0(atac_data@meta.data$dataset,'+', rownames(atac_data@meta.data)) 
  
  meta_data = read.table(human_atlas_metadata, header = T, sep = "\t")
  meta_data = meta_data[which(meta_data$Life.stage == "Adult"),]
  rownames(meta_data) = meta_data$cellID
  summary(atac_data@meta.data$cell_ID %in% meta_data$cellID)
  meta_data = meta_data[atac_data@meta.data$cell_ID,] 
  
  atac_data = atac_data[,which(meta_data$cell.type != "Naive T cell")] 
  meta_data = meta_data[which(meta_data$cell.type != "Naive T cell"),] 
  
  meta_data$cell.type2 = plyr::revalue(meta_data$cell.type,c("Fibroblast (General)" = "Fibroblasts","Fibroblast (Liver Adrenal)" = "Fibroblasts","Cardiac Fibroblasts" = "Fibroblasts","Fibroblast (Epithelial)" = "Fibroblasts",
                                                             "Fibroblast (Peripheral Nerve)" = "Fibroblasts","Fibroblast (Sk Muscle Associated)" = "Fibroblasts","Fibroblast (Gastrointestinal)" = "Fibroblasts",
                                                             "T Lymphocyte 1 (CD8+)" = "CD8_Tcells",
                                                             "Macrophage (General)" = "Macrophages","Macrophage (General,Alveolar)"  = "Macrophages", 
                                                             "Endothelial Cell (General) 1" = "Endothelial", "Endothelial Cell (Myocardial)" = "Endothelial","Lymphatic Endothelial Cell" = "Endothelial",
                                                             "Endothelial Cell (General) 2" = "Endothelial","Endothelial Cell (General) 3" = "Endothelial", "Blood Brain Barrier Endothelial Cell" = "Endothelial",
                                                             "Natural Killer T Cell" = "NK",
                                                             "T lymphocyte 2 (CD4+)" = "CD4_Tcells",
                                                             "Memory B Cell" = "Bcells",
                                                             "Alveolar Capillary Endothelial Cell" = "Endothelial","Endothelial (Exocrine Tissues)"  = "Endothelial"))
  
  table(meta_data$cell.type2)
  
  
  peaks_gr = Signac::StringToGRanges(rownames(atac_data@assays$peaks@data),sep = c("-","-"))
  region.length = peaks_gr@ranges@width
  names(region.length) = rownames(atac_data@assays$peaks@data)
  atac_data@assays$peaks@data <-  atac_data@assays$peaks@data / region.length[rownames(atac_data@assays$peaks@data)]
  
  atac_data@assays$peaks@data <- Matrix::t(Matrix::t(atac_data@assays$peaks@data) * 1e6 / (sparseMatrixStats::colSums2(atac_data@assays$peaks@data)))
  
  profiles = data.frame(row.names = rownames(atac_data@assays$peaks@data ))
  for(i in 1:length(unique(meta_data$cell.type2))) {
    c_type = as.vector(unique(meta_data$cell.type2))[i]
    cell_index = which(meta_data$cell.type2 %in% c_type)
    profiles = cbind(profiles, sparseMatrixStats::rowMeans2(atac_data@assays$peaks@data[, cell_index,drop = F]))
    colnames(profiles)[ncol(profiles)] = c_type
  }
  
  cell_types = c("Bcells","CD4_Tcells","CD8_Tcells","NK","Macrophages","Endothelial","Fibroblasts")
  markers_df_subset = markers_df[which(markers_df$cell_type %in% cell_types & markers_df$peak_id %in% markers), ]
  markers_df_subset = markers_df_subset %>% arrange(factor(cell_type, levels = c("Bcells","CD4_Tcells","CD8_Tcells","NK","Macrophages","Endothelial","Fibroblasts")))

  
  validation_sc_profile_heatmap <- pheatmap(
    as.matrix(profiles[markers_df_subset$peak_id, c("Bcells","CD4_Tcells","CD8_Tcells","NK","Macrophages","Endothelial","Fibroblasts")]),
    scale = "row",
    show_rownames = F,
    cluster_rows = F,
    cluster_cols = F,fontsize = 17, heatmap_legend_param = list(title = ""), main = "Validation in single-cell data"
  )
  
  pdf(paste0(fig_path,"fig2_panel_E.pdf"), w = 4, h = 6)
  print(validation_sc_profile_heatmap)
  dev.off()
}



# Libraries ---------------------------------------------------------------
library(rtracklayer)
library(readxl)
library(dplyr)
library(ggplot2)
library(cowplot)
library(GenomicRanges)
library(SummarizedExperiment)

# Parameters --------------------------------------------------------------
data_path <- as.character(commandArgs(TRUE)[1])  
fig_path <- as.character(commandArgs(TRUE)[2])

dir.create(fig_path)

# Breast dataset -----------------------------------------------------------

breast_data <- readRDS(paste0(data_path, "JFCR_BRCA_ATAC_42_se.rds"))

meta_data <- breast_data@colData
meta_data$cluster <- "CA-A"
meta_data$cluster[which(meta_data$Patient %in% c("P8","P30", "P54", "P39", "P11", "P12", "P3", "P22", "P14", "P26"))]  <- "CA-B"
meta_data$cluster[which(meta_data$Patient %in% c("P13","P28", "P38", "P19", "P65", "P24L"))]  <- "CA-C"

breast_data_counts <- as.matrix(SummarizedExperiment::assay(breast_data, "counts"))
rownames(breast_data_counts) <- paste0(seqnames(breast_data@rowRanges), "-", start(breast_data@rowRanges), "-", end(breast_data@rowRanges))
colnames(breast_data_counts) <- breast_data@colData$Patient

# Get TPM-like data 
peaks_gr <- Signac::StringToGRanges(rownames(breast_data_counts),sep = c("-","-"))
region.length <- peaks_gr@ranges@width
names(region.length) <- rownames(breast_data_counts)
breast_data_counts <- breast_data_counts / region.length[rownames(breast_data_counts)]
breast_data_counts <- t(t(breast_data_counts) * 1e6 / colSums(breast_data_counts))

epic_res <- EPICATAC::EPIC_ATAC(bulk = breast_data_counts, reference = EPICATAC::atacRef_TME, withOtherCells = T, ATAC = T, genome_version = "hg38") 

epic_output <- as.data.frame(epic_res$cellFractions)
epic_output$T_cells <- epic_output$CD4_Tcells + epic_output$CD8_Tcells
epic_output$Myeloid <- epic_output$Macrophages + epic_output$DCs + epic_output$Neutrophils

epic_output$samples <- rownames(epic_output)
head(epic_output)
epic_output$Group <- meta_data[epic_output$samples,"Subtype"] 
epic_output$Cluster <- meta_data[epic_output$samples,"cluster"] 
epic_output$Subtype <- plyr::revalue(epic_output$Group, c("L" = "ER+/HER2-", "TN" = "TNBC"))


# Generate main grid figure -----------------------------------------------

text_size <- 15
legend_text_size <- 14
my_theme <- theme(text = element_text(size = text_size, face = "bold"),
                  axis.text = element_text(size = text_size, face = "bold"),
                  legend.title = element_text(size = legend_text_size, face = "bold"),
                  legend.text = element_text(size = legend_text_size),
                  axis.title = element_text(size = text_size, face = "bold"),
                  strip.text = element_text(size = text_size, face = "bold"),
                  axis.text.x = element_text(angle=45, hjust=1, size = text_size))

plot_df <- tidyr:: gather(epic_output, cell_type, proportion, Macrophages:Myeloid, factor_key=TRUE)
plot_df <- plot_df[plot_df$cell_type %in% c("Endothelial", "Fibroblasts", "Myeloid", "Bcells", "T_cells", "NK"), ]
plot_df$cell_type <- factor(plot_df$cell_type, levels=c("Endothelial", "Fibroblasts", "Myeloid", "Bcells", "T_cells", "NK")) 

my_comparisons <- list( c("ER+/HER2-","TNBC") )
subtype_boxplots <- ggpubr::ggboxplot(plot_df[which(plot_df$cell_type %in% c(  "Bcells", "T_cells", "NK",  "Endothelial", "Fibroblasts", "Myeloid")),], x = "Subtype", y = "proportion", 
                                      fill = "Subtype",  facet.by = "cell_type", nrow = 1)+ 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.14))) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.format",paired = F, vjust =  -0.01, size = 4) + labs(x = "", y = "Proportions")+
  theme_bw(base_size = 12) + theme(legend.position = "none")+ theme_bw() + my_theme 


my_comparisons <- list( c("CA-A", "CA-B"),c("CA-B", "CA-C"),c("CA-A", "CA-C") )
cluster_boxplots <- ggpubr::ggboxplot(plot_df[which(plot_df$cell_type %in% c(  "Bcells", "T_cells", "NK",  "Endothelial", "Fibroblasts", "Myeloid")),], x = "Cluster", y = "proportion", 
                                      fill = "Cluster",  facet.by = "cell_type", nrow = 1)+ 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.14))) +
  ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.format",paired = F,  vjust =  -0.01, size = 4) + labs(x = "", y = "Proportions")+
  theme_bw(base_size = 12) + theme(legend.position = "none") + theme_bw() + my_theme


# Load CREs in hg19 -------------------------------------------------------
CREs <- as.data.frame(read_excel(paste0(data_path, "41523_2022_438_MOESM3_ESM.xlsx"), sheet = 2, skip = 1))

# Convert coordiantes to hg38 ---------------------------------------------

peaks_gr <- makeGRangesFromDataFrame(CREs, keep.extra.columns = T)
peaks_gr$id <- paste0(CREs$seqnames, ":", CREs$start, "-", CREs$end)
ch <- EPICATAC::liftover_chains$hg19_to_hg38_chain
GenomeInfoDb::seqlevelsStyle(ch) <- "UCSC"
hg38_regions <- unlist(rtracklayer::liftOver(peaks_gr, ch))

# Match bulk peaks with the CREs ------------------------------------------
cellTypeDAR <- lapply(unique(hg38_regions$cellType), function(x){
  res <- hg38_regions[hg38_regions$cellType == x]
  res <- subsetByOverlaps(res, hg38_regions[hg38_regions$cellType == "Epithelial"], invert = T)
  return(res)
})
names(cellTypeDAR) <- unique(hg38_regions$cellType)
cellTypeDAR_total <- GRangesList(cellTypeDAR) %>% unlist %>% unique

overlapsCRE <- lapply(cellTypeDAR, function(x) subsetByOverlaps(Signac::StringToGRanges(rownames(breast_data_counts)), x))


# Get mean signal of matched CREs -----------------------------------------

df <- lapply(cellTypeDAR[2:7], function(x){
  overlaps <- findOverlaps(rowRanges(breast_data), x)
  out <- colMeans(breast_data_counts[queryHits(overlaps),])
  return(out)
})
df <- as.data.frame(do.call(cbind, df))
df <- cbind(df, Subtype = as.character(breast_data@colData$Subtype))

epic_atac_res <- data.frame(sample = epic_output$samples,
                            Endothelial = epic_output$Endothelial,
                            Fibroblasts = epic_output$Fibroblasts,
                            T_cells = epic_output$CD4_Tcells + epic_output$CD8_Tcells,
                            Bcells = epic_output$Bcells,
                            Macrophages = epic_output$Macrophages,
                            NK = epic_output$NK,
                            DCs = epic_output$DCs,
                            Neutrophils = epic_output$Neutrophils,
                            Myeloid = epic_output$Macrophages + epic_output$DCs + epic_output$Neutrophils,
                            subtype = epic_output$Group)
epic_atac_res <- epic_atac_res %>% 
  tidyr::pivot_longer(Endothelial:Myeloid, names_to = 'Cell_types', values_to = 'EPIC_ATAC_prop')

original_res <- data.frame(sample = rownames(df),
                           Endothelial = df$Endothelial,
                           Fibroblasts = df$Fibroblast,
                           T_cells = df$T,
                           Bcells = df$B,
                           Macrophages = df$Myeloid,
                           Myeloid = df$Myeloid,
                           DCs = df$Myeloid,
                           Neutrophils = df$Myeloid,
                           NK = 0,
                           subtype = df$Subtype)
original_res <- original_res %>% 
  tidyr::pivot_longer(Endothelial:NK, names_to = 'Cell_types', values_to = 'Paper_values')

combined_df <- inner_join(epic_atac_res, original_res, 
                          by = c("sample", "subtype", "Cell_types"))

combined_df$Cell_types <- factor(combined_df$Cell_types, levels=c("Endothelial", "Fibroblasts", "Macrophages","Myeloid", "Bcells", "T_cells", "NK", "DCs", "Neutrophils"))

plot_subset <- combined_df[combined_df$Cell_types %in% c("Endothelial", "Fibroblasts","Myeloid", "Bcells", "T_cells", "NK"),]

paper_comp <- ggplot(data = plot_subset, aes(x = Paper_values, y = EPIC_ATAC_prop)) + geom_point() +
  ggpubr::stat_cor(method = "pearson",cor.coef.name = c("R"), size = 4) + facet_wrap(~Cell_types, nrow = 1,scales = "free_x") + 
  geom_smooth(method = "lm", se = T, color = "darkgray",) + 
  ylab("EPIC-ATAC (%)") + xlab("Original paper (Mean ATAC signal at CREs)") + theme_bw() + my_theme


boxplots_panel <- cowplot::plot_grid(paper_comp, NULL, subtype_boxplots + theme(legend.position = 'hidden'), 
                                     NULL, 
                                     cluster_boxplots + theme(legend.position = 'hidden'), 
                                     label_size = 20, ncol = 1, nrow = 5, rel_heights = c(0.9, 0.01,  1, 0.01, 0.9))

boxplots_panel 

pdf(paste0(fig_path, "Figure6.pdf"), w = 15, h = 10)
boxplots_panel
dev.off()


# Only Myeloid cells ------------------------------------------------------

plot_subset <- combined_df[combined_df$Cell_types %in% c("Macrophages", "DCs", "Neutrophils"),]

text_size = 15
legend_text_size = 14
my_theme <- theme(text = element_text(size = text_size, face = "bold"),
                  axis.text = element_text(size = text_size, face = "bold"),
                  legend.title = element_text(size = legend_text_size, face = "bold"),
                  legend.text = element_text(size = legend_text_size),
                  axis.title = element_text(size = text_size, face = "bold"),
                  strip.text = element_text(size = text_size, face = "bold"),
                  axis.text.x = element_text(size = text_size))

paper_comp <- ggplot(data = plot_subset, aes(x = Paper_values, y = EPIC_ATAC_prop)) + geom_point() +
  ggpubr::stat_cor(method = "pearson",cor.coef.name = c("R"), size = 4) + facet_wrap(~Cell_types, nrow = 1,scales = "free_x") + 
  geom_smooth(method = "lm", se = T, color = "darkgray",) + 
  ylab("EPIC-ATAC (%)") + xlab("Original paper (Mean ATAC signal at myeloid CREs)") + theme_bw() + my_theme

plot_df <- tidyr:: gather(epic_output, cell_type, proportion, Macrophages:Myeloid, factor_key=TRUE)
plot_df <- plot_df[plot_df$cell_type %in% c("Macrophages", "DCs", "Neutrophils"), ]
plot_df$cell_type <- factor(plot_df$cell_type, levels=c("Macrophages", "DCs", "Neutrophils")) 

my_comparisons <- list( c("ER+/HER2-","TNBC") )
subtype_boxplots <- ggpubr::ggboxplot(plot_df[which(plot_df$cell_type %in% c("Macrophages", "DCs", "Neutrophils")),], x = "Subtype", y = "proportion", 
                                      fill = "Subtype",  facet.by = "cell_type", nrow = 1)+ 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.14))) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.format",paired = F, vjust =  -0.01, size = 4) + labs(x = "", y = "Proportions")+
  theme_bw(base_size = 12) + theme(legend.position = "none")+ theme_bw() + my_theme 


my_comparisons <- list( c("CA-A", "CA-B"),c("CA-B", "CA-C"),c("CA-A", "CA-C") )
cluster_boxplots <- ggpubr::ggboxplot(plot_df[which(plot_df$cell_type %in% c("Macrophages", "DCs", "Neutrophils")),], x = "Cluster", y = "proportion", 
                                      fill = "Cluster",  facet.by = "cell_type", nrow = 1)+ 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.14))) +
  ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.format",paired = F,  vjust =  -0.01, size = 4) + labs(x = "", y = "Proportions")+
  theme_bw(base_size = 12) + theme(legend.position = "none") + theme_bw() + my_theme



boxplots_panel <- cowplot::plot_grid(paper_comp, NULL,
                                     subtype_boxplots + theme(legend.position = 'hidden'), 
                                     NULL, 
                                     cluster_boxplots + theme(legend.position = 'hidden'), 
                                     label_size = 20, ncol = 1, nrow = 5, rel_heights = c(1, 0.01,  1, 0.01, 1))

boxplots_panel 

pdf(paste0(fig_path, "Figure6_sup_Myeloid.pdf"), w = 10, h = 10)
boxplots_panel
dev.off()

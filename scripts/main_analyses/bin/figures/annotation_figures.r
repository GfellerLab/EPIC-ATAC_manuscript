library(ComplexHeatmap)
library(dplyr)
library(ggplot2)

# Parameters --------------------------------------------------------------
fig_path = as.character(commandArgs(TRUE)[1])
annotation_path = as.character(commandArgs(TRUE)[2])
dir.create(fig_path)

descriminative_colors = c("Bcells"="darkgreen","CD4_Tcells"="#3288bdff","CD8_Tcells"="#4daf4aff","DCs"="#984ea3ff","Dendritic"="#984ea3ff",
                          "Macrophages"="#ff7f00ff","Monocytes"="#ffd92fff", "NK"="#70c4b4ff",
                          "Neutrophils"="black","Endothelial"="#a0451fff","Fibroblasts"="gray")


# Load annotations --------------------------------------------------------

anno_output_df_PBMC = read.table(paste0(annotation_path, "PBMC_noSubtypes/GO_enrichment.txt"),
                            header = T, sep = "\t", quote = "", comment.char = "")

anno_output_df_Tumor = read.table(paste0(annotation_path, "TME_noSubtypes/GO_enrichment.txt"),
                                 header = T, sep = "\t", quote = "", comment.char = "")

anno_output_df_PBMC$profile = "PBMC"
anno_output_df_Tumor$profile = "TME"

anno_all = rbind(anno_output_df_PBMC, anno_output_df_Tumor)
anno_all = unique(anno_all[, c("Description", "cell_type", "profile")] )
openxlsx::write.xlsx(list(PBMC = anno_output_df_PBMC, Tumor = anno_output_df_Tumor, combined = anno_all), file = paste0(fig_path, "GO_enrichment.xlsx"))


# Plot heatmap ------------------------------------------------------------
anno_all = rbind(anno_output_df_PBMC, anno_output_df_Tumor)

manual_selection <- list("Bcells" = c("B cell receptor signaling pathway", 
                                      "regulation of cytosolic calcium ion concentration", 
                                      "B cell activation",
                                      "regulation of B cell proliferation", 
                                      "B cell mediated immunity", 
                                      "regulation of interleukin-10 production",
                                      "interleukin-4 production", 
                                      "negative regulation of kinase activity", 
                                      "nitric-oxide synthase biosynthetic process",
                                      "positive regulation of B cell mediated immunity"),
                         "CD4_Tcells" = c("regulation of T cell activation", 
                                          "T cell costimulation", "regulation of T cell proliferation", "T-helper 17 type immune response", 
                                          "interleukin-2 biosynthetic process", 
                                          "cellular response to interleukin-4", 
                                          "T cell receptor signaling pathway", 
                                          "regulation of cell killing"),
                         "CD8_Tcells" = c("T cell activation", "adaptive immune response", 
                                          "regulation of cell killing", "cell killing", "T cell mediated cytotoxicity", 
                                          "MyD88-independent toll-like receptor signaling pathway",
                                          "regulation of T cell receptor signaling pathway"),
                         "Endothelial" = c("blood vessel development", "regulation of angiogenesis", 
                                           "angiogenesis involved in wound healing", 
                                           "morphogenesis of an endothelium"),
                         "Fibroblasts" = c("regulation of necroptotic process", "beta-catenin destruction complex disassembly"),
                         "Macrophages" = c('tumor necrosis factor biosynthetic process', 
                                           "regulation of intrinsic apoptotic signaling pathway", 
                                           "cytokine production", "regulation of interleukin-12 production",
                                           "ERK1 and ERK2 cascade", 
                                           "macrophage chemotaxis", "macrophage migration"
                                           ),
                         "Neutrophils" = c("neutrophil degranulation", 
                                           "neutrophil activation",
                                           "neutrophil mediated immunity", 
                                           "glucan metabolic process"
                                           ),
                         "NK" = c("innate immune response", "response to tumor cell", 
                                  "natural killer cell mediated immunity","response to virus", 
                                  "defense response",
                                  "Fc receptor signaling pathway"
                                  )
)


go_res_subset <- anno_all[which(anno_all$Description %in% unlist(manual_selection)),]
go_res_subset <- go_res_subset[order(go_res_subset$cell_type),] 
go_res_subset$Significance <- -log10(go_res_subset$FDR)
go_res_subset$Description <- factor(go_res_subset$Description, levels=unique((go_res_subset$Description)[order(go_res_subset$cell_type)]))
go_res_subset$cell_type2 <- paste0(go_res_subset$cell_type, "-", go_res_subset$profile)
go_res_subset_wide <- tidyr::pivot_wider(go_res_subset[, c('Description', "Significance", "cell_type2")] , names_from = cell_type2, values_from = Significance)
plot_matrix <- t(as.matrix(go_res_subset_wide[,2:ncol(go_res_subset_wide)]))
plot_matrix[is.na(plot_matrix)] <- 0
colnames(plot_matrix) <- go_res_subset_wide$Description

plot_matrix <- plot_matrix[, unique(unlist(manual_selection))]
# identify common pathways
manual_selection[["common_pathways"]] <- colnames(plot_matrix)[which(colSums(plot_matrix != 0) > 1)]

cn = colnames(plot_matrix)
col_fun = circlize::colorRamp2(c(0, 6), c("white", "firebrick3"))
pdf(paste0(fig_path, "fig2_panel_G.pdf"), w = 15, h = 6)
ht <- Heatmap(plot_matrix, cluster_rows = F, cluster_columns = F, show_column_dend = F, name = "Significance", col = col_fun, # rev(brewer.pal(n = 7, name = "RdYlBu")),
        split = as.factor(sapply(rownames(plot_matrix), function(i) unlist(strsplit(i, "-"))[[1]])), row_title_rot = 0,
        column_names_rot = 45, row_names_side = "left", border = T, rect_gp = gpar(col = "black", lwd = 1),
        column_names_gp = grid::gpar(fontsize = 12,
                                     fontface =  sapply(colnames(plot_matrix), function(i){
                                       mapping <- as.integer(mapply(`%in%`, i, manual_selection))
                                       if (sum(mapping) > 1) {
                                         return("bold")
                                       } else{
                                         return("plain")
                                       }
                                       }),
                                     col =  sapply(colnames(plot_matrix), function(i){
                                       mapping <- as.integer(mapply(`%in%`, i, manual_selection))
                                       descriminative_colors[names(manual_selection)[which(mapping == 1)][1]]
                                     })),
        row_names_gp = grid::gpar(fontsize = 15, fontface = "bold", col = descriminative_colors[sapply(rownames(plot_matrix), function(i) unlist(strsplit(i, "-"))[[1]])] ),  row_title=NULL
)
draw(ht, padding = unit(c(10, 20, 2, 2), "mm"))
dev.off()



# Plot the piecharts ------------------------------------------------------

# marker peaks annotations (PBMC filtering)
chipSeeker_annot_path = paste0(annotation_path, "PBMC_noSubtypes/chipSeeker_output.txt")
chipSeeker_annotations_PBMC <- read.table(chipSeeker_annot_path, header = T, sep = "\t", quote = "")
rownames(chipSeeker_annotations_PBMC) <- chipSeeker_annotations_PBMC$peak_id

CBPs_enrichment_res_PBMC <- read.table(paste0(annotation_path, "PBMC_noSubtypes/combined_enrichment_res.txt"),
                                       header = T, sep = "\t")

# marker peaks annotations (TME filtering)
chipSeeker_annot_path = paste0(annotation_path, "TME_noSubtypes/chipSeeker_output.txt")
chipSeeker_annotations_TME <- read.table(chipSeeker_annot_path, header = T, sep = "\t", quote = "")
rownames(chipSeeker_annotations_TME) <- chipSeeker_annotations_TME$peak_id

CBPs_enrichment_res_TME <- read.table(paste0(annotation_path, "TME_noSubtypes/combined_enrichment_res.txt"),
                                      header = T, sep = "\t")

plot_piecharts <- function(peaks_annotations, combined_enrichment_res){

  bins <-  c(-3000000,-100000,-10000,-5000,-1000,0,1000,5000,10000,100000,3000000)
  peaks_annotations$distanceToTSS_categ <- cut(peaks_annotations$distanceToTSS, breaks = bins, right = F, include.lowest = T)
  table(peaks_annotations$distanceToTSS_categ)
  
  peaks_annotations$distanceToTSS_categ <- plyr::revalue(peaks_annotations$distanceToTSS_categ, c("[-3e+06,-1e+05)" = "<-100k",
                                                                                                  "[-1e+05,-1e+04)" = "-100k",
                                                                                                  "[-1e+04,-5e+03)" = "-10k",
                                                                                                  "[-5e+03,-1e+03)" = "-5k",
                                                                                                  "[-1e+03,0)" = "-1k",
                                                                                                  "[0,1e+03)" = "1k",
                                                                                                  "[1e+03,5e+03)" = "5k",
                                                                                                  "[5e+03,1e+04)" = "10k",
                                                                                                  "[1e+04,1e+05)" = "100k",
                                                                                                  "[1e+05,3e+06]" = ">100k"))
  
  distance_df <- peaks_annotations[,c("distanceToTSS_categ"), drop = F] 
  distance_df %>% 
    group_by(distanceToTSS_categ) %>% 
    dplyr::summarise(counts = n()) -> distance_df_counts
  distance_df_counts = as.data.frame(distance_df_counts)
  distance_df_counts$distanceToTSS_categ = factor(distance_df_counts$distanceToTSS_categ, levels = c("1k","5k","10k","100k",">100k","<-100k","-100k","-10k" ,"-5k","-1k"))
  
  ordered_df_distances =vector()
  for(i in c("1k","5k","10k","100k",">100k","<-100k","-100k","-10k" ,"-5k","-1k")){
    ordered_df_distances = rbind(ordered_df_distances, distance_df_counts[which(distance_df_counts$distanceToTSS_categ == i),] )
  } 
  
  
  cols <- c("Distal Intergenic" = "#33a02cff", "Downstream" = "#ff7f00ff", "Intron" = "#fdbf6fff", 
            "Exon" = "#d95f02ff", "5' UTR" = "#f21c39", "3' UTR" = "#fb9a99ff", "Promoter" = "#a6cee3ff",
            "<-100k" = "#dfbfe6ff",
            "-100k" = "#88dbd3ff",
            "-10k"= "#9bc48dff",
            "-5k" = "#dcceb4ff",
            "-1k" = "#3182bdff",
            "1k"= "#31688eff",
            "5k" = "#c7a76cff",
            "10k" = "#86b875ff",
            "100k" = "#39beb1ff",
            ">100k" = "#cd99d8ff")
  
  ordered_df_distances <- ordered_df_distances %>% 
    mutate(x = 1.5)
  ordered_df_distances$distanceToTSS_categ <- factor(ordered_df_distances$distanceToTSS_categ, levels = c("-1k","-5k","-10k" ,"-100k","<-100k", ">100k","100k","10k","5k","1k"))
  
  TSS_piechart <-ggplot(ordered_df_distances, aes(x = 1.5, y = counts, fill = distanceToTSS_categ)) +
    geom_col(color = "black") +
    geom_text(aes(label = counts),
              position = position_stack(vjust = 0.5), size = 7) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = cols, breaks=c("1k", "5k", "10k", "100k", ">100k", "-1k", "-5k", "-10k", "-100k", "<-100k" )) +
    xlim(c(0.2, 1.5 + 0.5)) +
    guides(fill = guide_legend(ncol=2, title="Distance to TSS")) +
    theme_void() + 
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position = 'left',
          legend.direction = "vertical",
          legend.title = element_text( size = 25, face ="bold"),
          legend.text = element_text( size = 25))
  TSS_piechart
  
  annot_df <- peaks_annotations[,c("annotation_type"), drop = F] 
  annot_df %>% 
    group_by(annotation_type) %>% 
    dplyr::summarise(counts = n()) -> annot_df_counts
  annot_df_counts = as.data.frame(annot_df_counts)
  
  annot_df_counts <- annot_df_counts %>% 
    mutate(x = 1.5)
  
  annotation_piechart <-ggplot(annot_df_counts, aes(x = 1.5, y = counts, fill = annotation_type)) +
    geom_col(color = "black") +
    geom_text(aes(label = counts),
              position = position_stack(vjust = 0.5), size = 7) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = cols) +
    xlim(c(0.2, 1.5 + 0.5)) +
    guides(fill = guide_legend(ncol=1, title="ChipSeeker annotation")) +
    theme_void() + 
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position = 'left',
          legend.direction = "vertical",
          legend.title = element_text( size = 25, face ="bold"),
          legend.text = element_text( size = 25))

  
  return(cowplot::plot_grid(
    TSS_piechart, NULL,
    annotation_piechart,
    labels = c("", "", ""), ncol = 1,  label_size = 20, rel_heights = c(0.9, 0.05, 0.9)
  ))
}

PBMCpiecharts <- plot_piecharts(peaks_annotations = chipSeeker_annotations_PBMC, combined_enrichment_res = CBPs_enrichment_res_PBMC)
TMEpiecharts <- plot_piecharts(peaks_annotations = chipSeeker_annotations_TME, combined_enrichment_res = CBPs_enrichment_res_TME)

pdf(paste0(fig_path, "fig2_panel_F.pdf"))
PBMCpiecharts
TMEpiecharts
dev.off()


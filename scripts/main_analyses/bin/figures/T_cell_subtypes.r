
# Libraries ---------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(Metrics)
# Parameters --------------------------------------------------------------

descriminative_colors = c("CD4 Tcells"="#3288bdff","CD8 Tcells"="#4daf4aff",
                          "Naive CD4"="#ff7f00ff","Memory + Helper CD4"="#660033","Tregs"="#f24d90ff", 
                          "Naive CD8"="#70c4b4ff","Non Naive CD8"="#ffd92fff")

benchmarking_output_path = as.character(commandArgs(TRUE)[1])
fig_path = as.character(commandArgs(TRUE)[2])
dir.create(fig_path)


# Source functions --------------------------------------------------------
deconvolution_scatter_plot <- function(deconv_data, method_name, bulk,xlim = NULL, ylim = NULL, ylab = "Estimated proportions", xlab = "True proportions", absolute = T, annot_text_size = 3, point_size = 15){
  plot_df <- deconv_data %>% dplyr::filter(method %in% method_name) %>% dplyr::filter(bulk_name %in% bulk) 
  
  p <- plot_df %>% 
    ggplot(aes(x = true_fraction, y = estimate)) +
    geom_abline(slope = 1, intercept = 0, color="grey", linewidth = 2) + 
    geom_point(aes(color = cell_type, shape = sample), size = point_size) + ylim(0, ylim) + xlim(0, xlim)+
    scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% plot_df$cell_type)]) + 
    scale_shape_manual(values = c(16,17,18,8,7,15,9,10,12,13,14,1,2))+
    facet_grid(~ bulk_name, scales = "free", switch = "y")+
    cowplot::background_grid(major="xy") +
    guides(col=guide_legend("Cell types"),
           shape=guide_legend("Sample")) + ylab(ylab) + xlab(xlab) + theme_bw() + 
    theme(legend.position="right") + 
    guides(color = guide_legend(ncol = 1,title = NULL),shape = "none")
  
  text_to_add = paste0("RMSE=", round(Metrics::rmse(plot_df$estimate, plot_df$true_fraction),3),#With cancer cells:   
                       "\n COR=", round(cor.test(plot_df$estimate, plot_df$true_fraction)$estimate,2))
  annotations <- data.frame(xpos = c(-Inf), ypos =  c(Inf),
                            annotateText = c(text_to_add),
                            hjustvar = c(-0.1),
                            vjustvar = c(1.5))
  p <- p + geom_text(data = annotations, aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText), size = annot_text_size)
  return(p)
} 

#################################################################
##         Predictions of T cell subtypes proportions          ##
#################################################################

axis_text_size = 16
text_size = 19
legend_text_size = 14
strip_text_size = 14

my_theme_scatter <- theme(
  axis.text = element_text(size = axis_text_size, face = "bold"),
  axis.title = element_text(size = text_size, face = "bold"),
  strip.text = element_text(size = strip_text_size, face = "bold", colour = "black"),
  strip.background = element_rect(fill = "grey"),
  axis.title.y = element_text(vjust = +1.5),
  axis.title.x = element_text(vjust = -0.5),
  legend.title = element_text(size = text_size, face = "bold"),
  legend.text = element_text(size = legend_text_size),
  legend.key.size = unit(1.5, 'lines')
) 

my_theme_barplot <- theme(
  axis.text = element_text(size = axis_text_size, face = "bold"),
  axis.title = element_text(size = text_size, face = "bold"),
  strip.text = element_text(size = strip_text_size, face = "bold", colour = "black"),
  strip.background = element_rect(fill = "grey"),
  axis.title.y = element_text(vjust = +1.5),
  axis.title.x = element_text(vjust = -0.5),
  legend.title = element_text(size = text_size, face = "bold"),
  legend.text = element_text(size = legend_text_size),
  legend.key.size = unit(1.5, 'lines')
) 

for(method in c("EPIC-ATAC", "quanTIseq", "ABIS", "CIBERSORTx_ourMarkers", "DeconPeaker_ourMarkers", "CIBERSORTx-Custom", "DeconPeaker-Custom")){
  if(method %in% c("EPIC-ATAC", "quanTIseq")){
    summary_predictions <- read.table(paste0(benchmarking_output_path,"summary_atac_predictions.txt"), header = T, sep = "\t")
  } else{
    summary_predictions <- read.table(paste0(benchmarking_output_path,"summary_atac_predictions-renorm.txt"), header = T, sep = "\t")
  } 
  
  summary_predictions_subset <- summary_predictions[which(summary_predictions$method == method),] 
  summary_predictions_subset <- summary_predictions_subset[which(summary_predictions_subset$cell_type %in% c("Tregs", "Memory + Helper CD4", "Naive CD4", "Non Naive CD8", "Naive CD8")), ]
  summary_predictions_subset$bulk_name[which(summary_predictions_subset$bulk_name %in% c( "PBMC pseudobulk", "PBMC experiment"))]  <- "PBMC (exp.+pseudo.)"
  
  # Per cell type correlations and RMSE values  
  cor_rmse_summary <- data.frame(cell_type = vector(), bulk = vector(), Correlation = vector(), RMSE = vector())
  for(bulk in c("Basal cell carcinoma", "PBMC (exp.+pseudo.)")){
    for(cell_type in unique(summary_predictions_subset$cell_type)){
      estimates = summary_predictions_subset[which(summary_predictions_subset$bulk_name == bulk & summary_predictions_subset$cell_type == cell_type), "estimate"]
      true_prop = summary_predictions_subset[which(summary_predictions_subset$bulk_name == bulk & summary_predictions_subset$cell_type == cell_type), "true_fraction"]
      cor_rmse_summary <- rbind(cor_rmse_summary, data.frame(cell_type = cell_type, bulk = bulk, Correlation = cor.test(estimates, true_prop)$estimate, RMSE = rmse(estimates, true_prop)))
    }
  }
  
  for(bulk in c("Basal cell carcinoma", "PBMC (exp.+pseudo.)")){
    plot_name = ifelse(grepl("PBMC", bulk), "PBMC", bulk)
    assign(x = paste0(gsub(" ", "_", plot_name ), "_Tcells_subtypes"), 
           value = deconvolution_scatter_plot(deconv_data = summary_predictions_subset,
                                              method_name = method,
                                              bulk = bulk,
                                              xlim = 1 , ylim = 1, point_size = 4, annot_text_size = 5, xlab = "True proportions", ylab = "Predictions"))
  }   
  
  
  if(method %in% c("EPIC-ATAC", "quanTIseq")){
    summary_predictions <- read.table(paste0(gsub("benchmarking_withSubtypes", "benchmarking", benchmarking_output_path), "summary_atac_predictions.txt"), header = T, sep = "\t")
  } else{
    summary_predictions <- read.table(paste0(gsub("benchmarking_withSubtypes", "benchmarking", benchmarking_output_path), "summary_atac_predictions-renorm.txt"), header = T, sep = "\t")
  } 
  summary_predictions_subset <- summary_predictions[which(summary_predictions$method == method),]
  summary_predictions_subset <- summary_predictions_subset[which(summary_predictions_subset$cell_type %in% c("CD4_Tcells", "CD8_Tcells")), ]
  summary_predictions_subset$bulk_name[which(summary_predictions_subset$bulk_name %in% c( "PBMC pseudobulk", "PBMC experiment"))]  <- "PBMC (exp.+pseudo.)"
  
  summary_predictions_subset$cell_type <- plyr::revalue(summary_predictions_subset$cell_type, c("CD4_Tcells" = "CD4 Tcells", "CD8_Tcells" = "CD8 Tcells"))
  
  # Per cell type correlations and RMSE values  
  for(bulk in c("Basal cell carcinoma", "PBMC (exp.+pseudo.)")){
    for(cell_type in unique(summary_predictions_subset$cell_type)){
      estimates = summary_predictions_subset[which(summary_predictions_subset$bulk_name == bulk & summary_predictions_subset$cell_type == cell_type), "estimate"]
      true_prop = summary_predictions_subset[which(summary_predictions_subset$bulk_name == bulk & summary_predictions_subset$cell_type == cell_type), "true_fraction"]
      cor_rmse_summary <- rbind(cor_rmse_summary, data.frame(cell_type = cell_type, bulk = bulk, Correlation = cor.test(estimates, true_prop)$estimate, RMSE = rmse(estimates, true_prop)))
    }
  }
  
  for(bulk in c("Basal cell carcinoma", "PBMC (exp.+pseudo.)")){
    plot_name = ifelse(grepl("PBMC", bulk), "PBMC", bulk)
    assign(x = paste0(gsub(" ", "_", plot_name ), "_Tcells"), 
           value = deconvolution_scatter_plot(deconv_data = summary_predictions_subset,
                                              method_name = method,
                                              bulk = bulk,
                                              xlim = 1 , ylim = 1,  point_size = 4, annot_text_size = 5, xlab = "True proportions", ylab = "Predictions"))
  }   

  
  cor_rmse_summary$bulk <- factor(cor_rmse_summary$bulk, levels = c("PBMC (exp.+pseudo.)", "Basal cell carcinoma"))
  
  
  cor_rmse_summary$cell_categ <- ifelse(cor_rmse_summary$cell_type %in% c("CD4 Tcells", "CD8 Tcells"), "Major types", "Subtypes")
  cor_rmse_summary$cell_type <- factor(cor_rmse_summary$cell_type, levels = c("Tregs", "Memory + Helper CD4", "Naive CD4", "Non Naive CD8", "Naive CD8", "CD4 Tcells", "CD8 Tcells"))
  cor_plot <- cor_rmse_summary %>% ggplot(aes(x = Correlation, y = cell_type, fill = cell_type )) +
    geom_bar(stat = "identity", position = "dodge", show.legend = F, width = 0.7) + labs( y = "") +
    facet_wrap(~ bulk, nrow = 2) +
    scale_fill_manual(values=descriminative_colors) +
    theme_linedraw() + my_theme_barplot
  
  
  scatter_PBMC <- cowplot::plot_grid(
    PBMC_Tcells + theme(legend.position = "none") + my_theme_scatter, NULL, 
    PBMC_Tcells_subtypes + theme(legend.position = "none") + my_theme_scatter,
    ncol = 3, nrow = 1, rel_widths = c(0.9, 0.05, 0.9)
  )
  
  scatter_BCC <- cowplot::plot_grid(
    Basal_cell_carcinoma_Tcells + theme(legend.position = "none") + my_theme_scatter, NULL, 
    Basal_cell_carcinoma_Tcells_subtypes + theme(legend.position = "none") + my_theme_scatter,
    ncol = 3, nrow = 1, rel_widths = c(0.9, 0.05, 0.9)
  )
  
  scatter_combined <- cowplot::plot_grid(
    scatter_PBMC, NULL, scatter_BCC,
    ncol = 1, nrow = 3, rel_heights = c(0.9, 0.05, 0.9)
  )
  
  legends <- cowplot::plot_grid(
    cowplot::get_legend(PBMC_Tcells + my_theme_scatter), NULL, cowplot::get_legend(PBMC_Tcells_subtypes + my_theme_scatter),
    ncol = 1, nrow = 3, rel_heights = c(0.9, 0.05, 0.9)
  )
  
  assign(x= paste0(gsub("-", "_", method), "_subtype_fig"), value = cowplot::plot_grid(
    scatter_combined, legends, NULL, cor_plot, ncol = 4,  
    rel_widths =  c(1, 0.4, 0.01, 0.7), labels = c("", "", "")
  ))
} 

pdf(paste0(fig_path, "Figure5_subtypes_benchmark.pdf"), w = 15, h = 7)
EPIC_ATAC_subtype_fig
dev.off()


pdf(paste0(fig_path, "Figure5_subtypes_benchmark_sup.pdf"), w = 15, h = 7)
for(method in c("quanTIseq", "ABIS", "CIBERSORTx_ourMarkers", "DeconPeaker_ourMarkers", "CIBERSORTx-Custom", "DeconPeaker-Custom")){
  print(cowplot::plot_grid(
    NULL, 
    get(paste0(gsub("-", "_", method), "_subtype_fig")),
    ncol = 1, nrow = 2, rel_heights = c(0.1, 0.9),
    labels = c(method), label_size = 30
  ))

} 
dev.off()

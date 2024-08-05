
# Libraries ---------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(Metrics)

# Functions --------------------------------------------------------
deconvolution_scatter_plot <- function(deconv_data, method_name, bulk,xlim = NULL, ylim = NULL, ylab = "Estimated proportions", xlab = "True proportions", absolute = T, annot_text_size = 3, point_size = 15){
  plot_df <- deconv_data %>% dplyr::filter(method %in% method_name) %>% dplyr::filter(bulk_name %in% bulk) 
  
  p <- plot_df %>% 
    ggplot(aes(x = true_fraction, y = estimate)) +
    geom_abline(slope = 1, intercept = 0, color="grey", linewidth = 2) + 
    geom_point(aes(color = cell_type, shape = sample), size = point_size) + ylim(0, ylim) + xlim(0, xlim)+
    scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% plot_df$cell_type)]) + 
    scale_shape_manual(values = c(16,17,18,8,7,15,9,10,12,13,14,1,2,0, 3, 4, 5, 6, 11, 19, 20, 21, 22, 23, 24))+
    facet_grid(~ bulk_name, scales = "free", switch = "y")+
    cowplot::background_grid(major="xy") +
    guides(col=guide_legend("Cell types"),
           shape=guide_legend("Sample")) + ylab(ylab) + xlab(xlab) + theme_bw() + 
    theme(legend.position="right") + 
    guides(color = guide_legend(ncol = 1,title = NULL),shape = "none")
  
  text_to_add = paste0("RMSE=", round(Metrics::rmse(plot_df$estimate, plot_df$true_fraction),3),  
                       " COR=", round(cor.test(plot_df$estimate, plot_df$true_fraction)$estimate,2))
  annotations <- data.frame(xpos = c(-Inf), ypos =  c(Inf),
                            annotateText = c(text_to_add),
                            hjustvar = c(-0.2),
                            vjustvar = c(2))
  p <- p + geom_text(data = annotations, aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText), size = annot_text_size)
  return(p)
} 
cor_rmse_barplots <- function(deconv_data, bulk, labels, my_theme, annot_text_size = 7){
  
  plot_df <- deconv_data %>% dplyr::filter(bulk_name %in% bulk) 
  
  plot_df$RMSE = plot_df$Correlation = 0
  plot_df$pvalue = ""
  
  for(method in unique(plot_df$method)){
    for(pred_type in unique(plot_df$predictions_type)){
      tmp <- plot_df[which(plot_df$method == method & plot_df$predictions_type == pred_type), ]
      plot_df$RMSE[which(plot_df$method == method & plot_df$predictions_type == pred_type)] <- round(rmse(tmp$true_fraction, tmp$estimate), 3)
      cor_test <- cor.test(tmp$true_fraction, tmp$estimate)
      if(cor_test$p.value < 0.001){pval = "***"} else if (cor_test$p.value < 0.01){pval = "**"} else if (cor_test$p.value <= 0.05){ pval = "*"} else{pval = ""} 
      
      plot_df$Correlation[which(plot_df$method == method & plot_df$predictions_type == pred_type)] <- round(cor_test$estimate, 3)
      plot_df$pvalue[which(plot_df$method == method & plot_df$predictions_type == pred_type)] <- pval
    } 
  }  
  
  plot_df <- unique(plot_df[, c("method", "bulk_name", "Correlation", "RMSE", "predictions_type", "pvalue")])
  plot_df$method <- sapply(plot_df$method, function(i) unlist(strsplit(i, "_"))[[1]])
  plot_df$method <- factor(plot_df$method, levels = labels)
  
  
  rmse_plot_withUnc <- plot_df[which(plot_df$predictions_type == "With Uncharacterized"),] %>% 
    ggplot(aes(x = reorder(method, RMSE), y = RMSE, fill = factor(ifelse(method == "EPIC-ATAC", "EPIC-ATAC", "Other")) )) +
    geom_bar(stat = "identity", width = 0.8, show.legend = F) + labs(x = "Method") + 
    scale_fill_manual(values=c("red","#A5A5A5")) +
    facet_grid(~ bulk_name, scales = "free", switch = "y")+
    geom_text(aes(y = 0, label = method), hjust = 0, nudge_y = 0.001, angle = 90, color = "black", fontface = "bold", size = annot_text_size) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + theme_linedraw() + my_theme
  
  cor_plot_withUnc <- plot_df[which(plot_df$predictions_type == "With Uncharacterized"),] %>% 
    ggplot(aes(x = reorder(method, -Correlation), y = Correlation, fill = factor(ifelse(method == "EPIC-ATAC", "EPIC-ATAC", "Other")) )) +
    geom_bar(stat = "identity", width = 0.8, show.legend = F) + labs(x = "Method") + 
    scale_fill_manual(values=c("red","#A5A5A5")) +
    facet_grid(~ bulk_name, scales = "free", switch = "y")+
    geom_text(aes(y = 0, label = method), hjust = 0, nudge_y = 0.01, angle = 90, color = "black", fontface = "bold", size = annot_text_size) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + theme_linedraw() + my_theme
  
  cp_withUnc <- cowplot::plot_grid(
    cor_plot_withUnc, NULL, rmse_plot_withUnc, ncol = 3, nrow = 1, label_size = 20,
    rel_widths = c(0.9, 0.05, 0.9)
  )
  
  rmse_plot_withoutUnc<- plot_df[which(plot_df$predictions_type == "Without Uncharacterized"),] %>% 
    ggplot(aes(x = reorder(method, RMSE), y = RMSE, fill = factor(ifelse(method == "EPIC-ATAC", "EPIC-ATAC", "Other")) )) +
    geom_bar(stat = "identity", width = 0.8, show.legend = F) + labs(x = "Method") + 
    scale_fill_manual(values=c("red","#A5A5A5")) +
    facet_grid(~ bulk_name, scales = "free", switch = "y")+
    geom_text(aes(y = 0, label = method), hjust = 0, nudge_y = 0.001, angle = 90, color = "black", fontface = "bold", size = annot_text_size) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + theme_linedraw() + my_theme
  
  cor_plot_withoutUnc <- plot_df[which(plot_df$predictions_type == "Without Uncharacterized"),] %>% 
    ggplot(aes(x = reorder(method, -Correlation), y = Correlation, fill = factor(ifelse(method == "EPIC-ATAC", "EPIC-ATAC", "Other")) )) +
    geom_bar(stat = "identity", width = 0.8, show.legend = F) + labs(x = "Method") + 
    scale_fill_manual(values=c("red","#A5A5A5")) +
    facet_grid(~ bulk_name, scales = "free", switch = "y")+
    geom_text(aes(y = 0, label = method), hjust = 0, nudge_y = 0.01, angle = 90, color = "black", fontface = "bold", size = annot_text_size) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + theme_linedraw() + my_theme 
  
  cp_withoutUnc <- cowplot::plot_grid(
    cor_plot_withoutUnc, NULL, rmse_plot_withoutUnc, ncol = 3, nrow = 1, label_size = 20,
    rel_widths = c(0.9, 0.05, 0.9)
  )
  
  
  return(list(cp_withUnc, cp_withoutUnc, rmse_plot_withUnc, cor_plot_withUnc, rmse_plot_withoutUnc, cor_plot_withoutUnc))
  
} 
benchmark_scatter_plots <- function(summary_predictions, bulk, methods, PBMC = T, text_size = 8, axis_size = 7){
  
  summary_predictions_subset <- summary_predictions[which(summary_predictions$bulk_name == bulk &
                                                            summary_predictions$method %in% methods),]
  summary_predictions_subset$method <- factor(summary_predictions_subset$method, levels = methods)
  
  if(PBMC){
    cors <- plyr::ddply(summary_predictions_subset, c("method"), summarise, cor = round(cor(true_fraction, estimate), 2))
    rmse <- plyr::ddply(summary_predictions_subset, c("method"), summarise, cor = round(Metrics::rmse(true_fraction, estimate), 2))
  }else{
    cors <- plyr::ddply(summary_predictions_subset, c("predictions_type","method"), summarise, cor = round(cor(true_fraction, estimate), 2))
    rmse <- plyr::ddply(summary_predictions_subset, c("predictions_type","method"), summarise, cor = round(Metrics::rmse(true_fraction, estimate), 2))
  }
  
  p <- summary_predictions_subset |> mutate(cell_type := factor(cell_type)) |> mutate(method := factor(method)) |>
    ggplot(aes(x = true_fraction, y = estimate)) +
    geom_point(aes(color = cell_type, shape = sample), size = 2) +
    xlim(0, 1) + ylim(0, 1)
  
  if(PBMC){
    p <- p + facet_grid(~ method, scales = "free", switch = "y")
  }else{
    p <- p + facet_grid(predictions_type ~ method, scales = "free", switch = "y")
  }
  
  p <- p + geom_text(data=cors, aes(label=paste("cor=", cor, sep="")), x=0.2, y=0.95)+
    geom_text(data=rmse, aes(label=paste("rmse=", cor, sep="")), x=0.22, y=0.85)+
    scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% summary_predictions$cell_type)]) +
    scale_shape_manual(values = c(16, 17, 18, 8, 7, 15, 9, 10, 12, 13, 14, 1, 2, 0, 3, 4, 5, 6, 11, 19, 20, 21, 22, 23, 24)) +
    geom_abline(slope = 1, intercept = 0, color="grey") +
    cowplot::background_grid(major="xy") +
    guides(col=guide_legend("Cell types"),
           shape=guide_legend("Sample")) +
    ylab("Predictions") + xlab("True proportions") + theme_bw() +
    theme(axis.text = element_text(size = axis_size, face = "bold"),
          legend.title = element_text(size = text_size, face = "bold"),
          legend.text = element_text(size = text_size),
          axis.title = element_text(size = text_size, face = "bold"),
          strip.text = element_text(size = text_size, face = "bold"))+
    guides(color = guide_legend(ncol = 1,title = NULL), shape = "none")
  
  
  return(p)
}

# Parameters --------------------------------------------------------------

descriminative_colors = c("Bcells"="darkgreen","CD4_Tcells"="#3288bdff","CD8_Tcells"="#4daf4aff","Dendritic"="#984ea3ff",  
                          "Macrophages"="#ff7f00ff","Monocytes"="#ffd92fff","Myeloid"="#a0451fff", "DCs + Macrophages"="#a0451fff", "NKcells"="#70c4b4ff",
                          "NK + T cells"="#bea0ccff", "Tcells"="#660033","Uncharacterized"="#d90017ff", "Neutrophils"="black","Endothelial"="#ffd92fff","Fibroblasts"="gray",
                          "Naive_CD4_Tcells"="#bea0ccff","Non_Naive_CD4_Tcells"="#99FFFF","Tregs"="#3288bdff", "Naive_CD8_Tcells"="#660033","Non_Naive_CD8_Tcells"="#4daf4aff")


benchmarking_output_path = as.character(commandArgs(TRUE)[1])
bulk_path = as.character(commandArgs(TRUE)[2])
fig_path = as.character(commandArgs(TRUE)[3])

dir.create(fig_path)
fig_individual_path = paste0(fig_path, "individual_figures/")
dir.create(fig_individual_path)

#################################################################
##          Individual scatter plots: Fig 3 and Fig 4          ##
#################################################################

datasets <- c("HTAN", "Basal cell carcinoma", "Gynecological cancers", "PBMC pseudobulk", "PBMC experiment")
methods <- c("EPIC-ATAC", "quanTIseq", "ABIS", "DeconPeaker_ourMarkers", "DeconPeaker-Custom", "CIBERSORTx_ourMarkers", "CIBERSORTx-Custom")   
axis_text_size = 25
text_size = 30
legend_text_size = 20
my_theme <- theme(axis.text = element_text(size = axis_text_size, face = "bold"),
                  legend.title = element_text(size = text_size, face = "bold"),
                  legend.text = element_text(size = legend_text_size),
                  axis.title = element_text(size = text_size, face = "bold"),
                  strip.text = element_text(size = text_size, face = "bold"),
                  axis.title.y = element_text(vjust = +1.5),
                  axis.title.x = element_text(vjust = -0.5),
                  legend.key.size = unit(1.5, 'lines')) 


# Load predictions summary without renormalization excluding the uncharacterized cells

summary_predictions = read.table(paste0(benchmarking_output_path, "summary_atac_predictions.txt"), header = T, sep = "\t")
# summary_predictions = read.table("summary_atac_predictions.txt", header = T, sep = "\t")
for(method in methods){ #unique(summary_predictions$method) 
  print(method)
  for(bulk in datasets){
    scatter_name <- paste0(gsub(" ", "_", bulk), "_", gsub(" |-", "_", method))
    if(grepl("PBMC", bulk)){
      annot_size = 7
      xlim = max(summary_predictions$true_fraction[which(summary_predictions$method == method & summary_predictions$bulk_name == bulk)]) + 0.05
      ylim = max(summary_predictions$estimate[which(summary_predictions$method == method & summary_predictions$bulk_name == bulk)]) + 0.05
    }else{
      annot_size = 5
      xlim = 1
      ylim = 1
    }
    assign(x = scatter_name,
           value = deconvolution_scatter_plot(deconv_data = summary_predictions,
                                              method_name = method,
                                              bulk = bulk,
                                              xlim = xlim , ylim = ylim, point_size = 4, 
                                              annot_text_size = annot_size, xlab = "True proportions", ylab = "Predictions"))
  }
}


# Load predictions summary with renormalization excluding the uncharacterized cell
# summary_predictions = read.table("summary_atac_predictions-renorm.txt", header = T, sep = "\t")
summary_predictions = read.table(paste0(benchmarking_output_path,"summary_atac_predictions-renorm.txt"), header = T, sep = "\t")

for(method in methods){
  print(method)
  for(bulk in datasets){
    scatter_name <- paste0(gsub(" ", "_", bulk), "_", gsub(" |-", "_", method), "_renorm")
    if(grepl("PBMC", bulk)){
      annot_size = 7
      xlim = max(summary_predictions$true_fraction[which(summary_predictions$method == method & summary_predictions$bulk_name == bulk)]) + 0.05
      ylim = max(summary_predictions$estimate[which(summary_predictions$method == method & summary_predictions$bulk_name == bulk)]) + 0.05
    }else{
      annot_size = 5
      xlim = 1
      ylim = 1
    }
    assign(x = scatter_name,
           value = deconvolution_scatter_plot(deconv_data = summary_predictions,
                                              method_name = method,
                                              bulk = bulk,
                                              xlim = xlim , ylim = ylim, point_size = 4, annot_text_size = annot_size, xlab = "True proportions", ylab = "Predictions"))
    
  }
}


# Get scatter plots for each cohort of the HTAN data ----------------------

summary_predictions = read.table(paste0(benchmarking_output_path, "summary_atac_predictions.txt"), header = T, sep = "\t")

htan_predictions <- summary_predictions[which(summary_predictions$bulk_name == "HTAN" & summary_predictions$method == "EPIC-ATAC"),]

cors <- plyr::ddply(htan_predictions, c("bulk_name_original"), summarise, cor = round(cor(true_fraction, estimate), 2))
rmse <- plyr::ddply(htan_predictions, c("bulk_name_original"), summarise, cor = round(Metrics::rmse(true_fraction, estimate), 2))

p <- htan_predictions |> mutate(cell_type := factor(cell_type)) |> mutate(bulk_name_original := factor(bulk_name_original)) |>
  ggplot(aes(x = true_fraction, y = estimate)) + scale_shape_manual(values = c(16, 17, 18, 8, 7, 15, 9, 10, 12, 13, 14, 1, 2, 0, 21, 22, 23, 24, 3, 4, 5, 6, 11, 19, 20)) +
  geom_point(aes(color = cell_type, shape = sample), size = 2) + xlim(0,1)  + ylim(0,1) 
p <- p + facet_wrap(~ bulk_name_original, scales = "fixed", nrow = 2)

axis_size = 7
text_size = 8
p <- p + geom_text(data=cors, aes(label=paste("cor=", cor, sep="")), x=0.2, y=0.95)+
  geom_text(data=rmse, aes(label=paste("rmse=", cor, sep="")), x=0.22, y=0.85)+
  scale_color_manual(values = descriminative_colors[which(names(descriminative_colors) %in% htan_predictions$cell_type)]) +
  scale_shape_manual(values = c(16, 17, 18, 8, 7, 15, 9, 10, 12, 13, 14, 1, 2, 0, 21, 22, 23, 24, 3, 4, 5, 6, 11, 19, 20)) +
  geom_abline(slope = 1, intercept = 0, color="grey") +
  cowplot::background_grid(major="xy") +
  guides(col=guide_legend("Cell types"),
         shape=guide_legend("Sample")) +
  ylab("Predictions") + xlab("True proportions") + theme_bw() +
  theme(axis.text = element_text(size = axis_size, face = "bold"),
        legend.title = element_text(size = text_size, face = "bold"),
        legend.text = element_text(size = text_size),
        axis.title = element_text(size = text_size, face = "bold"),
        strip.text = element_text(size = text_size, face = "bold"))+
  guides(color = guide_legend(ncol = 1,title = NULL), shape = "none")

pdf(paste0(fig_path, "Figure4_sup1.pdf"), w = 11, h = 5)
p
dev.off()

#################################################################
##                RMSE and correlation barplots                ##
#################################################################

summary_predictions_abs <- read.table(paste0(benchmarking_output_path,"summary_atac_predictions.txt"), header = T, sep = "\t")
summary_predictions_abs <- summary_predictions_abs[which(summary_predictions_abs$method %in% c("EPIC-ATAC","quanTIseq","ABIS","DeconPeaker_ourMarkers","DeconPeaker-Custom","CIBERSORTx_ourMarkers","CIBERSORTx-Custom")),]
summary_predictions_abs$predictions_type <- "With Uncharacterized"

summary_predictions_relative <- read.table(paste0(benchmarking_output_path,"summary_atac_predictions-renorm.txt"), header = T, sep = "\t")
summary_predictions_relative <- summary_predictions_relative[which(summary_predictions_relative$method %in% c("EPIC-ATAC","quanTIseq","ABIS","DeconPeaker_ourMarkers","DeconPeaker-Custom","CIBERSORTx_ourMarkers","CIBERSORTx-Custom")),]
summary_predictions_relative$predictions_type <- "Without Uncharacterized"
summary_predictions <- rbind(summary_predictions_abs, summary_predictions_relative)

# Get Average RMSE 
summary_predictions_subset <- summary_predictions_abs[which(summary_predictions_abs$method == "EPIC-ATAC" & summary_predictions_abs$predictions_type == "With Uncharacterized"),]
RMSE <- cor <- vector()
for(cell_type in unique(summary_predictions_subset$cell_type)){
  tmp <- summary_predictions_subset[which(summary_predictions_subset$cell_type== cell_type), ]
  RMSE <- c(RMSE, round(rmse(tmp$true_fraction, tmp$estimate), 3))
  cor_test <- cor.test(tmp$true_fraction, tmp$estimate)
  cor <- c(cor, round(cor_test$estimate, 3))
}
names(RMSE) <- names(cor) <- unique(summary_predictions_subset$cell_type)
mean(RMSE)


axis_text_size = 30
text_size = 35

my_theme <- theme(
  axis.text = element_text(size = axis_text_size, face = "bold"),
  axis.title = element_text(size = text_size, face = "bold"),
  strip.text = element_text(size = text_size, face = "bold", colour = "black"),
  axis.text.x = element_blank(),  
  axis.ticks.x = element_blank(),
  strip.background = element_rect(fill = "grey"),
  axis.title.y = element_text(vjust = +1.5),
  axis.title.x = element_text(vjust = -0.5)
) 


for(bulk in datasets){
  barplots <- cor_rmse_barplots(deconv_data = summary_predictions[which(summary_predictions$method %in% c("EPIC-ATAC", "quanTIseq", "DeconPeaker_ourMarkers", "CIBERSORTx_ourMarkers", "ABIS", "DeconPeaker-Custom", "CIBERSORTx-Custom")),] ,
                                bulk = bulk,
                                labels = c("EPIC-ATAC", "quanTIseq", "CIBERSORTx", "DeconPeaker", "ABIS", "DeconPeaker-Custom" , "CIBERSORTx-Custom"), 
                                my_theme = my_theme, annot_text_size = 6)
  
  assign(x = paste0(gsub(" ", "_", bulk), "_barPlot_withUnc"), value = barplots[[1]])
  assign(x = paste0(gsub(" ", "_", bulk), "_barPlot_withoutUnc"), value = barplots[[2]])
  
  assign(x = paste0(gsub(" ", "_", bulk), "_RMSEbarPlot_withUnc"), value = barplots[[3]])
  assign(x = paste0(gsub(" ", "_", bulk), "_CORbarPlot_withUnc"), value = barplots[[4]])
  assign(x = paste0(gsub(" ", "_", bulk), "_RMSEbarPlot_withoutUnc"), value = barplots[[5]])
  assign(x = paste0(gsub(" ", "_", bulk), "_CORbarPlot_withoutUnc"), value = barplots[[6]])
  
} 


#################################################################
##                        Figure 3 grid                        ##
#################################################################
axis_text_size = 19
text_size = 25
legend_text_size = 15
strip_text_size = 19
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
  axis.text.x = element_blank(),  
  axis.ticks.x = element_blank(),
  strip.background = element_rect(fill = "grey"),
  axis.title.y = element_text(vjust = +1.5),
  axis.title.x = element_text(vjust = -0.5),
  legend.title = element_text(size = text_size, face = "bold"),
  legend.text = element_text(size = legend_text_size),
  legend.key.size = unit(1.5, 'lines'),
  panel.grid.minor = element_line(size = 0.2,colour = "darkgray"), panel.grid.major = element_line(size = 0.2,colour = "darkgray")
) 


PBMC_benchmark_figure <- cowplot::plot_grid(
  PBMC_experiment_EPIC_ATAC + theme(legend.position = c(0.8,0.31), legend.background = element_rect(size=0.3, linetype="solid", colour ="black")) + xlab("True proportions") + my_theme_scatter, 
  NULL,
  PBMC_experiment_CORbarPlot_withUnc + my_theme_barplot, NULL, PBMC_experiment_RMSEbarPlot_withUnc + my_theme_barplot,  
  NULL, NULL, NULL, NULL, NULL,
  PBMC_pseudobulk_EPIC_ATAC + theme(legend.position = c(0.8,0.31), legend.background = element_rect(size=0.3, linetype="solid", colour ="black")) + xlab("True proportions") + my_theme_scatter, 
  NULL,
  PBMC_pseudobulk_CORbarPlot_withUnc + my_theme_barplot, NULL, PBMC_pseudobulk_RMSEbarPlot_withUnc + my_theme_barplot,
  ncol = 5, nrow = 3, 
  rel_heights = c(0.9, 0.1, 0.9), rel_widths = c(0.6, 0.01, 0.4, 0.01, 0.4)
)


pdf(paste0(fig_path, "Figure3_grid.pdf"), w = 13, h = 10)
PBMC_benchmark_figure
dev.off()


#################################################################
##                        Figure 4 grid                        ##
#################################################################

axis_text_size = 17
text_size = 16
legend_text_size = 15

my_theme_scatter <- theme(
  axis.text = element_text(size = axis_text_size, face = "bold"),
  axis.title = element_text(size = text_size, face = "bold"),
  strip.text = element_text(size = text_size, face = "bold", colour = "black"),
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
  strip.text = element_text(size = text_size, face = "bold", colour = "black"),
  axis.text.x = element_blank(),  
  axis.ticks.x = element_blank(),
  strip.background = element_rect(fill = "grey"),
  axis.title.y = element_text(vjust = +1.5),
  axis.title.x = element_text(vjust = -0.5),
  legend.title = element_text(size = text_size, face = "bold"),
  legend.text = element_text(size = legend_text_size),
  legend.key.size = unit(1.5, 'lines'),
  panel.grid.minor = element_line(size = 0.2,colour = "darkgray"), panel.grid.major = element_line(size = 0.2,colour = "darkgray")
) 


Basal_cell_carcinoma_scatter <- cowplot::plot_grid(
  Basal_cell_carcinoma_EPIC_ATAC + theme(legend.position = 'none')+ my_theme_scatter, NULL, 
  cowplot::get_legend(Basal_cell_carcinoma_EPIC_ATAC + guides(color = guide_legend(ncol = 1, title = NULL))+ my_theme_scatter),
  ncol = 3, nrow = 1, rel_widths =  c(0.7, 0.01, 0.3)
)

Gynecological_cancers_scatter <- cowplot::plot_grid(
  Gynecological_cancers_EPIC_ATAC + theme(legend.position = 'none')+ my_theme_scatter, NULL, 
  cowplot::get_legend(Gynecological_cancers_EPIC_ATAC + guides(color = guide_legend(ncol = 1, title = NULL))+ my_theme_scatter),
  ncol = 3, nrow = 1, rel_widths =  c(0.7, 0.01, 0.3)
)

HTAN_scatter <- cowplot::plot_grid(
  HTAN_EPIC_ATAC + theme(legend.position = 'none')+ my_theme_scatter, NULL, 
  cowplot::get_legend(HTAN_EPIC_ATAC + guides(color = guide_legend(ncol = 1, title = NULL))+ my_theme_scatter),
  ncol = 3, nrow = 1, rel_widths =  c(0.7, 0.01, 0.3)
)


pdf(paste0(fig_path, "Figure4_grid.pdf"), w = 20, h = 12)
cowplot::plot_grid(
  Basal_cell_carcinoma_scatter, NULL,
  Basal_cell_carcinoma_CORbarPlot_withUnc + my_theme_barplot, NULL, Basal_cell_carcinoma_RMSEbarPlot_withUnc + my_theme_barplot, NULL,
  Basal_cell_carcinoma_CORbarPlot_withoutUnc + my_theme_barplot, NULL, Basal_cell_carcinoma_RMSEbarPlot_withoutUnc + my_theme_barplot, 
  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
  Gynecological_cancers_scatter, NULL,
  Gynecological_cancers_CORbarPlot_withUnc + my_theme_barplot, NULL, Gynecological_cancers_RMSEbarPlot_withUnc + my_theme_barplot, NULL,
  Gynecological_cancers_CORbarPlot_withoutUnc + my_theme_barplot, NULL, Gynecological_cancers_RMSEbarPlot_withoutUnc + my_theme_barplot, 
  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
  HTAN_scatter, NULL,
  HTAN_CORbarPlot_withUnc + my_theme_barplot, NULL, HTAN_RMSEbarPlot_withUnc + my_theme_barplot, NULL,
  HTAN_CORbarPlot_withoutUnc + my_theme_barplot, NULL, HTAN_RMSEbarPlot_withoutUnc + my_theme_barplot, 
  ncol = 9, nrow = 5, 
  rel_heights = c(0.9, 0.1, 0.9, 0.1, 0.9), rel_widths = c(0.7, 0.01, 0.4, 0.01, 0.4, 0.01, 0.4, 0.01, 0.4)
)
dev.off()



# Sup figure all scatter plots -----------------------------------------------------------


# PBMC pseudobulks 

PBMC_pseudobulk_benchmark_scatter <- benchmark_scatter_plots(summary_predictions_relative, bulk = "PBMC pseudobulk", methods = c("EPIC-ATAC", "quanTIseq", "DeconPeaker_ourMarkers", "CIBERSORTx_ourMarkers", "ABIS"),
                                                             PBMC = T, text_size = 8, axis_size = 7)
PBMC_pseudobulk_benchmark_scatter_newRefs <- benchmark_scatter_plots(summary_predictions_relative, bulk = "PBMC pseudobulk", methods = c("CIBERSORTx-Custom", "DeconPeaker-Custom"), PBMC = T, text_size = 8, axis_size = 7)

# PBMC_experiment 
PBMC_experiment_benchmark_scatter <- benchmark_scatter_plots(summary_predictions_relative, bulk = "PBMC experiment", methods = c("EPIC-ATAC", "quanTIseq", "DeconPeaker_ourMarkers", "CIBERSORTx_ourMarkers", "ABIS"),
                                                             PBMC = T, text_size = 8, axis_size = 7)
PBMC_experiment_benchmark_scatter_newRefs <- benchmark_scatter_plots(summary_predictions_relative, bulk = "PBMC experiment", methods = c("CIBERSORTx-Custom", "DeconPeaker-Custom"),PBMC = T, text_size = 8, axis_size = 7)


# BCC pseudobulks 
Basal_cell_carcinoma_benchmark_scatter <- benchmark_scatter_plots(summary_predictions, bulk = "Basal cell carcinoma", methods = c("EPIC-ATAC", "quanTIseq", "DeconPeaker_ourMarkers", "CIBERSORTx_ourMarkers", "ABIS"),
                                                                  PBMC = F, text_size = 8, axis_size = 7)
Basal_cell_carcinoma_benchmark_scatter_newRefs <- benchmark_scatter_plots(summary_predictions, bulk = "Basal cell carcinoma", methods = c("CIBERSORTx-Custom", "DeconPeaker-Custom"), PBMC = F, text_size = 8, axis_size = 7)


# Gynecological_cancers 
Gynecological_cancers_benchmark_scatter <- benchmark_scatter_plots(summary_predictions, bulk = "Gynecological cancers", methods = c("EPIC-ATAC", "quanTIseq", "DeconPeaker_ourMarkers", "CIBERSORTx_ourMarkers", "ABIS"),
                                                                   PBMC = F, text_size = 8, axis_size = 7)
Gynecological_cancers_benchmark_scatter_newRefs <- benchmark_scatter_plots(summary_predictions, bulk = "Gynecological cancers", methods = c("CIBERSORTx-Custom", "DeconPeaker-Custom"), PBMC = F, text_size = 8, axis_size = 7)

# HTAN samples
HTAN_benchmark_scatter <- benchmark_scatter_plots(summary_predictions, bulk = "HTAN", methods = c("EPIC-ATAC", "quanTIseq", "DeconPeaker_ourMarkers", "CIBERSORTx_ourMarkers", "ABIS"),
                                                  PBMC = F, text_size = 8, axis_size = 7)
HTAN_benchmark_scatter_newRefs <- benchmark_scatter_plots(summary_predictions, bulk = "HTAN", methods = c("CIBERSORTx-Custom", "DeconPeaker-Custom"), PBMC = F, text_size = 8, axis_size = 7)



# Save sup figures  
PBMC_benchmark_sup1 <- cowplot::plot_grid(cowplot::get_legend(PBMC_experiment_benchmark_scatter + theme(legend.position="bottom") + guides(color = guide_legend(nrow = 1,title = NULL))), NULL,
                                          PBMC_experiment_benchmark_scatter + theme(legend.position = 'hidden'), NULL,
                                          PBMC_pseudobulk_benchmark_scatter + theme(legend.position = 'hidden'),
                                          label_size = 20, ncol = 1, nrow = 5, rel_heights = c(0.1, 0.05, 0.9, 0.05, 0.9))

cancer_benchmark_sup1 <- cowplot::plot_grid(Basal_cell_carcinoma_benchmark_scatter + theme(legend.position = 'hidden'),
                                            cowplot::get_legend(Basal_cell_carcinoma_benchmark_scatter + theme(legend.position="bottom") + guides(color = guide_legend(nrow = 1,title = NULL))), NULL,
                                            Gynecological_cancers_benchmark_scatter + theme(legend.position = 'hidden'),
                                            cowplot::get_legend(Gynecological_cancers_benchmark_scatter + theme(legend.position="bottom") + guides(color = guide_legend(nrow = 1,title = NULL))), NULL,
                                            HTAN_benchmark_scatter + theme(legend.position = 'hidden'),
                                            cowplot::get_legend(HTAN_benchmark_scatter + theme(legend.position="bottom") + guides(color = guide_legend(nrow = 1,title = NULL))),
                                            label_size = 20, ncol = 1, nrow = 8, rel_heights = c(0.9, 0.1, 0.05, 0.9, 0.1, 0.05, 0.9, 0.1))


PBMC_benchmark_sup2 <- cowplot::plot_grid(cowplot::get_legend(PBMC_experiment_benchmark_scatter_newRefs + theme(legend.position="bottom") + guides(color = guide_legend(nrow = 2,title = NULL))), NULL,
                                          PBMC_experiment_benchmark_scatter_newRefs + theme(legend.position = 'hidden'), NULL,
                                          PBMC_pseudobulk_benchmark_scatter_newRefs + theme(legend.position = 'hidden'),
                                          label_size = 20, ncol = 1, nrow = 5, rel_heights = c(0.1, 0.05, 0.9, 0.05, 0.9))

cancer_benchmark_sup2 <- cowplot::plot_grid(Basal_cell_carcinoma_benchmark_scatter_newRefs + theme(legend.position = 'hidden'),
                                            cowplot::get_legend(Basal_cell_carcinoma_benchmark_scatter_newRefs + theme(legend.position="bottom") + guides(color = guide_legend(nrow = 2,title = NULL))), NULL,
                                            Gynecological_cancers_benchmark_scatter_newRefs + theme(legend.position = 'hidden'),
                                            cowplot::get_legend(Gynecological_cancers_benchmark_scatter_newRefs + theme(legend.position="bottom") + guides(color = guide_legend(nrow = 2,title = NULL))),NULL,
                                            HTAN_benchmark_scatter_newRefs + theme(legend.position = 'hidden'),
                                            cowplot::get_legend(HTAN_benchmark_scatter_newRefs + theme(legend.position="bottom") + guides(color = guide_legend(nrow = 2,title = NULL))),
                                            label_size = 20, ncol = 1, nrow = 8, rel_heights = c(0.9, 0.1, 0.05, 0.9, 0.1, 0.05, 0.9, 0.1))

pdf(paste0(fig_path, "Figure3_sup2A.pdf"), w = 11, h = 5)
PBMC_benchmark_sup1
dev.off()

pdf(paste0(fig_path, "Figure4_sup2.pdf"), w = 11, h = 14)
cancer_benchmark_sup1
dev.off()

pdf(paste0(fig_path, "Figure3_sup2B.pdf"), w = 5, h = 5.5)
PBMC_benchmark_sup2
dev.off()

pdf(paste0(fig_path, "Figure4_sup3.pdf"), w = 5, h = 15)
cancer_benchmark_sup2
dev.off()



##################################################################
##                    Benchmarking boxplots                     ##
##################################################################

my_theme <- theme_bw() + theme(text = element_text(size = 12, face = "bold"),
                               axis.text = element_text(size = 12, face = "bold"),
                               legend.title = element_text(size = 12, face = "bold"),
                               legend.text = element_text(size = 12),
                               axis.title = element_text(size = 12, face = "bold"),
                               strip.text = element_text(size = 12, face = "bold"),
                               axis.text.x = element_text(angle=45, hjust=1, size = 12))

# For PBMC samples ---------------------------------------------------------

summary_predictions_abs <- read.table(paste0(benchmarking_output_path,"summary_atac_predictions.txt"), header = T, sep = "\t")
summary_predictions_abs <- summary_predictions_abs[which(summary_predictions_abs$bulk_name %in% c("PBMC pseudobulk","PBMC experiment") &
                                                           summary_predictions_abs$method %in% c("EPIC-ATAC","quanTIseq","ABIS","DeconPeaker_ourMarkers","DeconPeaker-Custom","CIBERSORTx_ourMarkers","CIBERSORTx-Custom")),]
summary_predictions_abs$RMSE = summary_predictions_abs$Correlation = 0
summary_predictions_abs = summary_predictions_abs[-which(summary_predictions_abs$method == "quanTIseq" & summary_predictions_abs$cell_type == "otherCells"),]
for(method in unique(summary_predictions_abs$method)){
  for(cell_type in unique(summary_predictions_abs$cell_type)[which(unique(summary_predictions_abs$cell_type) != "otherCells")] ){
    
    tmp <- summary_predictions_abs[which(summary_predictions_abs$method == method &
                                           summary_predictions_abs$cell_type == cell_type), ]
    
    summary_predictions_abs$RMSE[which(summary_predictions_abs$method == method &
                                         summary_predictions_abs$cell_type == cell_type)] <- round(rmse(tmp$true_fraction, tmp$estimate), 3)
    
    cor_test <- cor.test(tmp$true_fraction, tmp$estimate)
    summary_predictions_abs$Correlation[which(summary_predictions_abs$method == method &
                                                summary_predictions_abs$cell_type == cell_type)] <- round(cor_test$estimate, 3)
  }
  
}
summary_predictions_abs <- unique(summary_predictions_abs[, c("method", "bulk_name", "Correlation", "RMSE", "cell_type")])

my_comparisons <- list( c("EPIC-ATAC", "quanTIseq"), c("EPIC-ATAC", "ABIS"),c("EPIC-ATAC", "DeconPeaker_ourMarkers"),  c("EPIC-ATAC", "DeconPeaker-Custom"),
                        c("EPIC-ATAC", "CIBERSORTx_ourMarkers"),c("EPIC-ATAC", "CIBERSORTx-Custom"))

PBMC_boxplots_cor <- ggpubr::ggboxplot(summary_predictions_abs, x = "method", y = "Correlation")+ geom_jitter(position = "dodge", size = 2, aes(color = cell_type)) +
  scale_color_manual(values = descriminative_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.14))) +
  ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.format", paired = T,  vjust =  0.2, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                             symbols = c("****", "***", "**", "*", "ns"))) +
  labs(x = "", y = "Correlation")+
  my_theme + theme(legend.position = "none")

PBMC_boxplots_rmse <- ggpubr::ggboxplot(summary_predictions_abs, x = "method", y = "RMSE")+ geom_jitter(position = "dodge", size = 2, aes(color = cell_type)) +
  scale_color_manual(values = descriminative_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.14))) +
  ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.format", paired = T,  vjust =  0.2,symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                            symbols = c("****", "***", "**", "*", "ns"))) +
  labs(x = "", y = "RMSE")+
  my_theme + theme(legend.position = "none")


# For cancer pseudobulks --------------------------------------------------

# Only considering the case without cancer cells and both datasets combined

summary_predictions_relative <- read.table(paste0(benchmarking_output_path,"summary_atac_predictions-renorm.txt"), header = T, sep = "\t")
summary_predictions_relative <- summary_predictions_relative[which(summary_predictions_relative$method %in% c("EPIC-ATAC", "quanTIseq", "ABIS", "DeconPeaker_ourMarkers", 
                                                                                                              "DeconPeaker-Custom", "CIBERSORTx_ourMarkers", "CIBERSORTx-Custom") &
                                                                     summary_predictions_relative$bulk_name %in% c("Basal cell carcinoma", "Gynecological cancers", "HTAN")),]
summary_predictions_relative$RMSE = summary_predictions_relative$Correlation = 0

for(bulk in c("Basal cell carcinoma", "Gynecological cancers", "HTAN")){
  summary_predictions <- summary_predictions_relative[which(summary_predictions_relative$bulk_name == bulk),]
  
  for(method in unique(summary_predictions$method)){
    for(cell_type in unique(summary_predictions$cell_type[which(summary_predictions$method == method)] )){
      
      tmp <- summary_predictions[which(summary_predictions$method == method &
                                         summary_predictions$cell_type == cell_type), ]
      if(anyNA(tmp$estimate)){
        summary_predictions_relative$RMSE[which(summary_predictions_relative$method == method &
                                                  summary_predictions_relative$bulk_name == bulk &
                                                  summary_predictions_relative$cell_type == cell_type)] <- NA
        
        summary_predictions_relative$Correlation[which(summary_predictions_relative$method == method &
                                                         summary_predictions_relative$bulk_name == bulk &
                                                         summary_predictions_relative$cell_type == cell_type)] <- NA
      } else{
        summary_predictions_relative$RMSE[which(summary_predictions_relative$method == method &
                                                  summary_predictions_relative$bulk_name == bulk &
                                                  summary_predictions_relative$cell_type == cell_type)] <- round(rmse(tmp$true_fraction, tmp$estimate), 3)
        
        cor_test <- cor.test(tmp$true_fraction, tmp$estimate)
        summary_predictions_relative$Correlation[which(summary_predictions_relative$method == method &
                                                         summary_predictions_relative$bulk_name == bulk &
                                                         summary_predictions_relative$cell_type == cell_type)] <- round(cor_test$estimate, 3)
      } 
      
    }
    
  }
  
}

summary_predictions_relative <- unique(summary_predictions_relative[, c("method", "bulk_name", "Correlation", "RMSE", "cell_type")])

my_comparisons <- list( c("EPIC-ATAC", "quanTIseq"), c("EPIC-ATAC", "ABIS"), c("EPIC-ATAC", "DeconPeaker_ourMarkers"),  c("EPIC-ATAC", "DeconPeaker-Custom"),
                        c("EPIC-ATAC", "CIBERSORTx_ourMarkers"), c("EPIC-ATAC", "CIBERSORTx-Custom"))

boxplots_cor <- ggpubr::ggboxplot(summary_predictions_relative, x = "method", y = "Correlation")+ geom_jitter(position = "dodge", size = 2, aes(color = cell_type, shape = bulk_name)) +
  scale_color_manual(values = descriminative_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.14))) +
  ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.format", paired = T,  vjust =  0.2,symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                            symbols = c("****", "***", "**", "*", "ns"))) +
  labs(x = "", y = "Correlation")+
  my_theme + theme(legend.position = "none")

boxplots_rmse <- ggpubr::ggboxplot(summary_predictions_relative, x = "method", y = "RMSE")+ geom_jitter(position = "dodge", size = 2, aes(color = cell_type, shape = bulk_name)) +
  scale_color_manual(values = descriminative_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.14))) +
  ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.format", paired = T,  vjust =  0.2,symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                            symbols = c("****", "***", "**", "*", "ns"))) +
  labs(x = "", y = "RMSE")+
  my_theme + theme(legend.position = "none")

PBMC_benchmark_sup3 <- cowplot::plot_grid(PBMC_boxplots_cor,  NULL,
                                          PBMC_boxplots_rmse,
                                          label_size = 20, ncol = 1, nrow = 3, rel_heights = c(0.9, 0.05, 0.9))

cancer_benchmark_sup3 <- cowplot::plot_grid(boxplots_cor,  NULL,
                                            boxplots_rmse,
                                            label_size = 20, ncol = 1, nrow = 3, rel_heights = c(0.9, 0.01, 0.9))

pdf(paste0(fig_path, "Figure3_sup2C.pdf"), w=5, h = 10)
PBMC_benchmark_sup3
dev.off()


boxplots_rmse_noCancer <- boxplots_rmse

# Considering the case with cancer cells and both datasets combined

summary_predictions_abs <- read.table(paste0(benchmarking_output_path,"summary_atac_predictions.txt"), header = T, sep = "\t")
summary_predictions_abs <- summary_predictions_abs[which(summary_predictions_abs$method %in% c("EPIC-ATAC", "quanTIseq", "ABIS", "DeconPeaker_ourMarkers", 
                                                                                               "DeconPeaker-Custom", "CIBERSORTx_ourMarkers", "CIBERSORTx-Custom") &
                                                           summary_predictions_abs$bulk_name %in% c("Basal cell carcinoma", "Gynecological cancers", "HTAN")),]
summary_predictions_abs$RMSE = summary_predictions_abs$Correlation = 0

for(bulk in c("Basal cell carcinoma", "Gynecological cancers", "HTAN")){
  summary_predictions <- summary_predictions_abs[which(summary_predictions_abs$bulk_name == bulk),]
  
  for(method in unique(summary_predictions$method)){
    for(cell_type in unique(summary_predictions$cell_type[which(summary_predictions$method == method)] )){
      
      tmp <- summary_predictions[which(summary_predictions$method == method &
                                         summary_predictions$cell_type == cell_type), ]
      if(anyNA(tmp$estimate)){
        summary_predictions_abs$RMSE[which(summary_predictions_abs$method == method &
                                             summary_predictions_abs$bulk_name == bulk &
                                             summary_predictions_abs$cell_type == cell_type)] <- NA
        
        summary_predictions_abs$Correlation[which(summary_predictions_abs$method == method &
                                                    summary_predictions_abs$bulk_name == bulk &
                                                    summary_predictions_abs$cell_type == cell_type)] <- NA
      } else{
        summary_predictions_abs$RMSE[which(summary_predictions_abs$method == method &
                                             summary_predictions_abs$bulk_name == bulk &
                                             summary_predictions_abs$cell_type == cell_type)] <- round(rmse(tmp$true_fraction, tmp$estimate), 3)
        
        cor_test <- cor.test(tmp$true_fraction, tmp$estimate)
        summary_predictions_abs$Correlation[which(summary_predictions_abs$method == method &
                                                    summary_predictions_abs$bulk_name == bulk &
                                                    summary_predictions_abs$cell_type == cell_type)] <- round(cor_test$estimate, 3)
      } 
      
    }
    
  }
  
}

summary_predictions_abs <- unique(summary_predictions_abs[, c("method", "bulk_name", "Correlation", "RMSE", "cell_type")])

my_comparisons <- list( c("EPIC-ATAC", "quanTIseq"), c("EPIC-ATAC", "ABIS"), c("EPIC-ATAC", "DeconPeaker_ourMarkers"),  c("EPIC-ATAC", "DeconPeaker-Custom"),
                        c("EPIC-ATAC", "CIBERSORTx_ourMarkers"), c("EPIC-ATAC", "CIBERSORTx-Custom"))

boxplots_rmse_withCancer <- ggpubr::ggboxplot(summary_predictions_abs, x = "method", y = "RMSE", outlier.shape = NA)+ geom_jitter(position = "dodge", size = 2, aes(color = cell_type, shape = bulk_name)) +
  scale_color_manual(values = descriminative_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.14))) +
  ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.format", paired = T,  vjust =  0.2,symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                            symbols = c("****", "***", "**", "*", "ns"))) +
  labs(x = "", y = "RMSE", color = "Cell types", shape = "Dataset")+
  my_theme 


pdf(paste0(fig_path, "Figure4_sup4.pdf"), w=8, h = 10)

plt_no_legend <- cowplot::plot_grid(boxplots_rmse_withCancer + theme(legend.position = "none"),  NULL,
                                    boxplots_rmse_noCancer ,
                                    label_size = 20, ncol = 1, nrow = 3, rel_heights = c(0.9, 0.01, 0.9))
cowplot::plot_grid(plt_no_legend,  NULL,
                   cowplot::get_legend(boxplots_rmse_withCancer) ,
                   label_size = 20, ncol = 3, nrow = 1, rel_widths = c(0.8, 0.05, 0.2))
dev.off()


##################################################################
##                  Benchmarking per cell-type                  ##
##################################################################
axis_text_size = 11; text_size = 9
pdf(paste0(fig_individual_path, "banchmarking_per_cell_type.pdf"), w = 14, h = 5)
for(bulk in c("PBMC", "Basal cell carcinoma", "Gynecological cancers", "HTAN")){
  summary_predictions = read.table(paste0(benchmarking_output_path, "summary_atac_predictions.txt"), header = T, sep = "\t")
  summary_predictions_subset <- summary_predictions[grepl(bulk, summary_predictions$bulk_name), ]
  summary_predictions_subset <- summary_predictions_subset[which(summary_predictions_subset$bulk_name != "PBMC_multiome_simulated" &
                                                                   !summary_predictions_subset$cell_type %in% c("Ery","CLP","GMP","CMP","MEP","MPP","LMPP","HSC")),] 

  rmse_matrix = cor_matrix = cor_text = matrix(nrow = 9, ncol = length(unique(summary_predictions_subset$cell_type)),
                                               dimnames = list(c("EPIC-ATAC", "quanTIseq", "DeconPeaker_ourMarkers", "DeconPeaker-Custom", "DeconPeaker_originalMarkers", "CIBERSORTx_ourMarkers", "CIBERSORTx-Custom", "ABIS", "MCPcounter"),
                                                               unique(summary_predictions_subset$cell_type)))
  for(method in rownames(rmse_matrix)){
    if(method %in% c("EPIC-ATAC", "quanTIseq", "MCPcounter")){
      summary_predictions = read.table(paste0(benchmarking_output_path, "summary_atac_predictions.txt"), header = T, sep = "\t")
      summary_predictions_subset <- summary_predictions[grepl(bulk, summary_predictions$bulk_name), ]
      summary_predictions_subset <- summary_predictions_subset[which(summary_predictions_subset$bulk_name != "PBMC_multiome_simulated" &
                                                                       !summary_predictions_subset$cell_type %in% c("Ery","CLP","GMP","CMP","MEP","MPP","LMPP","HSC")),] 
      
    } else{
      summary_predictions = read.table(paste0(benchmarking_output_path, "summary_atac_predictions-renorm.txt"), header = T, sep = "\t")
      summary_predictions_subset <- summary_predictions[grepl(bulk, summary_predictions$bulk_name), ]
      summary_predictions_subset <- summary_predictions_subset[which(summary_predictions_subset$bulk_name != "PBMC_multiome_simulated" &
                                                                       !summary_predictions_subset$cell_type %in% c("Ery","CLP","GMP","CMP","MEP","MPP","LMPP","HSC")),] 
      
    }
    
    for(cell_type in unique(summary_predictions_subset$cell_type)){
      summary_predictions_method <- summary_predictions_subset[which(summary_predictions_subset$method == method & summary_predictions_subset$cell_type == cell_type), ]
      if(nrow(summary_predictions_method) <= 1 ){
        rmse_matrix[method, cell_type] = NA
        cor_matrix[method, cell_type] = NA
      } else{
        if(method %in% c("MCPcounter")){
          rmse_matrix[method, cell_type] = NA
        } else{
          rmse_matrix[method, cell_type] = round(Metrics::rmse(summary_predictions_method$estimate, summary_predictions_method$true_fraction),2)
        }
        if(sum(summary_predictions_method$true_fraction)!=0 & sum(summary_predictions_method$estimate, na.rm = T)!=0 & length(summary_predictions_method$estimate[which(!is.na(summary_predictions_method$estimate))] ) >= 3){
          correlation_test = cor.test(summary_predictions_method$estimate, summary_predictions_method$true_fraction, alternative = "two.sided", method = "pearson")
          if(is.na(correlation_test$estimate)){
            cor_matrix[method, cell_type] = NA
          } else{
            if(correlation_test$p.value < 0.001){pval = "***"} else if (correlation_test$p.value < 0.01){pval = "**"} else if (correlation_test$p.value <= 0.05){ pval = "*"} else{pval = ""}
            cor_text[method, cell_type] = paste0(round(correlation_test$estimate,2), "\n", pval)
            cor_matrix[method, cell_type] = round(correlation_test$estimate,2)
          }
          
        } else{
          cor_matrix[method, cell_type] = NA
          
        }
        
      }
    }
  }
  
  if(grepl("PBMC", bulk)){
    cor_matrix = cor_matrix[, -which(colnames(cor_matrix) == "otherCells")]
    rmse_matrix = rmse_matrix[, -which(colnames(rmse_matrix) == "otherCells")]
  }
  col_fun <- colorRampPalette(c("#FC2232","#FFD85B","white"))(256)
  cn = colnames(rmse_matrix)
  ht_rmse <- Heatmap(rmse_matrix, col = col_fun, show_row_names = T, row_names_gp = gpar(fontsize = axis_text_size), column_names_gp = gpar(fontsize = axis_text_size),
                     show_column_names = F, row_names_side = "left",
                     bottom_annotation = HeatmapAnnotation(
                       text = anno_text(cn, rot = 45)
                     ),
                     row_title= NULL, column_title = "RMSE", cluster_rows  = F, cluster_columns = F,#  row_labels = spread_format$method,
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(rmse_matrix[i, j], x, y, gp = gpar(fontsize = text_size))
                     },
                     na_col = "snow3", row_split = c(1,2,3,3,3,4,4,5,6),# 1:nrow(rmse_matrix),
                     column_split = 1:ncol(rmse_matrix), border = TRUE ,heatmap_legend_param = list(title = "RMSE"))
  
  col_fun <- colorRampPalette(c("white","#FFD85B","#FC2232"))(256) #"#4A9FD4","white","#FC2232"
  cn = colnames(cor_matrix)
  ht_cor <- Heatmap(cor_matrix, col = col_fun, show_row_names = T, row_names_gp = gpar(fontsize = axis_text_size), column_names_gp = gpar(fontsize = axis_text_size),
                    show_column_names = F, row_names_side = "left",
                    bottom_annotation = HeatmapAnnotation(
                      text = anno_text(cn, rot = 45)
                    ),
                    row_title= NULL, column_title = "Correlation", cluster_rows  = F, cluster_columns = F,#  row_labels = spread_format$method,
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(cor_text[i, j], x, y, gp = gpar(fontsize = text_size))
                    },
                    na_col = "snow3", row_split = c(1,2,3,3,3,4,4,5,6), column_split = 1:ncol(cor_matrix), border = TRUE ,heatmap_legend_param = list(title = "Correlation"))
  
  combined_heatmap <- ht_cor + ht_rmse
  draw(combined_heatmap, column_title = bulk, ht_gap = unit(2,"cm"))
  assign(x = paste0(gsub(" ", "_", bulk), "_cor_rmse_perCellType"),
         value = draw(combined_heatmap, ht_gap = unit(2,"cm")))
  
}
dev.off()

PBMC_benchmark_sup4 <-cowplot::plot_grid(grid.grabExpr(ComplexHeatmap::draw(PBMC_cor_rmse_perCellType)),
                                         label_size = 20, ncol = 1, nrow = 1)

cancer_benchmark_sup4 <-cowplot::plot_grid(grid.grabExpr(ComplexHeatmap::draw(Basal_cell_carcinoma_cor_rmse_perCellType)), NULL, 
                                           grid.grabExpr(ComplexHeatmap::draw(Gynecological_cancers_cor_rmse_perCellType)), NULL,
                                           grid.grabExpr(ComplexHeatmap::draw(HTAN_cor_rmse_perCellType)),
                                           label_size = 20, ncol = 1, nrow = 5, rel_heights = c(0.9, 0.05, 0.9, 0.05, 0.9), labels = c("A", "", "B", "", "C"))

pdf(paste0(fig_path, "Figure5_sup1.pdf"), w = 10, h = 5)
PBMC_benchmark_sup4
dev.off()

pdf(paste0(fig_path, "Figure5_sup2.pdf"), w = 10, h = 14)
cancer_benchmark_sup4
dev.off()

#################################################################
##                          Figure 5                           ##
#################################################################
summary_predictions <- read.table(paste0(benchmarking_output_path,"summary_atac_predictions.txt"), header = T, sep = "\t")

# consider only EPIC-ATAC
summary_predictions <- summary_predictions[which(summary_predictions$method %in% c("EPIC-ATAC")),]
unique(summary_predictions$bulk_name)
unique(summary_predictions$bulk_name_original)

summary_predictions$bulk_name[summary_predictions$bulk_name == "HTAN"] <- summary_predictions$bulk_name_original[summary_predictions$bulk_name == "HTAN"]
summary_predictions <- summary_predictions[which(!summary_predictions$bulk_name %in% c("PBMC_multiome_simulated")),]
summary_predictions$RMSE = summary_predictions$Correlation = 0


for(bulk in c("PBMC", "Basal cell carcinoma", "Gynecological cancers", "HNSCC","BRCA","CRC","CESC","SKCM","PDAC", "OV")){  
  
  
  cell_types <- unique(summary_predictions$cell_type[grepl(bulk, summary_predictions$bulk_name)])
  for(cell_type in cell_types){
    
    tmp <- summary_predictions[which(grepl(bulk, summary_predictions$bulk_name) &
                                       summary_predictions$cell_type == cell_type), ]
    if(anyNA(tmp$estimate)){
      summary_predictions$RMSE[which(grepl(bulk, summary_predictions$bulk_name) &
                                       summary_predictions$cell_type == cell_type)] <- NA
      
      summary_predictions$Correlation[which(grepl(bulk, summary_predictions$bulk_name) &
                                              summary_predictions$cell_type == cell_type)] <- NA
    } else{
      summary_predictions$RMSE[which(grepl(bulk, summary_predictions$bulk_name) &
                                       summary_predictions$cell_type == cell_type)] <- round(rmse(tmp$true_fraction, tmp$estimate), 3)
      
      cor_test <- cor.test(tmp$true_fraction, tmp$estimate)
      summary_predictions$Correlation[which(grepl(bulk, summary_predictions$bulk_name) &
                                              summary_predictions$cell_type == cell_type)] <- round(cor_test$estimate, 3)
    } 
    
    
  }
}
summary_predictions$bulk_name[summary_predictions$bulk_name %in% c("PBMC pseudobulk", "PBMC experiment")] <- "PBMC"
summary_predictions <- unique(summary_predictions[, c("bulk_name", "Correlation", "RMSE", "cell_type", "true_fraction")])

bulk_colors  <- c("Basal cell carcinoma" = '#4daf4aff', "Gynecological cancers" = '#3288bdff', "PBMC" = '#57C3F3',"HNSCC" = 'darkgreen', "BRCA" = '#E59CC4',
                  "CRC" = '#712820', "CESC" = '#8C549C', "SKCM" = 'black', "PDAC" = '#ff7f00ff', "OV" = "#d90017ff",
                  'yellow', '#712820', '#C1E6F3', '#6778AE', '#91D0BE', '#68A180', '#968175')[1:length(unique(summary_predictions$bulk_name))]

names(bulk_colors) <- unique(summary_predictions$bulk_name)

my_theme <- theme_bw() + theme(text = element_text(size = 12, face = "bold"),
                               axis.text = element_text(size = 12, face = "bold"),
                               legend.title = element_text(size = 12, face = "bold"),
                               legend.text = element_text(size = 12),
                               axis.title = element_text(size = 12, face = "bold"),
                               strip.text = element_text(size = 12, face = "bold"),
                               axis.text.x = element_text(angle=45, hjust=1, size = 12))


cell_type_boxplots_cor <- ggpubr::ggboxplot(summary_predictions, x = "cell_type", y = "Correlation")+ geom_jitter(position = "dodge", size = 2, aes(color = bulk_name)) +
  scale_color_manual(values = bulk_colors) +
  labs(x = "", y = "Correlation", color = "") +
  my_theme 


cell_type_boxplots_rmse <- ggpubr::ggboxplot(summary_predictions, x = "cell_type", y = "RMSE")+ geom_jitter(position = "dodge", size = 2, aes(color = bulk_name)) +
  scale_color_manual(values = bulk_colors) +
  labs(x = "", y = "RMSE", color = "") +
  my_theme 

cell_type_boxplots_prop <- ggpubr::ggboxplot(summary_predictions, x = "cell_type", y = "true_fraction")+ 
  geom_jitter(position = "dodge", size = 2, aes(color = bulk_name)) +
  scale_color_manual(values = bulk_colors) +
  labs(x = "", y = "True proportions", color = "") +
  my_theme 

pdf(paste0(fig_path, "Figure5_cell_type_boxplots.pdf"), w = 12, h = 7)

df_long <- tidyr::pivot_longer(summary_predictions, cols = c(Correlation, RMSE, true_fraction), names_to = "metric", values_to = "value")

ggpubr::ggboxplot(df_long, x = "cell_type", y = "value")+ geom_jitter(position = "dodge", size = 2, aes(color = bulk_name)) +
  scale_color_manual(values = bulk_colors) +
  facet_wrap(~ metric, scales = "free_y", strip.position = "left", nrow = 3, labeller = as_labeller(c("true_fraction" = "True proportions", "RMSE" = "RMSE", "Correlation" = "Correlation"))) +
  theme(strip.placement = "outside") + 
  labs(x = "", y = "", color = "") +
  my_theme 

dev.off()



##################################################################
##                        Time benchmark                        ##
##################################################################

time_files <- list.files(paste0(benchmarking_output_path, "time/"), 
                         pattern = "time", full.names = T)
summary_df=vector()
for(file in time_files){
  lines <- readLines(file)
  user <- gsub("user\t", "",lines[grepl("user\t", lines)])
  sys <- gsub("sys\t", "",lines[grepl("sys\t", lines)])
  
  get_time <- function(time_output){
    matches <- regmatches(time_output, regexec("([0-9]+)m([0-9.]+)s", time_output))
    minutes <- as.numeric(matches[[1]][2])
    seconds <- as.numeric(matches[[1]][3])
    total_seconds <- minutes * 60 + seconds
    return(total_seconds)
  }
  bulk_data_name <- strsplit(gsub("-time.txt","",basename(file)), "-")[[1]][1]
  method <- strsplit(gsub("-time.txt","",basename(file)), "-")[[1]][2]
  
  if(method != "EPIC"){
    # Load bulk 
    load(paste0(bulk_path, bulk_data_name, "_bulk.Rdata"))
    bulks_data_raw = bulks_data$bulk[, colnames(bulks_data$obs)[which(colnames(bulks_data$obs) != "cell_type")], drop = F]
    
    if(bulk_data_name %in% c("Gynecological_cancers","PBMC_multiome", "PBMC_multiome_simulated", "PBMC_experiment", "HNSCC","BRCA","CRC","CESC","SKCM","PDAC", "OV", "CEAD", "GBM", "UCEC")){
      genome_version = "hg38" 
    }else{
      genome_version = "hg19" 
    }   
    
    if(grepl(pattern = 'PBMC', bulk_data_name)){
      ref = EPICATAC::atacRef_PBMC
    }else{
      ref = EPICATAC::atacRef_TME
    }  
    
    execution_time <- system.time({
      # Match peaks from bulk and our reference 
      if (genome_version == "hg19") {
        lifted_matrix = EPICATAC:::run_liftOver(bulk_matrix = bulks_data_raw, from = genome_version)
      } else {
        lifted_matrix = bulks_data_raw
      }
      matched_bulks_data = EPICATAC:::match_peaks(bulk_matrix = lifted_matrix, profile_features = rownames(ref$refProfiles))
    })
    
    time_to_add <- execution_time[1] + execution_time[2]
  }else{
    time_to_add <- 0
  }
  
  
  
  summary_df <- rbind(summary_df, data.frame(dataset = bulk_data_name,
                                             method = method,
                                             CPU_time = get_time(user) + get_time(sys) + time_to_add
  )
  )
}

summary_df$method <- plyr::revalue(summary_df$method, c("EPIC" = "EPIC-ATAC",
                                                        "quantiSeq" = "quanTIseq",
                                                        "MPC_counter" = "MCPcounter",
                                                        "deconpeaker_OM" = "DeconPeaker_originalMarkers",
                                                        "DeconPeaker_ourRef" = "DeconPeaker-Custom",
                                                        "CIBERSORTx_ourRef" = "CIBERSORTx-Custom"))
pdf(paste0(fig_path, "time_benchmark.pdf"), w=6, h = 5)
ggplot(summary_df[which(summary_df$method != "DeconPeaker_originalMarkers"),], aes(x=method, y=CPU_time)) + 
  geom_boxplot() + ylab("CPU time (seconds)") + xlab("") + theme_bw() + 
  theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1))
dev.off()
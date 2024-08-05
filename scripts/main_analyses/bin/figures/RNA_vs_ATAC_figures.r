
# RNA vs ATAC -------------------------------------------------------------
benchmarking_output_path = as.character(commandArgs(TRUE)[1])
fig_path = as.character(commandArgs(TRUE)[2])
dir.create(fig_path)


bulk_names = c("PBMC_multiome_simulated", "HNSCC","BRCA","CRC","CESC","SKCM","PDAC", "OV") 


# Load RNA_deconvolution res
rna_files <- list.files(paste0(benchmarking_output_path, "RNA_deconvolution/"), pattern = "RNA_deconvolution_summary.txt", full.names = T)
rna_summary = vector()
for(file in rna_files){
  rna_summary = rbind(rna_summary, read.table(file, header = T, sep = "\t"))
  
}
rna_summary$bulk_name_original = rna_summary$bulk_name
rna_summary$bulk_name[which(rna_summary$bulk_name %in% c("HNSCC","BRCA","CRC","CESC","SKCM","PDAC", "OV", "GBM", "UCEC", "CEAD"))] = "HTAN"
rna_summary$cell_type <- plyr::revalue(rna_summary$cell_type, c("Other" = "Uncharacterized", "Myeloid" = "DCs + Macrophages",
                                                                "Non_Naive_CD4_Tcells" = "Memory + Helper CD4",
                                                                "Naive_CD4_Tcells" = "Naive CD4",
                                                                "Non_Naive_CD8_Tcells" = "Non Naive CD8",
                                                                "Naive_CD8_Tcells" = "Naive CD8",
                                                                "NK" = "NKcells",
                                                                "DCs" = "Dendritic"))
ga_summary = read.table(paste0(benchmarking_output_path, "summary_GA_predictions.txt"), header = T, sep = "\t")


modality_comp = data.frame(modality = vector(), method = vector(), sample = vector(), cor = vector(), rmse = vector(), bulk = vector())
method ="EPIC-ATAC"
  
summary_predictions = read.table(paste0(benchmarking_output_path, "summary_atac_predictions.txt"), header = T, sep = "\t")

atac_predictions = summary_predictions[which(summary_predictions$bulk_name_original %in% bulk_names), ]

tmp_atac = atac_predictions[grepl(method, atac_predictions$method),]
tmp_rna = rna_summary[grepl(unlist(strsplit(method, "_|-"))[[1]], rna_summary$method),]
tmp_ga = ga_summary[grepl(unlist(strsplit(method, "_|-"))[[1]], ga_summary$method),]

tmp_atac$id = paste0(tmp_atac$bulk_name_original, "_", tmp_atac$sample, "_", tmp_atac$cell_type)
colnames(tmp_atac)[c(3,5)] = c("estimate_ATAC", "true_fraction_ATAC")
tmp_rna$id = paste0(tmp_rna$bulk_name_original, "_", tmp_rna$sample, "_", tmp_rna$cell_type)
colnames(tmp_rna)[c(3,5)] = c("estimate_RNA", "true_fraction_RNA")
tmp_ga$id = paste0(tmp_ga$bulk_name_original, "_", tmp_ga$sample, "_", tmp_ga$cell_type)
colnames(tmp_ga)[c(3,5)] = c("estimate_GA", "true_fraction_GA")

merged_predictions = merge(tmp_atac, tmp_rna[,c("estimate_RNA", "true_fraction_RNA", "id")] )
merged_predictions = merge(merged_predictions, tmp_ga[,c("estimate_GA", "true_fraction_GA", "id")] )

for(bulk in unique(merged_predictions$bulk_name_original)){
  for(sample in unique(merged_predictions$sample[merged_predictions$bulk_name_original == bulk])){
    print(sample)
    merged_predictions_subset = merged_predictions[which(merged_predictions$sample == sample & merged_predictions$bulk_name_original == bulk),]
    
    # add atac res
    modality_comp = rbind(modality_comp, data.frame(modality = "ATAC",
                                                    method = method,
                                                    sample = sample,
                                                    cor = round(cor.test(merged_predictions_subset$true_fraction_ATAC, merged_predictions_subset$estimate_ATAC)$estimate,3),
                                                    rmse = round(Metrics::rmse(merged_predictions_subset$true_fraction_ATAC, merged_predictions_subset$estimate_ATAC),3),
                                                    bulk = bulk))
    
    # add rna res
    modality_comp = rbind(modality_comp, data.frame(modality = "RNA",
                                                    method = method,
                                                    sample = sample,
                                                    cor = round(cor.test(merged_predictions_subset$true_fraction_RNA, merged_predictions_subset$estimate_RNA)$estimate,3),
                                                    rmse = round(Metrics::rmse(merged_predictions_subset$true_fraction_RNA, merged_predictions_subset$estimate_RNA),3),
                                                    bulk = bulk))
    
    # add GA res
    modality_comp = rbind(modality_comp, data.frame(modality = "GA",
                                                    method = method,
                                                    sample = sample,
                                                    cor = round(cor.test(merged_predictions_subset$true_fraction_GA, merged_predictions_subset$estimate_GA)$estimate,3),
                                                    rmse = round(Metrics::rmse(merged_predictions_subset$true_fraction_GA, merged_predictions_subset$estimate_GA),3),
                                                    bulk = bulk))
    
  }
}



modality_comp$bulk[which(modality_comp$bulk != "PBMC_multiome_simulated")] <- "HTAN"

for(bulk_name in c("PBMC_multiome_simulated", "HTAN")){
  rna_atac_comp_subset = modality_comp[which(modality_comp$bulk == bulk_name), ]
  rna_atac_comp_subset$method = plyr::revalue(rna_atac_comp_subset$method, c("EPIC-ATAC" = "EPIC", 
                                                                             "DeconPeaker_ourMarkers" = "DeconPeaker", 
                                                                             "CIBERSORTx_ourMarkers" = "CIBERSORTx"))
  
  library(ggplot2)
  size = 16
  my_comparisons <- list( c("ATAC", "RNA"),c("ATAC", "GA"),c("GA", "RNA") )
  rna_atac_comp_subset$method <- factor(rna_atac_comp_subset$method, levels = c("EPIC")) 
  rna_atac_comp_subset <- rna_atac_comp_subset[rna_atac_comp_subset$method == "EPIC", ] 
  cor_p <- ggpubr::ggboxplot(rna_atac_comp_subset, x = "modality", y = "cor",
                             fill = "grey",  facet.by = "method", nrow = 1)+ 
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.14))) + 
    ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.format",paired = F,  vjust = -0.1, size = 4) + labs(x = "", y = "Correlation")+
    theme(legend.position = "none") + theme_bw() + theme(text = element_text(size = size, face = "bold"),
                                                         axis.text = element_text(size = size, face = "bold"),
                                                         legend.title = element_text(size = size, face = "bold"),
                                                         legend.text = element_text(size = size),
                                                         axis.title = element_text(size = size, face = "bold"),
                                                         strip.text = element_text(size = size, face = "bold"),
                                                         axis.text.x = element_text(angle=45, hjust=1, size = size)) 
  
  cor_p
  
  
  rmse_p <- ggpubr::ggboxplot(rna_atac_comp_subset, x = "modality", y = "rmse",
                              fill = "grey",  facet.by = "method", nrow = 1)+ 
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.14))) + 
    ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.format",paired = F,  vjust = 0.1, size = 4) + labs(x = "", y = "RMSE")+
    theme(legend.position = "none") + theme_bw() + theme(text = element_text(size = size, face = "bold"),
                                                         axis.text = element_text(size = size, face = "bold"),
                                                         legend.title = element_text(size = size, face = "bold"),
                                                         legend.text = element_text(size = size),
                                                         axis.title = element_text(size = size, face = "bold"),
                                                         strip.text = element_text(size = size, face = "bold"),
                                                         axis.text.x = element_text(angle=45, hjust=1, size = size)) 
  
  rmse_p
  
  
  pdf(paste0(fig_path, bulk_name, "_Figure7.pdf"), w = 7, h = 4)
  print(cowplot::plot_grid(cor_p , NULL, rmse_p,
                           label_size = 20, ncol = 3, rel_widths = c(0.9, 0.05, 0.9)))
  dev.off()
}





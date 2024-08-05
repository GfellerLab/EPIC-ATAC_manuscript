# Libraries ---------------------------------------------------------------
library(readxl)
library(ggplot2)
library(Metrics)

# Functions ---------------------------------------------------------------
get_pred_truth = function(deconv_res, bulk_proportions, cells_matching, method = "EPIC"){
  head(deconv_res)
  
  deconv_comp = data.frame(sample = vector(),
                           obs = vector(),
                           cell = vector(),
                           pred = vector())
  
  for(cell_type in unique(cells_matching$Grouping_name_bulk)){
    
    matching_rows = cells_matching[which(cells_matching[, "Grouping_name_bulk"] == cell_type), , drop = F] 
    
    prediction_columns = unique(matching_rows[, paste0("Cell_types_pred_",method)])
    prediction_columns = unique(as.vector(unlist(sapply(prediction_columns, function(i) unlist(strsplit(i,","))))))
    
    if(!"NA" %in% prediction_columns){
      cell_type_bulk_proportions = rowSums(bulk_proportions[, matching_rows[, "Cell_type_bulk"], drop = F])
      cell_type_predictions = rowSums(deconv_res[, prediction_columns[which(prediction_columns %in% colnames(deconv_res))] , drop =F])
    } else{
      cell_type_bulk_proportions = rowSums(bulk_proportions[, matching_rows[, "Cell_type_bulk"], drop = F])
      cell_type_predictions = rep(0, nrow(deconv_res))
    } 
    
    deconv_comp = rbind(
      deconv_comp,
      data.frame(
        sample = rownames(bulk_proportions),
        obs = cell_type_bulk_proportions,
        cell = rep(cell_type, nrow(bulk_proportions)),
        pred = cell_type_predictions
      )
    )
    
  } 
  
  considered_cell_types = unique(cells_matching[, paste0("Cell_types_pred_",method)])
  considered_cell_types = unique(as.vector(unlist(sapply(considered_cell_types, function(i) unlist(strsplit(i,","))))))
  
  cell_types_absent_in_bulk = colnames(deconv_res)[!colnames(deconv_res) %in% considered_cell_types]
  if (length(cell_types_absent_in_bulk) != 0) {
    for (cell_type in cell_types_absent_in_bulk) {
      deconv_comp = rbind(
        deconv_comp,
        data.frame(
          sample = rownames(bulk_proportions),
          obs = rep(0, nrow(deconv_res)),
          cell = rep(cell_type, nrow(bulk_proportions)),
          pred = deconv_res[, cell_type]
        )
      )
    }
  } 
  
  return(deconv_comp)
}

# Parameters --------------------------------------------------------------
benchmarking_output_path = as.character(commandArgs(TRUE)[1])
ga_res_path = as.character(commandArgs(TRUE)[2])
path_cell_matching = as.character(commandArgs(TRUE)[3])

# Gather predictions for each test dataset and each method ----------------

bulk_names = c("PBMC_multiome_simulated", "HNSCC","BRCA","CRC","CESC","SKCM","PDAC", "OV") 

summary_predictions = data.frame(cell_type = vector(),
                                 sample = vector(),
                                 estimate = vector(),
                                 method = vector(),
                                 true_fraction = vector(),
                                 bulk_name = vector())      
for(bulk in bulk_names){
  
  # Load true proportions ---------------------------------------------------
  bulk_proportions = read.table(paste0(benchmarking_output_path, "intermediate_files/", bulk, "/", bulk, "_bulk_proportions.txt"), header = T)
  rownames(bulk_proportions) = bulk_proportions$cell_type
  bulk_proportions = bulk_proportions[, c(1:(ncol(bulk_proportions) - 2)), drop = F]
  head(bulk_proportions)
  bulk_proportions = as.data.frame(t(bulk_proportions))
  
  # Read cells matching -----------------------------------------------------
  cells_matching = as.data.frame(read_excel(paste0(path_cell_matching)))
  head(cells_matching)
  cells_matching = cells_matching[which(cells_matching$Data_type == paste0("Bulk_", bulk)), ]
  
  # EPIC results ------------------------------------------------------------
  deconv_res = read.table(paste0(ga_res_path, bulk, "_EPIC_GA_predictions.txt"), header = T)
  deconv_res = deconv_res[, c(1:(ncol(deconv_res) - 1)), drop = F]
  
  deconv_comp = get_pred_truth(deconv_res, bulk_proportions, cells_matching, method = "withOtherCells")
  summary_predictions = rbind(summary_predictions, data.frame(cell_type = deconv_comp$cell,
                                                              sample = deconv_comp$sample,
                                                              estimate = deconv_comp$pred,
                                                              method = "EPIC_GA",
                                                              true_fraction = deconv_comp$obs,
                                                              bulk_name = rep(bulk, nrow(deconv_comp))))   
  
} 



summary_predictions$bulk_name_original = summary_predictions$bulk_name

summary_predictions$bulk_name[which(summary_predictions$bulk_name %in% c("HNSCC","BRCA","CRC","CESC","SKCM","PDAC", "OV", "GBM", "UCEC", "CEAD"))] = "HTAN"
summary_predictions$cell_type <- plyr::revalue(summary_predictions$cell_type, c("Other" = "Uncharacterized", "Myeloid" = "DCs + Macrophages", 
                                                                                "Non_Naive_CD4_Tcells" = "Memory + Helper CD4",
                                                                                "Naive_CD4_Tcells" = "Naive CD4", 
                                                                                "Non_Naive_CD8_Tcells" = "Non Naive CD8", 
                                                                                "Naive_CD8_Tcells" = "Naive CD8",
                                                                                "NK" = "NKcells",
                                                                                "DCs" = "Dendritic"))

# Save the predictions summary --------------------------------
write.table(summary_predictions, file = paste0("summary_GA_predictions.txt"), quote = F, row.names = F, col.names = T, sep = "\t")


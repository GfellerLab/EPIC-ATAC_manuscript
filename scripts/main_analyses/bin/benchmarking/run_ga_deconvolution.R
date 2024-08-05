library(EPIC)
library(plyr)

# Parameters --------------------------------------------------------------
bulk_path = as.character(commandArgs(TRUE)[1])
bulk_data_name = gsub("_bulk_GA.rdata", "", basename(bulk_path))

# Load counts -------------------------------------------------------------
load(bulk_path)
input_matrix = bulks_data_ga$bulk[, colnames(bulks_data_ga$obs)[which(colnames(bulks_data_ga$obs) != "cell_type")], drop = F]

epic_new_res = EPIC(
  bulk = input_matrix,
  reference = ifelse(grepl("PBMC", bulk_data_name), "BRef", "TRef"),
  withOtherCells = ifelse(grepl("PBMC", bulk_data_name), F, T)
)

epic_output = as.data.frame(epic_new_res$cellFractions)
epic_output$samples = rownames(epic_output)
colnames(epic_output) <- plyr::revalue(colnames(epic_output), c("CAFs" = "Fibroblasts", "NKcells" = "NK"))

write.table(epic_output, file = paste0(bulk_data_name, "_EPIC_GA_predictions.txt"), quote = F, row.names = F, col.names = T, sep = "\t")








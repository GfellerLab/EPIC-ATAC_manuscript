

# Parameters --------------------------------------------------------------

profile_path <- as.character(commandArgs(TRUE)[1])

# Save the reference profiles with major cell-types ---------------------------------------------
load(paste0(profile_path, "PBMC_profile_noSubtypes.Rdata"))
PBMC_ATAC_Ref_Major <- list()
PBMC_ATAC_Ref_Major[["refProfiles"]] <- profiles_list$refProfiles
PBMC_ATAC_Ref_Major[["refProfiles.var"]] <- profiles_list$refProfiles.var

PBMC_ATAC_Ref_Major[["sigPeaks"]] <- profiles_list$marker_peaks

markers_list <- sapply(colnames(PBMC_ATAC_Ref_Major[["refProfiles"]]), function(i){ gsub(":","-", markers_df$peak_id)[which(markers_df$cell_type == i & 
                                                                                                               gsub(":","-", markers_df$peak_id) %in% PBMC_ATAC_Ref_Major[["sigPeaks"]])]  })
PBMC_ATAC_Ref_Major[["markers"]] <- markers_list

load(paste0(profile_path, "TME_profile_noSubtypes.Rdata"))
Tumor_ATAC_Ref_Major <- list()
Tumor_ATAC_Ref_Major[["refProfiles"]] <- profiles_list$refProfiles
Tumor_ATAC_Ref_Major[["refProfiles.var"]] <- profiles_list$refProfiles.var
Tumor_ATAC_Ref_Major[["sigPeaks"]] <- profiles_list$marker_peaks

markers_list <- sapply(colnames(Tumor_ATAC_Ref_Major[["refProfiles"]]), function(i){ gsub(":","-", markers_df$peak_id)[which(markers_df$cell_type == i & 
                                                                                                                gsub(":","-", markers_df$peak_id) %in% Tumor_ATAC_Ref_Major[["sigPeaks"]])]  })
Tumor_ATAC_Ref_Major[["markers"]] <- markers_list

save(PBMC_ATAC_Ref_Major, file = paste0("PBMC_ATAC_Ref_Major.rda"))
save(Tumor_ATAC_Ref_Major, file = paste0("Tumor_ATAC_Ref_Major.rda"))

# Save the reference profiles with subtypes ---------------------------------------------

load(paste0(profile_path, "PBMC_profile_withSubtypes.Rdata"))
PBMC_ATAC_Ref_Major
PBMC_ATAC_Ref_Subtypes <- list()
PBMC_ATAC_Ref_Subtypes[["refProfiles"]] <- profiles_list$refProfiles
PBMC_ATAC_Ref_Subtypes[["refProfiles.var"]] <- profiles_list$refProfiles.var

PBMC_ATAC_Ref_Subtypes[["sigPeaks"]] <- profiles_list$marker_peaks


markers_list <- sapply(colnames(PBMC_ATAC_Ref_Subtypes[["refProfiles"]]), function(i){ gsub(":","-", markers_df$peak_id)[which(markers_df$cell_type == i &
                                                                                                                  gsub(":","-", markers_df$peak_id) %in% PBMC_ATAC_Ref_Subtypes[["sigPeaks"]])]  })
PBMC_ATAC_Ref_Subtypes[["markers"]] <- markers_list

load(paste0(profile_path, "TME_profile_withSubtypes.Rdata"))
Tumor_ATAC_Ref_Subtypes <- list()
Tumor_ATAC_Ref_Subtypes[["refProfiles"]] <- profiles_list$refProfiles
Tumor_ATAC_Ref_Subtypes[["refProfiles.var"]] <- profiles_list$refProfiles.var


Tumor_ATAC_Ref_Subtypes[["sigPeaks"]] <- profiles_list$marker_peaks
markers_list <- sapply(colnames(Tumor_ATAC_Ref_Subtypes[["refProfiles"]]), function(i){ gsub(":","-", markers_df$peak_id)[which(markers_df$cell_type == i &
                                                                                                                   gsub(":","-", markers_df$peak_id) %in% Tumor_ATAC_Ref_Subtypes[["sigPeaks"]])]  })
Tumor_ATAC_Ref_Subtypes[["markers"]] <- markers_list

save(PBMC_ATAC_Ref_Subtypes, file = paste0("PBMC_ATAC_Ref_Subtypes.rda"))
save(Tumor_ATAC_Ref_Subtypes, file = paste0("Tumor_ATAC_Ref_Subtypes.rda"))


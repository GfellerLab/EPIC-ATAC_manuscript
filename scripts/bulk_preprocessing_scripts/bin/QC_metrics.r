##################################################################
##                       Extract QC stats                       ##
##################################################################


# Parameters --------------------------------------------------------------

path_to_preprocessed_data <- as.character(commandArgs(TRUE)[1]) 

# Gather all stats --------------------------------------------------------

stats_files <- list.files(path_to_preprocessed_data, pattern = "stats.tsv", recursive = T, full.names = T)
summary_df <- data.frame(sample=vector(),NRF=vector(),PBC1=vector(),PBC2=vector(),TSS_score=vector(),FRiP=vector(),Peak_count=vector(),Alignment_rate=vector())
for(i in stats_files){
  sample <- unlist(strsplit(i,"/"))
  sample <- sample[grepl("SRR|EGAN",sample)]
  stats_df <- read.table(i,header = F,sep="\t")
  summary_df <- rbind(summary_df,data.frame(sample=sample,
                                           NRF=as.numeric(stats_df$V2[which(stats_df$V1 == "NRF")]),
                                           PBC1=as.numeric(stats_df$V2[which(stats_df$V1 == "PBC1")]),
                                           PBC2=as.numeric(stats_df$V2[which(stats_df$V1 == "PBC2")]),
                                           TSS_score=as.numeric(stats_df$V2[which(stats_df$V1 == "TSS_score")]),
                                           FRiP=as.numeric(stats_df$V2[which(stats_df$V1 == "FRiP")]),
                                           Peak_count=as.numeric(stats_df$V2[which(stats_df$V1 == "Peak_count")]),
                                           Alignment_rate=as.numeric(stats_df$V2[which(stats_df$V1 == "Alignment_rate")])))
}
head(summary_df)

write.table(summary_df, file = "summary_QC_stats.txt", quote = F, col.names = T, row.names = F, sep = "\t")


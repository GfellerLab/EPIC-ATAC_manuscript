
# Libraries ---------------------------------------------------------------


res_path=as.character(commandArgs(TRUE)[1]) 
sample_pattern=as.character(commandArgs(TRUE)[2]) 
annotation_path = as.character(commandArgs(TRUE)[3]) 

stats_files=list.files(res_path,pattern = "stats.tsv",recursive = T,full.names = T)
samples_annot = read.table(annotation_path,header = T,sep="\t")

summary_df = data.frame(sample=vector(),NRF=vector(),PBC1=vector(),PBC2=vector(),TSS_score=vector(),FRiP=vector(),Peak_count=vector(),Alignment_rate=vector())
for(i in stats_files){
  sample = unlist(strsplit(i,"/"))
  sample = sample[grepl(sample_pattern,sample)]
  stats_df = read.table(i,header = F,sep="\t")
  summary_df = rbind(summary_df,data.frame(sample=sample,
                                           NRF=as.numeric(stats_df$V2[which(stats_df$V1 == "NRF")]),
                                           PBC1=as.numeric(stats_df$V2[which(stats_df$V1 == "PBC1")]),
                                           PBC2=as.numeric(stats_df$V2[which(stats_df$V1 == "PBC2")]),
                                           TSS_score=as.numeric(stats_df$V2[which(stats_df$V1 == "TSS_score")]),
                                           FRiP=as.numeric(stats_df$V2[which(stats_df$V1 == "FRiP")]),
                                           Peak_count=as.numeric(stats_df$V2[which(stats_df$V1 == "Peak_count")]),
                                           Alignment_rate=as.numeric(stats_df$V2[which(stats_df$V1 == "Alignment_rate")])))
}
head(summary_df)
head(samples_annot)
merged_data = merge(summary_df,samples_annot,by.x = "sample",by.y = "SRA_ID")


table(merged_data$cell_type2)
filtered_samples = merged_data[which(merged_data$TSS_score >= 5),]

table(filtered_samples$cell_type2)
dim(filtered_samples)
merged_data$to_keep = F
merged_data$to_keep[which(merged_data$sample %in% filtered_samples$sample)] = T

write.table(merged_data,file='summary_stats.txt',quote = F,row.names = F,col.names = T,sep="\t")


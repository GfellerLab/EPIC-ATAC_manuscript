
# Determine overlap between markers identified in the different subsampling 
markers_all = vector()
for(i in 1:10){
  markers_df = read.table(paste0("markers_limma-cv",i,".txt"), header = T, sep = '\t')

  markers_df$subsampling = i
  markers_all = rbind(markers_all,markers_df)
} 

summary_frequencies = data.frame()
redundant_peaks = data.frame(cell_type = vector(), peak_id = vector())
for(cell_type in unique(markers_all$cell_type)){
  markers_subset = markers_all[which(markers_all$cell_type == cell_type),]
  peaks_freq = as.data.frame( table(markers_subset$peak_id))
  peaks_freq$Freq2 = peaks_freq$Freq
  peaks_freq$Freq2[which(peaks_freq$Freq < 3)] = "< 3"
  peaks_freq$Freq2[which(peaks_freq$Freq >= 7)] = ">= 7"
  peaks_freq$Freq2[which(peaks_freq$Freq >= 3 & peaks_freq$Freq < 7)] = ">=3 & <7"
  
  redundant_peaks = rbind(redundant_peaks,data.frame(cell_type = rep(cell_type, length(peaks_freq$Var1[which(peaks_freq$Freq2 != "< 3")])),
                                                     peak_id = as.vector(peaks_freq$Var1[which(peaks_freq$Freq2 != "< 3")])) )
  recurrence_freq = as.data.frame( table(peaks_freq$Freq2))
  recurrence_freq$cell_type = cell_type
  summary_frequencies = rbind(summary_frequencies, recurrence_freq)
} 

library(ggplot2)
pdf(paste0("markers_freq_limma.pdf"))
ggplot(summary_frequencies, aes(fill=Var1, y=Freq, x=cell_type)) + 
  geom_bar(position="stack", stat="identity") +
  geom_text(aes(label=Freq),position = position_stack(0.5))
dev.off()

write.table(redundant_peaks, file = paste0("markers_limma-consensus.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

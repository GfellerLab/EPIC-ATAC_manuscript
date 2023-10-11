# Libraries ---------------------------------------------------------------
library(dplyr)
library(GenomicRanges)

# Parameters --------------------------------------------------------------
counts_path = as.character(commandArgs(TRUE)[1])
samples_metadata_path = as.character(commandArgs(TRUE)[2]) 
QC_path = as.character(commandArgs(TRUE)[3]) 

# Load the count matrix ------------------------------
counts = read.table(counts_path, skip = 1, header = T)
rownames_save = counts$Geneid
counts = counts[, c(7:ncol(counts))]
rownames(counts) = rownames_save
colnames(counts) = sapply(colnames(counts), function(i) unlist(strsplit(i, '_sort_dedup.bam'))[[1]])
dim(counts)

# Load metadata ------------------------------

# Load samples metadata
metadata = read.table(samples_metadata_path, header = T, sep = "\t")
# Load QC stats and keep the good quality samples
QCs_summary = read.table(QC_path, header = T, sep = "\t")
metadata = merge(QCs_summary[,c(1:8,15)], metadata, by.y = "SRA_ID", by.x = "sample")

metadata = metadata[which(metadata$to_keep == T),] 
sample_counts = metadata %>% group_by(groups, study) %>% count()

counts = counts[,which(colnames(counts) %in% metadata$sample)]
counts = counts[grep("chrY|chrUn|random",rownames(counts),invert = T),]
dim(counts)

rownames(metadata) = metadata$sample
metadata = metadata[colnames(counts),] 

peaks_gr = as(rownames(counts),"GRanges")
region.length = peaks_gr@ranges@width
names(region.length) = rownames(counts)
x <-  counts  / region.length[rownames(counts)]
tpm.mat <- t( t(x) * 1e6 / colSums(x) )
head(colSums(tpm.mat))


# Save raw counts ---------------------------------------------------------
write.table(counts, file = "raw_counts.txt", quote = F, row.names = T, col.names = T, sep = "\t") 

# Save tpm counts ---------------------------------------------------------
write.table(tpm.mat, file = "tpm.txt", quote = F, row.names = T, col.names = T, sep = "\t") 

# Save ordered metadata ---------------------------------------------------
write.table(metadata, file = "metadata.txt", quote = F, row.names = F, col.names = T, sep = "\t") 


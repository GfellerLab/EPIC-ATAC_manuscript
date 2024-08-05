
# Libraries ---------------------------------------------------------------
library(dplyr)
library(ChIPseeker)
cache_dir <- "/tmp/R_cache/BiocFileCache"
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
library(JASPAR2022)
library(chipenrich)
library(ggplot2)
set.seed(1234)

# Functions ---------------------------------------------------------------
run_chipSeeker_annotation = function(peaks_gr, genome_version = "hg38", output_dir = NULL){
  if(genome_version == "hg38"){
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  } else if(genome_version == "hg19"){
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  } 
  anno = ChIPseeker::annotatePeak(
    peaks_gr,
    TxDb = txdb,
    tssRegion=c(-1000, 1000),
    addFlankGeneInfo = TRUE,
    overlap = "all", #any nearest gene is reported regardless of the overlap with the TSS 
    annoDb = "org.Hs.eg.db"
  )
  
  # Save annotation 
  anno_df = as.data.frame(anno)
  head(anno_df)
  
  anno_df$annotation_type = anno_df$annotation
  anno_df$annotation_type[grep("Intron",anno_df$annotation_type)]="Intron"
  anno_df$annotation_type[grep("Exon",anno_df$annotation_type)]="Exon"
  anno_df$annotation_type[grep("Promoter",anno_df$annotation_type)]="Promoter"
  anno_df$annotation_type[grep("Downstream",anno_df$annotation_type)]="Downstream"
  
  table(anno_df$annotation_type)
  
  pdf(paste0(output_dir, "ChipSeeker_annoBar.pdf"), w=7, h=3)
  print(ChIPseeker::plotAnnoBar(anno))
  dev.off()
  
  pdf(paste0(output_dir, "ChipSeeker_dTSS.pdf"), w=7, h=3)
  print(ChIPseeker::plotDistToTSS(anno))
  dev.off()
  
  pdf(paste0(output_dir, "ChipSeeker_annoPie.pdf"))
  print(ChIPseeker::plotAnnoPie(anno))
  dev.off()
  
  return(anno_df)
} 

get_signac_annotation <- function(background_peaks, pfm, genome_version = "hg38"){
  
  # Generate signac object from the reference profiles ----------------------
  chrom_assay <- Signac::CreateChromatinAssay(
    counts =  data.frame(row.names = background_peaks, cell1 = rep(1,length(background_peaks)), cell2 = rep(1,length(background_peaks))), min.cells = 0, min.features = 0,
    sep = c("-", "-"),
    ranges = Signac::StringToGRanges(background_peaks)
  )
  
  signac_data <- Seurat::CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks"
  )
  
  if(genome_version == "hg38"){
    # add motif information
    signac_data <- Signac::AddMotifs(
      object = signac_data,
      genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
      pfm = pfm
    )
  } else{
    # add motif information
    signac_data <- Signac::AddMotifs(
      object = signac_data,
      genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
      pfm = pfm
    )
  } 

  return(signac_data)
} 

run_monaLisa <- function(peaks_gr, bin_annotation, pwm, genome_version, background, nb_cores = 3){
  peak_bins <- as.factor(bin_annotation)
  print(table(peak_bins))
  
  # Resize regions to avoid length biases
  peaks_gr <- IRanges:::trim(GenomicRanges::resize(peaks_gr,
                                                   width = median(rtracklayer::width(peaks_gr)),
                                                   fix = "center")
  )
  
  if (genome_version == "hg38") {
    bs <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  } else if (genome_version == "hg19") {
    bs <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  }
  lmrseqs <- Biostrings::getSeq(bs, peaks_gr)
  print(monaLisa::plotBinDiagnostics(seqs = lmrseqs, bins = peak_bins, aspect = "GCfrac"))
  print(monaLisa::plotBinDiagnostics(seqs = lmrseqs, bins = peak_bins, aspect = "dinucfreq"))
  
  se <- monaLisa::calcBinnedMotifEnrR(seqs = lmrseqs,
                                      bins = peak_bins,
                                      pwmL = pwm,
                                      BPPARAM = BiocParallel::MulticoreParam(nb_cores),
                                      test="binomial", 
                                      background = "genome",
                                      genome.regions = background,
                                      genome = bs, genome.oversample = 10)
  
  p_adj_mat <- SummarizedExperiment::assay(se, "negLog10Padj")
  print(dim(p_adj_mat))
  
  # select strongly enriched motifs
  sel <- apply(SummarizedExperiment::assay(se, "negLog10Padj"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > 1.3
  seSel <- se[sel, ]
  
  significance_matrix <- p_adj_mat[sel,] 
  
  # Remove duplicated motifs before renaming the rows 
  dup_motifs <- which(duplicated(seSel@elementMetadata@listData$motif.name))
  significance_matrix <- significance_matrix[-dup_motifs, ] 
  seSel <- seSel[-dup_motifs, ]
  rownames(significance_matrix) <- seSel@elementMetadata@listData$motif.name
  
  summary_monaLisa <- as.data.frame(data.table::data.table(target = rownames(significance_matrix), cell_type = colnames(significance_matrix), significance = c(significance_matrix)))
  return(list(seSel = seSel, significance_matrix = significance_matrix, summary_monaLisa = summary_monaLisa))
} 



# Parameters --------------------------------------------------------------
profile_path <- as.character(commandArgs(TRUE)[1]) 
remap_catalog <- as.character(commandArgs(TRUE)[2]) 
jaspar_database <- as.character(commandArgs(TRUE)[3]) 
set.seed(1234)

# Load our markers ------------------------------------------------------------
major_groups <- ifelse(grepl("withSubtypes", profile_path), "withSubtypes", "noSubtypes")
if(grepl("PBMC", profile_path)){
  output_path <- paste0("PBMC_", major_groups, "/")
  fig_name <- paste0("PBMC_", major_groups)
}else{
  output_path <- paste0("TME_", major_groups, "/")
  fig_name <- paste0("TME_", major_groups)
} 
dir.create(output_path)
load(profile_path)

markers <- profiles_list$marker_peaks

markers_df <- markers_df[which(markers_df$peak_id %in% markers & markers_df$cell_type %in% colnames(profiles_list$refProfiles)), ]
gr_markers <- Signac::StringToGRanges(markers_df$peak_id)
gr_markers$cell_type <- markers_df$cell_type
gr_markers$peak_id <- markers_df$peak_id

# Run ChipSeeker annotation -----------------------------------------------
anno_df <- run_chipSeeker_annotation(peaks_gr = gr_markers,
                                     genome_version = "hg38",
                                     output_dir = output_path)

write.table(anno_df, file = paste0(output_path, "chipSeeker_output.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

# Run Remap annotation ------------------------------------------------------
catalog_nr <- ReMapEnrich::bedToGranges(remap_catalog)
catalog_nr <- catalog_nr[grepl("macrophage|CD8|T-cell|nk|endothelial|monocyte|CD4|T-cell|Th1|CD4-pos|neutrophil|MV4-11|MOLM|NB4|SET-2|HL-60|OCI-AML|THP-1|Kasumi-1|B-cell|fibroblast",catalog_nr$id),] 

biotype_matching <- c("Macrophages" = "macrophage",
                      "CD8_Tcells" = "CD8|T-cell",
                      "Naive_CD8_Tcells" = "CD8|T-cell",
                      "Non_Naive_CD8_Tcells" = "CD8|T-cell",
                      "NK" = "nk",
                      "Endothelial" = "endothelial",
                      "Monocytes" = "monocyte",
                      "CD4_Tcells" = "CD4|T-cell|Th1|CD4-pos",
                      "Naive_CD4_Tcells" = "CD4|T-cell|Th1|CD4-pos",
                      "Non_Naive_CD4_Tcells" = "CD4|T-cell|Th1|CD4-pos",
                      "Tregs" = "CD4|T-cell|Th1|CD4-pos",
                      "Neutrophils" = "neutrophil",
                      "DCs" = "MV4-11|MOLM|NB4|SET-2|HL-60|OCI-AML|THP-1|Kasumi-1",
                      "Bcells" = "B-cell",
                      "Fibroblasts" = "fibroblast")

# Remap enrichment 
remap_enrichment_summary <- vector()
for(cell_type in unique(anno_df$cell_type)){
  print(cell_type)
  enrichment.df <- ReMapEnrich::enrichment(GenomicRanges::makeGRangesFromDataFrame(anno_df[which(anno_df$cell_type == cell_type),]), 
                                           catalog_nr, 
                                           universe = Signac::StringToGRanges(rownames(profiles_list$refProfiles)), 
                                           nCores = 3, shuffles = 10)
  enrichment.df$cell_type = cell_type
  enrichment.df <- enrichment.df[which(enrichment.df$q.value < 0.05), ] 
  
  remap_enrichment_summary <- rbind(remap_enrichment_summary, enrichment.df)
} 
remap_summary <- remap_enrichment_summary %>% tidyr::separate(col = category, into = c("target", "biotype"), sep = ":", remove = F)
write.table(remap_summary, file = paste0(output_path, "ReMapEnrich_output.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

remap_summary_restricted <- vector()
for(cell_type in unique(remap_summary$cell_type)){
  remap_summary_restricted <- rbind(remap_summary_restricted, remap_summary[which(remap_summary$cell_type == cell_type & grepl(biotype_matching[cell_type], remap_summary$biotype)),])
} 
write.table(remap_summary_restricted, file = paste0(output_path, "ReMapEnrich_restricted_output.txt"), quote = F, row.names = F, col.names = T, sep = "\t")


rm(catalog_nr)
gc()

# TF assignment -----------------------------------------------------------

# Get a list of motif position frequency matrices from the JASPAR database 
jasparDb <- new("JASPAR2022", db=paste0(jaspar_database))
pfm <- TFBSTools::getMatrixSet(
  x = jasparDb,
  opts = list(collection = "CORE", tax_group = 'vertebrates', species = "Homo sapiens", all_versions = FALSE)
)

signac_annotation <- get_signac_annotation(background_peaks = rownames(profiles_list$refProfiles), pfm = pfm,
                                           genome_version = "hg38")

enriched_motifs_signif <- list()
for(cell_type in unique(anno_df$cell_type)){
  print(cell_type)
  cell_type_peaks <- anno_df$peak_id[which(anno_df$cell_type == cell_type)]
  
  enriched.motifs <- Signac::FindMotifs(
    object = signac_annotation,
    features = cell_type_peaks
  )
  
  enriched_motifs_signif[[cell_type]] = as.data.frame(enriched.motifs[which(enriched.motifs$p.adjust < 0.05),])
  print(dim(enriched_motifs_signif[[cell_type]]))
  
  if(nrow(enriched_motifs_signif[[cell_type]]) != 0){
    motif_list = rownames(enriched_motifs_signif[[cell_type]])
    motif_names = unlist(signac_annotation@assays$peaks@motifs@motif.names[motif_list])
    
    enriched_motifs_signif[[cell_type]] = cbind(enriched_motifs_signif[[cell_type]],data.frame(cell_type = rep(cell_type,nrow(enriched_motifs_signif[[cell_type]]))))
    
  } 
  
}
enriched_motifs_signif_df <- do.call(rbind, enriched_motifs_signif)

write.table(enriched_motifs_signif_df, file = paste0(output_path, "SignacEnrichment_output.txt"), quote = F, row.names = F, col.names = T, sep = "\t")


# MonaLisa analysis -------------------------------------------------------
pfm <- TFBSTools::getMatrixSet(
  x = jasparDb,
  opts = list(collection = "CORE", tax_group = 'vertebrates', species = "Homo sapiens", matrixtype = "PWM", all_versions = FALSE)
)
monaLisa_res <- run_monaLisa(peaks_gr = gr_markers, 
                             bin_annotation = anno_df$cell_type, 
                             pwm  = pfm,
                             genome_version = "hg38", 
                             background = Signac::StringToGRanges(rownames(profiles_list$refProfiles)), 
                             nb_cores = 2)
monaLisa_summary <- monaLisa_res$summary_monaLisa
monaLisa_summary <- monaLisa_summary[which(monaLisa_summary$significance > 1.3),] 

write.table(monaLisa_summary, file = paste0(output_path, "MonaLisaEnrichment_output.txt"), quote = F, row.names = F, col.names = T, sep = "\t")


# Combined all enrichment results -----------------------------------------
enriched_motifs_signif_df$significance = -log10(enriched_motifs_signif_df$p.adjust)
signac_subset <- enriched_motifs_signif_df[,c("motif.name","cell_type","significance")]
colnames(signac_subset)[1]  <- c('target')
signac_subset$method <- "Signac"

remap_subset <- remap_summary_restricted[, c('target',"cell_type","q.significance")] 
colnames(remap_subset)  <- c('target',"cell_type","significance")
remap_subset$method <- "Remap"

monaLisa_summary$method <- "MonaLisa"

combined_enrichment_res <- rbind(signac_subset, remap_subset, monaLisa_summary)
write.table(combined_enrichment_res, file = paste0(output_path, "combined_enrichment_res.txt"), quote = F, row.names = F, col.names = T, sep = "\t")



# Run chipenrich ----------------------------------------------------------

geneset_param = "reactome"
locusDef = "nearest_tss"
enrichment_method = "chipenrich"

peak_annotation_df <- anno_df
anno_output = list()
anno_output_df = vector()
peaks_annotation_chipenrich = vector()
for(cell_type in unique(peak_annotation_df$cell_type)){
  print(cell_type)
  
  cell_type_markers <- peak_annotation_df[which(peak_annotation_df$cell_type == cell_type), "peak_id"]
  peaks_ids <- as.data.frame(do.call(rbind, strsplit(cell_type_markers, "[:-]")))
  rownames(peaks_ids) <- cell_type_markers
  colnames(peaks_ids) <- c("chr", "start", "end")
  
  results <- chipenrich::chipenrich(peaks = peaks_ids, genome = 'hg38', genesets = "GOBP" , locusdef = locusDef, 
                                    qc_plots = F, out_name = NULL, n_cores = 4)
  
  
  signif_pathways <- results$results[which(results$results$FDR < 0.05),]
  peaks_assignment <- results$peaks
  peaks_assignment$peak_id <- paste0(peaks_assignment$chr, "-", peaks_assignment$peak_start, "-", peaks_assignment$peak_end)
  peaks_assignment$cell_type = cell_type
  genes_names_matching <- unique(peaks_assignment[, c("gene_id", "gene_symbol")])
  rownames(genes_names_matching) <- genes_names_matching$gene_id 
  if(nrow(signif_pathways) != 0){
    signif_pathways$genes_symbols <- sapply(1:nrow(signif_pathways), function(i){
      entrez_ids <- unlist(strsplit(signif_pathways$Geneset.Peak.Genes[i], ", "))
      gene_symbols <- genes_names_matching[entrez_ids, "gene_symbol"] 
      return(paste(gene_symbols, collapse = ", "))
    })
    
    signif_pathways$peaks <- sapply(1:nrow(signif_pathways), function(i){
      entrez_ids <- unlist(strsplit(signif_pathways$Geneset.Peak.Genes[i], ", "))
      gene_symbols <- genes_names_matching[entrez_ids, "gene_symbol"] 
      peaks_id <- sapply(gene_symbols, function(g) paste(peaks_assignment[which(peaks_assignment$gene_symbol == g), "peak_id"], collapse = "|") )
      return(paste(paste0(gene_symbols, ":", peaks_id), collapse = ", "))
    })
  } 
  anno_output[[cell_type]] = cbind(signif_pathways, data.frame(cell_type = rep(cell_type, nrow(signif_pathways))) )
  
  anno_output_df = rbind(anno_output_df, anno_output[[cell_type]])
  
  peaks_annotation_chipenrich = rbind(peaks_annotation_chipenrich, peaks_assignment)
  
  gc()
}
write.table(anno_output_df,
            file = paste0(output_path, "GO_enrichment.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

write.table(peaks_annotation_chipenrich,
            file = paste0(output_path, "Nearest_gene.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

go_res_top = vector()
for(cell_type in unique(anno_output_df$cell_type)){
  subset = anno_output_df[which(anno_output_df$cell_type == cell_type),] 
  go_res_top = rbind(go_res_top, subset[1:min(20, nrow(subset)),] )
  
} 
go_res_top$Significance <- -log10(go_res_top$FDR)
go_res_top$Description <- factor(go_res_top$Description, levels=unique((go_res_top$Description)[order(go_res_top$cell_type)]))

# pdf(paste0("fig2_pathways_",fig_name,".pdf"), w = 15, h = 14)
# p <- ggplot(go_res_top , aes(cell_type, Description)) + geom_tile(aes(fill = Significance), colour = "white") + scale_fill_gradient(low = "yellow", high = "red")
# text_size=15
# print(p + theme_bw() + 
#         theme(axis.text.x = element_text(angle=45, hjust=1, size = text_size),
#               axis.text.y = element_text( size = text_size),
#               axis.title.y = element_text( size = text_size),
#               text = element_text( size = text_size)))
# dev.off()


library(PEPATACr)
SRA_table_path <- as.character(commandArgs(TRUE)[1])
pep <- as.character(commandArgs(TRUE)[2])
chrom_sizes_path = as.character(commandArgs(TRUE)[3])
min_samples <- as.numeric(commandArgs(TRUE)[4]) #2
stats_path <- as.character(commandArgs(TRUE)[5])
min_score   <- as.numeric(commandArgs(TRUE)[6])
group_id = "groups"
organism <- as.character(commandArgs(TRUE)[7])

# functions ----------------------------------------------------------------

consensusPeaks_custom <- function(sample_table, summary_dir, results_subdir, assets,
                                  min_score=5, min_olap=1) {    
  
  # Produce summary output directory (if needed)
  dir.create(summary_dir, showWarnings = FALSE)
  
  setDT(sample_table)[assets[asset == 'chrom_sizes', ],
                      c_path := i.path, on = 'sample_name']
  
  # generate paths to peak files
  sample_table[,peak_files:=.((file.path(
    results_subdir,
    sample_table$sample_name,
    paste0("peak_calling_", sample_table$genome),
    paste0(sample_table$sample_name,
           "_peaks_normalized.narrowPeak"))))]
  
  #Only keep samples with valid peak files
  file_list   <- sample_table$peak_files
  file_exists <- character()
  for (i in 1:length(file_list)) {
    if(file.exists(file.path(file_list[i]))) {
      file_exists <- append(file_exists, file.path(file_list[i]))
    }
  }
  files <- data.table(peak_files=file_exists)
  consensus_peak_files = list()
  if (nrow(files) == 0) {
    return(consensus_peak_files)
  }
  
  sample_table <- unique(
    sample_table[sample_table$peak_files %in% files$peak_files,])
  
  # Need to group by genome, then create a consensus list by genome!
  st_list = splitDataTable(sample_table, "study")
  
  study_list=list()
  for (g in names(st_list)) {
    if (nrow(st_list[[g]]) == 1) {
      err_msg = paste0("Found only a single valid peak file for ",
                       g, ".")
      warning(err_msg)
      next
    }
    if (nrow(st_list[[g]]) == 0) {
      warning("Unable to find any valid peak files.")
      warning("Confirm peak files exist for your samples.")
      next
    }
    c_path <- unique(sample_table[study == g, c_path])
    if (file.exists(c_path)) {
      c_size <- fread(c_path)
      colnames(c_size) <- c("chr", "size")
    } else {
      warning("Unable to load the chromosome sizes file.")
      warning(paste0("Confirm that ", c_path,
                     " is present before continuing."))
      final <- NULL
    }
    message(paste0("Calculating ", g, " consensus peak set from ",
                   nrow(st_list[[g]]), " samples..."))
    min_samples_1 = round(nrow(st_list[[g]])/2)
    if(nrow(st_list[[g]]) == 2){
      min_samples_1 = 2
    }
    study_list[[g]] <- collapsePeaks(st_list[[g]], c_size,
                                     min_samples=min_samples_1, min_score, min_olap)
    study_list[[g]]$file=g
    
    
  }
  
  min_samples=length(study_list)
  final=collapsePeaks_custom(study_list, c_size,
                             min_samples, min_score, min_olap)
  
  if (!is.null(final)) {
    # save consensus peak set
    file_name   <- paste0("_hg38_consensusPeaks.narrowPeak")
    output_file <- file.path(summary_dir,
                             paste0(project_name, file_name))
    fwrite(final, output_file, sep="\t", col.names=FALSE)
    consensus_peak_files <- c(consensus_peak_files, output_file)
    rm(final)
    invisible(gc())
  } else {
    warning("Unable to produce a consensus peak file.")
    warning("Check that individual peak files exist for your samples.")
  }
  
  return(consensus_peak_files)
}

collapsePeaks_custom <- function(study_list, chr_sizes,
                                 min_samples=2, min_score=5, min_olap=1) {
  # create combined peaks
  peaks <- rbindlist(study_list, idcol="file")
  if (ncol(peaks) == 7) {
    colnames(peaks) <- c("file", "chr", "start", "end",
                         "name", "score", "strand")
  } else if (ncol(peaks) == 12) {
    colnames(peaks) <- c('file',"chr", "start", "end",
                         "name", "score", "strand",
                         "signalValue", "pValue", "qValue", "peak","file")
  } else {
    warning(paste0("Peak files did not contain a recognizable number", 
                   " of columns (", ncol(peaks), ")"))
    rm(peaks)
    final <- data.table(chr=character(),
                        start=integer(),
                        end=integer(),
                        name=character(),
                        score=numeric(),
                        strand=character(),
                        signalValue=numeric(),
                        pValue=numeric(),
                        qValue=numeric(),
                        peak=integer())
    return(final)
  }
  setkey(peaks, chr, start, end)
  # keep highest scored peaks
  # split by chromosome to minimize memory requirements
  peaks_by_chr   <- split(peaks, peaks$chr)
  hit_aggregator <- function(x) {
    #message(paste0("x: ", unique(x$chr)))  # DEBUG
    peaksGR <- makeGRangesFromDataFrame(x, keep.extra.columns=FALSE)
    hitsGR  <- suppressWarnings(
      findOverlaps(peaksGR, peaksGR,
                   ignore.strand=TRUE, minoverlap=min_olap))
    hits    <- data.table::data.table(xid=queryHits(hitsGR),
                                      yid=subjectHits(hitsGR))
    setkey(hits, xid)
    scores  <- data.table(index=rep(1:nrow(x)), score=x$score)
    setkey(scores, index)
    out     <- hits[scores, nomatch=0]
    keep    <- out[out[,.I[which.max(score)],by=yid]$V1]
    indices <- unique(keep$xid)
    reduced <- x[indices,]
    reduced[start < 0, start := 0]
    return(reduced)
  }
  final <- rbindlist(lapply(peaks_by_chr, hit_aggregator))
  
  # can't extend past chromosome
  for (i in nrow(chr_sizes)) {
    final[chr == chr_sizes$chr[i] & end > chr_sizes$size[i],
          end := chr_sizes$size[i]]
  }
  
  # identify reproducible peaks
  # peaks[,group := sample_table$sample_name[file]]
  peaks[,file:=NULL]
  final[,file:=NULL]
  peak_list <- splitDataTable(peaks, "file")
  rm(peaks)
  invisible(gc())
  invisible(sapply(peak_list, countReproduciblePeaks, peak_DT=final))
  
  # keep peaks present in 2 or more individual peak sets
  # keep peaks with score per million >= 5
  final <- final[count >= min_samples & score >= min_score,]
  final[,count := NULL]
  final[,file:=NULL]
  return(final)
}
collapsePeaks <- function(sample_table, chr_sizes,
                          min_samples=2, min_score=5, min_olap=1) {
  # create combined peaks
  peaks <- rbindlist(lapply(sample_table$peak_files, fread), idcol="file")
  if (ncol(peaks) == 7) {
    colnames(peaks) <- c("file", "chr", "start", "end",
                         "name", "score", "strand")
  } else if (ncol(peaks) == 11) {
    colnames(peaks) <- c("file", "chr", "start", "end",
                         "name", "score", "strand",
                         "signalValue", "pValue", "qValue", "peak")
  } else {
    warning(paste0("Peak files did not contain a recognizable number", 
                   " of columns (", ncol(peaks), ")"))
    rm(peaks)
    final <- data.table(chr=character(),
                        start=integer(),
                        end=integer(),
                        name=character(),
                        score=numeric(),
                        strand=character(),
                        signalValue=numeric(),
                        pValue=numeric(),
                        qValue=numeric(),
                        peak=integer())
    return(final)
  }
  setkey(peaks, chr, start, end)
  # keep highest scored peaks
  # split by chromosome to minimize memory requirements
  peaks_by_chr   <- split(peaks, peaks$chr)
  hit_aggregator <- function(x) {
    #message(paste0("x: ", unique(x$chr)))  # DEBUG
    peaksGR <- makeGRangesFromDataFrame(x, keep.extra.columns=FALSE)
    hitsGR  <- suppressWarnings(
      findOverlaps(peaksGR, peaksGR,
                   ignore.strand=TRUE, minoverlap=min_olap))
    hits    <- data.table::data.table(xid=queryHits(hitsGR),
                                      yid=subjectHits(hitsGR))
    setkey(hits, xid)
    scores  <- data.table(index=rep(1:nrow(x)), score=x$score)
    setkey(scores, index)
    out     <- hits[scores, nomatch=0]
    keep    <- out[out[,.I[which.max(score)],by=yid]$V1]
    indices <- unique(keep$xid)
    reduced <- x[indices,]
    reduced[start < 0, start := 0]
    return(reduced)
  }
  final <- rbindlist(lapply(peaks_by_chr, hit_aggregator))
  
  # can't extend past chromosome
  for (i in nrow(chr_sizes)) {
    final[chr == chr_sizes$chr[i] & end > chr_sizes$size[i],
          end := chr_sizes$size[i]]
  }
  
  # identify reproducible peaks
  peaks[,group := sample_table$sample_name[file]]
  peaks[,file:=NULL]
  final[,file:=NULL]
  peak_list <- splitDataTable(peaks, "group")
  rm(peaks)
  invisible(gc())
  invisible(sapply(peak_list, countReproduciblePeaks, peak_DT=final))
  
  # keep peaks present in 2 or more individual peak sets
  # keep peaks with score per million >= 5
  final <- final[count >= min_samples & score >= min_score,]
  final[,count := NULL]
  return(final)
}

consensusPeaks <- function(sample_table, summary_dir, results_subdir, assets,
                           min_samples=2, min_score=5, min_olap=1) {    
  
  # Produce summary output directory (if needed)
  dir.create(summary_dir, showWarnings = FALSE)
  
  setDT(sample_table)[assets[asset == 'chrom_sizes', ],
                      c_path := i.path, on = 'sample_name']
  
  # generate paths to peak files
  sample_table[,peak_files:=.((file.path(
    results_subdir,
    sample_table$sample_name,
    paste0("peak_calling_", sample_table$genome),
    paste0(sample_table$sample_name,
           "_peaks_normalized.narrowPeak"))))]
  
  #Only keep samples with valid peak files
  file_list   <- sample_table$peak_files
  file_exists <- character()
  for (i in 1:length(file_list)) {
    if(file.exists(file.path(file_list[i]))) {
      file_exists <- append(file_exists, file.path(file_list[i]))
    }
  }
  files <- data.table(peak_files=file_exists)
  consensus_peak_files = list()
  if (nrow(files) == 0) {
    return(consensus_peak_files)
  }
  
  sample_table <- unique(
    sample_table[sample_table$peak_files %in% files$peak_files,])
  
  # Need to group by genome, then create a consensus list by genome!
  st_list = splitDataTable(sample_table, "genome")
  
  for (g in names(st_list)) {
    if (nrow(st_list[[g]]) == 1) {
      err_msg = paste0("Found only a single valid peak file for ",
                       g, ".")
      warning(err_msg)
      next
    }
    if (nrow(st_list[[g]]) == 0) {
      warning("Unable to find any valid peak files.")
      warning("Confirm peak files exist for your samples.")
      next
    }
    c_path <- unique(sample_table[genome == g, c_path])
    if (file.exists(c_path)) {
      c_size <- fread(c_path)
      colnames(c_size) <- c("chr", "size")
    } else {
      warning("Unable to load the chromosome sizes file.")
      warning(paste0("Confirm that ", c_path,
                     " is present before continuing."))
      final <- NULL
    }
    message(paste0("Calculating ", g, " consensus peak set from ",
                   nrow(st_list[[g]]), " samples..."))
    final <- collapsePeaks(st_list[[g]], c_size,
                           min_samples, min_score, min_olap)
    
    if (!is.null(final)) {
      # save consensus peak set
      file_name   <- paste0("_", g,"_consensusPeaks.narrowPeak")
      output_file <- file.path(summary_dir,
                               paste0(project_name, file_name))
      fwrite(final, output_file, sep="\t", col.names=FALSE)
      consensus_peak_files <- c(consensus_peak_files, output_file)
      rm(final)
      invisible(gc())
    } else {
      warning("Unable to produce a consensus peak file.")
      warning("Check that individual peak files exist for your samples.")
    }
  }
  
  return(consensus_peak_files)
}
countReproduciblePeaks <- function(peak_list, peak_DT) {
  setkey(peak_DT, chr, start, end)
  setkey(peak_list, chr, start, end)
  hits <- foverlaps(peak_list, peak_DT,
                    by.x=c("chr", "start", "end"),
                    type="any", which=TRUE, nomatch=0)
  # track the number of overlaps of final peak set peaks
  if (!"count" %in% colnames(peak_DT)) {
    peak_DT[hits$yid, count := 1]
    peak_DT[is.na(get("count")), ("count") := 0]
  } else {
    peak_DT[hits$yid, count := get("count") + 1] 
  }
}


# get consensus peaks ------------------------------------------------------

summary_SRA_info = read.table(SRA_table_path, header = T, sep="\t")
if(length(which(colnames(summary_SRA_info) == "file"))!=0){
  summary_SRA_info = summary_SRA_info[, -which(colnames(summary_SRA_info) == "file")] 
} 
summary_SRA_info$sample_name = summary_SRA_info$SRA_ID
summary_SRA_info$protocol = rep("ATAC",nrow(summary_SRA_info))
summary_SRA_info$organism = rep(organism,nrow(summary_SRA_info))
summary_SRA_info$read1 = sapply(summary_SRA_info$sample_name, function(i) paste0(i,"_R1"))
summary_SRA_info$read2 = sapply(summary_SRA_info$sample_name, function(i) paste0(i,"_R2"))
summary_SRA_info$read_type = rep("paired",nrow(summary_SRA_info))

summary_stats = read.table(stats_path, header = T, sep="\t")
head(summary_SRA_info)
head(summary_stats)
print(head(summary_stats$to_keep))
summary_SRA_info = summary_SRA_info[which(summary_SRA_info$sample_name %in% summary_stats$sample[which(summary_stats$to_keep==T)]),]


for(group in unique(summary_SRA_info[,group_id])){
  samples = summary_SRA_info$sample_name[which(summary_SRA_info[,group_id] == group)]
  min_samples = round(length(samples)/2)

  print(paste0("group: ",group))
  
  write.table(summary_SRA_info[which(summary_SRA_info[,group_id] == group),],file = "samples_annotations.csv",row.names = F,col.names = T,quote = F,sep=",")
  
  # Change project name in the pepatac config
  
  system(command = paste0('sed \'2s/.*/name: ',group,'/\' ',pep,' > ',group,'_config.yaml') )
  pep = paste0(group,'_config.yaml')
  
  # Load the project metadata
  prj <- pepr::Project(pep)
  project_name    <- pepr::config(prj)$name
  project_samples <- pepr::sampleTable(prj)$sample_name
  sample_table    <- data.table::data.table(
    sample_name=pepr::sampleTable(prj)$sample_name,
    study=sapply(pepr::sampleTable(prj)$sample_name,function(i) summary_SRA_info$study[which(summary_SRA_info$SRA_ID==i)]),
    genome=pepr::sampleTable(prj)$genome)
    
  
  print(project_samples)
  print(sample_table)
  
  # Specify file locations
  output_dir  <- "./"
  results_dir <- file.path(output_dir,"results/")
  print(results_dir)
  summary_dir <- file.path(output_dir, "summary/")
  print(summary_dir)
  # Produce output directory (if needed)
  dir.create(summary_dir, showWarnings = FALSE)
  
  
  # Identify which chrom_sizes file was used from the project assets
  assets <- createAssetsSummary(project_samples, results_dir)
  print(assets)
  assets$path[which(assets$asset=="chrom_sizes")] = rep(chrom_sizes_path,length(assets$path[which(assets$asset=="chrom_sizes")] ))
  
  # Generate consensus peaks and write to project output directory
  peak_filepath <- consensusPeaks_custom(sample_table, summary_dir, results_dir, assets, min_score)
  
}

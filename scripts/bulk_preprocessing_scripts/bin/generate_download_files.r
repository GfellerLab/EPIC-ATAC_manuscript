

# Parameters --------------------------------------------------------------
project_path = "/mnt/c/Aurelie/postdoc_UNIL/EPIC_ATAC/"
SRA_table_path = as.character(commandArgs(TRUE)[1]) 
res_path=as.character(commandArgs(TRUE)[2]) 
EGA_data=as.logical(commandArgs(TRUE)[3]) 

# Generate download files -------------------------------------------------

summary_SRA_info = read.table(SRA_table_path, header = T, sep="\t")

for(id in unique(summary_SRA_info$SRA_ID)){
  if(grepl("SRR", id)){
    if(!file.exists(paste0(res_path,"/",id,"/PEPATAC_commands.sh"))){
      write.table(summary_SRA_info[which(summary_SRA_info$SRA_ID == id),"SRA_ID",drop=F],
                  file = paste0(id,".txt"),
                  row.names = F,col.names = F,quote = F)
      
    }
  }
  
  if(grepl("EGA", id)){
    if(!file.exists(paste0(res_path,"/",id,"/PEPATAC_commands.sh"))){
      files_id = unlist(strsplit(summary_SRA_info[which(summary_SRA_info$SRA_ID == id),"file"],","))
      write.table(data.frame(files_id = files_id),
                  file = paste0(id,".txt"),
                  row.names = F,col.names = F,quote = F)
      
    }
  }

}


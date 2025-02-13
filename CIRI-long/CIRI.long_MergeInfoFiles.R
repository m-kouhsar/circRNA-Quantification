
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

####################################################################
read.info <- function(info.path,info.id){
  data_ = read.table(file = paste0(info.path,"/",info.id,".info") , header = F,sep = '\t')
  
  info_ = as.data.frame(str_split(data_$V9 , pattern = ";" , simplify = T))
  info_  <- info_[,-10]
  names(info_) = c("circ_id" , "splice_site" , "equivalent_seq"  , "circ_type" , "circ_len" ,"isoform" ,"gene_id" ,"gene_name","gene_type")
  info_$circ_id <- str_remove(info_$circ_id , pattern = "circ_id ")
  info_$splice_site <- str_remove(info_$splice_site , pattern = "splice_site ")
  info_$equivalent_seq <- str_remove(info_$equivalent_seq , pattern = "equivalent_seq ")
  info_$circ_type <- str_remove(info_$circ_type , pattern = "circ_type ")
  info_$circ_len <- str_remove(info_$circ_len , pattern = "circ_len ")
  info_$isoform <- str_remove(info_$isoform , pattern = "isoform ")
  info_$gene_id <- str_remove(info_$gene_id , pattern = "gene_id ")
  info_$gene_name <- str_remove(info_$gene_name , pattern = "gene_name ")
  info_$gene_type <- str_remove(info_$gene_type , pattern = "gene_type ")
  data_ = data_[,c(-2,-3,-8,-9)]
  names(data_) <- c("chrom","start","end",info.id,"strand")
  data_ <- cbind.data.frame(data_ , info_)
  data_ <- data_[,c(6,1:3,5,7:14,4)]
  return(data_)
}
####################################################################
arguments <- commandArgs(T)

collapse.dir <- trimws(arguments[1])
OutPrefix <- trimws(arguments[2])

if(is.na(OutPrefix)){
  OutPrefix <- ""
}

collapse.path = list.dirs(path = collapse.dir , full.names = T , recursive = F)
collapse.id = list.dirs(path = collapse.dir , full.names = F , recursive = F)

results <- read.info(info.path = collapse.path[1],info.id = collapse.id[1])

for (i in 2:length(collapse.path)) {
  
  new_data <- read.info(info.path = collapse.path[i],info.id = collapse.id[i])
  results <- full_join(results, new_data, by = names(new_data)[1:13] )
  
}

results <- mutate_all(results,~replace_na(., 0))

write.table(results , file = paste(OutPrefix , "CIRIlong.circInfo.txt" ,sep = ifelse(OutPrefix == "" , "",".")),
            row.names = F,col.names = T,sep = '\t',quote = F)



library(stringr)

arguments <- commandArgs(T)

fastq.dir <- trimws(arguments[1])
call.dir <- trimws(arguments[2])
collapse.dir <- trimws(arguments[3])
OutPrefix <- trimws(arguments[4])

if(is.na(OutPrefix)){
  OutPrefix <- ""
  }

fastq.file <- list.files(path = fastq.dir , pattern = ".fastq" , full.names = F , recursive = F)
fastq.id <- str_remove(fastq.file , pattern = ".fastq")

call.id = list.dirs(path = call.dir , full.names = F , recursive = F)
call.path = list.dirs(path = call.dir , full.names = T , recursive = F)

collapse.id = list.dirs(path = collapse.dir , full.names = F , recursive = F)
collapse.path = list.dirs(path = collapse.dir , full.names = T , recursive = F)

results <- data.frame(Sample = fastq.id , fastq = fastq.file , Called = (fastq.id %in% call.id),
                      Collapsed = (fastq.id %in% collapse.id), 
                      Call.cand_circ.fa = F , Call.low_confidence.fa = F,
                      Call.json = F , Call.log = F , Collapse.expression = F , 
                      Collapse.isoforms = F , Collapse.info = F , Collapse.reads = F , 
                      Collapse.log = F)

for (i in 1:length(call.path)) {
  
  f = list.files(path = call.path[i] , full.names = F , recursive = F)
  
  index = which(results$Sample == call.id[i])
  
  results$Call.cand_circ.fa[index] = any(str_detect(f , pattern = ".cand_circ.fa"))
  results$Call.low_confidence.fa[index] = any(str_detect(f , pattern = ".low_confidence.fa"))
  results$Call.json[index] = any(str_detect(f , pattern = ".json"))
  results$Call.log[index] = any(str_detect(f , pattern = ".log"))
}

for (i in 1:length(collapse.path)) {
  
  f = list.files(path = collapse.path[i] , full.names = F , recursive = F)
  
  index = which(results$Sample == collapse.id[i])
  
  results$Collapse.expression[index] = any(str_detect(f , pattern = ".expression"))
  results$Collapse.isoforms[index] = any(str_detect(f , pattern = ".isoform"))
  results$Collapse.info[index] = any(str_detect(f , pattern = ".info"))
  results$Collapse.reads[index] = any(str_detect(f , pattern = ".reads"))
  results$Collapse.log[index] = any(str_detect(f , pattern = ".log"))
}

index = apply(results[,3:13] , 1 , all )
results = results[!index , ]
write.csv(results , file = paste(OutPrefix , "CIRIlong.ResultsCheck.csv" , 
                                 sep = ifelse(OutPrefix == "" , "",".")),row.names = F )

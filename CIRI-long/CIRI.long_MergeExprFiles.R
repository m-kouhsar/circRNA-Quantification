
suppressMessages(library(stringr))
suppressMessages(library(dplyr))

arguments <- commandArgs(T)

collapse.dir <- trimws(arguments[1])
OutPrefix <- trimws(arguments[2])

if(is.na(OutPrefix)){
  OutPrefix <- ""
}

collapse.path = list.dirs(path = collapse.dir , full.names = T , recursive = F)
collapse.id = list.dirs(path = collapse.dir , full.names = F , recursive = F)

results = read.table(file = paste0(collapse.path[1],"/",collapse.id[1],".expression") , header = T,sep = '\t')

for (i in 2:length(collapse.path)) {
  data_ = read.table(file = paste0(collapse.path[i],"/",collapse.id[i],".expression") , header = T,sep = '\t')
  results <- full_join(results , data_ , by= "circ_ID",) 
}

write.table(results , file = paste(OutPrefix , "CIRIlong.circExpr.txt" ,sep = ifelse(OutPrefix == "" , "",".")),
            row.names = F,col.names = T,sep = '\t',quote = F)

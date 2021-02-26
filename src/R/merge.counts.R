library(data.table)
library(tidyr)
library(dplyr)
library(glue)

#------------------------------------------------------------------------------
#   Read directory with files previously created to grep read counts
#------------------------------------------------------------------------------
files.list <- list.files("/home/jp/eclipse-workspace/Toulouse/data/results/",pattern="(*)-differential-gene-selection.txt$")
final.file <- "/home/jp/eclipse-workspace/Toulouse/data/results/JP_final.read.tsv"

#Check its existence
if (file.exists(final.file)) {
  #Delete file if it exists
  file.remove(final.file)
}
for (file in files.list){
  
  print(file)
  asbolutepath2file <- glue("/home/jp/eclipse-workspace/Toulouse/data/results/{file}")
  print(asbolutepath2file)
  
  if (!exists("dataset") ){
    dataset <- fread(asbolutepath2file,data.table=F)
    dataset <- subset(dataset,select=-c(logFC,logCPM,LR,PValue,FDR))
    print(head(dataset))
    
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <- fread(asbolutepath2file,data.table=F)
    temp_dataset <- subset(temp_dataset,select=-c(logFC,logCPM,LR,PValue,FDR))
    dataset <- inner_join(dataset,temp_dataset, by =c("genes"),keep=FALSE)
    rm(temp_dataset)
  }
  
}

colnames(dataset) <- gsub(".x",'', colnames(dataset), fixed=TRUE)
head (dataset)
write.table(dataset,file=final.file,quote=F,row.names=F,sep="\t")

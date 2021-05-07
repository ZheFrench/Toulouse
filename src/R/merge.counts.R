#################################################################
#
# date: April 02, 2021
# platform: Ubuntu 10.04
# R.version : 4.0.3
# author: Villemin Jean-Philippe
# team: Bioinformatique et biologie des systèmes du cancer : J. Colinge 
# Institute : IRCM
#
# merge.counts.R
# Usage : 
# 
# merge.counts.R
# 
# Description : 
#
# 
#
#################################################################

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(glue))

# ==========   Settings =============================

base.dir <- "/data/villemin/code/Toulouse-rnaseq/"

cellline = "h4"

#------------------------------------------------------------------------------
#   Read directory with files previously created to grep read counts
#------------------------------------------------------------------------------
files.list <- list.files(glue("{base.dir}/data/results/{cellline}/"),pattern="(*)-differential-gene-selection.txt$")

#output
final.file <- glue("{base.dir}/data/results/{cellline}_final.read.tsv")


#Check its existence
if (file.exists(final.file)) {
  #Delete file if it exists
  file.remove(final.file)
}
for (file in files.list){
  
  print(file)
  asbolutepath2file <- glue("{base.dir}/data/results/{cellline}/{file}")
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
#colnames(dataset) <- gsub("*.right",'', colnames(dataset), fixed=TRUE)
#colnames(dataset) <- gsub(".x",'', colnames(dataset), fixed=TRUE)

clean.dataset <- dataset[grep(".y$", names(dataset), invert = TRUE)]
colnames(clean.dataset) <- gsub(".x",'', colnames(clean.dataset), fixed=TRUE)
head(clean.dataset)

write.table(clean.dataset,file=final.file,quote=F,row.names=F,sep="\t")


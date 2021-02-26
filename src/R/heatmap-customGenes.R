library(data.table)
library(gplots)
library(ComplexHeatmap)
library(gtools)
library(dplyr)
library(tidyr)
library(glue)

pc9.dormance <- fread("../../data/results/JP-PC9_Erlo_DTC-PC9_CT-differential-gene-selection.txt",data.table=F)
pc9.expansion   <- fread("../../data/results/JP-PC9_3_DTEC_75d-PC9_Erlo_DTC-differential-gene-selection.txt",data.table=F)

h4.dormance <- fread("../../data/results/JP-H4006_Erlo_DTC-H4006_CT_24h-differential-gene-selection.txt",data.table=F)
h4.expansion   <- fread("../../data/results/JP-H4006_Erlo_DTEC-H4006_Erlo_DTC-differential-gene-selection.txt",data.table=F)

h3.control <- fread("../../data/results/JP-H3255_Erlo_21d-H3255_NT-differential-gene-selection.txt",data.table=F)

my.path <- "../../data/custom.geneslist.txt"
#my.path <- "../../data/rna.txt"
#my.path <- "../../data/dna.txt"
#my.path <- "../../data/chrom.txt"

my.genes <- fread(my.path,data.table=F)
filename.genes.list <- gsub(".txt","",basename(my.path))

print(filename.genes.list)

# list to character
genes <- unlist(my.genes, use.names=FALSE)

data.pc9 <- inner_join(select(pc9.dormance, -c(3:last_col())), select(pc9.expansion, -c(3:last_col())), by = c("genes"),keep=FALSE)
colnames(data.pc9) <- gsub('logFC.x', 'pc9.dormance', colnames(data.pc9), fixed=TRUE)
colnames(data.pc9) <- gsub('logFC.y', 'pc9.expansion', colnames(data.pc9), fixed=TRUE)

data.h4 <- inner_join(select(h4.dormance, -c(3:last_col())), select(h4.expansion, -c(3:last_col())), by =c("genes"),keep=FALSE) 
colnames(data.h4) <- gsub('logFC.x', 'h4.dormance', colnames(data.h4), fixed=TRUE)
colnames(data.h4) <- gsub('logFC.y', 'h4.expansion', colnames(data.h4), fixed=TRUE)

data.h4.pc9 <- inner_join(data.pc9,data.h4, by =c("genes"),keep=FALSE)

data.h4.pc9.h3 <- inner_join(select(h3.control, -c(3:last_col())) , data.h4.pc9 , by =c("genes"),keep=FALSE) 
colnames(data.h4.pc9.h3) <- gsub('logFC', 'h3.control', colnames(data.h4.pc9.h3), fixed=TRUE)

write.table(data.h4.pc9.h3[data.h4.pc9.h3$genes %in% genes,],file=glue("../../data/results/matrice_heatmap_{filename.genes.list}_patterns.txt"),quote=F,row.names=F,sep="\t")

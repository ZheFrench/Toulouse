library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(optparse)
library(tidyr)
library(dplyr)
library(viridis)
library(GSVA)
library(fgsea)
library(glue)

#https://github.com/rcastelo/GSVA
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6500180/pdf/pnas.201818210.pdf
#https://www.rdocumentation.org/packages/GSVA/versions/1.20.0/topics/gsva
#https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html#gsea-analysis
#https://www.programmersought.com/article/55084356499/
#ssgsea.norm=TRUE,

#------------------------------------------------------------------------------
#  Functions 
#------------------------------------------------------------------------------

matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[,1]
  m
}

#------------------------------------------------------------------------------
#   ssGsea
#------------------------------------------------------------------------------

final.file <- "/home/jp/eclipse-workspace/Toulouse/data/results/JP_final.read.tsv"
dataframe.expression <- fread(final.file,data.table=F)

C2<- read.gmt("/home/jp/eclipse-workspace/database/MSigDB/c2.all.v7.2.symbols.gmt")
C5<- read.gmt("/home/jp/eclipse-workspace/database/MSigDB/c5.go.v7.2.symbols.gmt")

# Use Fsgea Function
h.All <- gmtPathways("/home/jp/eclipse-workspace/database/MSigDB/h.all.v7.2.symbols.gmt")

# Show the first few pathways, and within those, show only the first few genes. 
h.All %>% head() %>% lapply(head)

matrix.expression <-matrix.please(dataframe.expression)

ssGseaScore <- gsva(matrix.expression, h.All , method=c( "ssgsea"), kcdf=c( "Poisson"),ssgsea.norm=FALSE,
     abs.ranking=FALSE, min.sz=1, max.sz=Inf,parallel.sz=4,verbose=TRUE)
     #no.bootstraps=0, bootstrap.percent = .632, parallel.sz=0, 
     #parallel.type="SOCK", mx.diff=TRUE, 
     #tau=switch(method, gsva=1, ssgsea=0.25, NA),
     #kernel=TRUE, ssgsea.norm=TRUE, verbose=TRUE
head(ssGseaScore)
write.table(ssGseaScore,file=final.file,quote=F,row.names=F,sep="\t")


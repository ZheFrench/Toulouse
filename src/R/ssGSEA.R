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
library(rstatix)

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6500180/pdf/pnas.201818210.pdf
#https://www.rdocumentation.org/packages/GSVA/versions/1.20.0/topics/gsva
#https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html#gsea-analysis
#https://www.programmersought.com/article/55084356499/


#------------------------------------------------------------------------------
#  Functions 
#------------------------------------------------------------------------------

matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[,1]
  m
}

filter.byWilcox<-function(mydataframe) {
  
  dt.ssgeascore.wtnames <- subset (mydataframe, select=-c(names))
  
  dt.ssgeascore.transposed <- as.data.frame(t(dt.ssgeascore.wtnames))
  
  #H3255_Erlo_21d_1
  #H3255_Erlo_21d_2
  #H3255_Erlo_21d_3	
  #H3255_NT1	
  #H3255_NT2
  #H3255_NT3
  #H4006_Erlo_DTC_1
  #H4006_Erlo_DTC_2
  #H4006_Erlo_DTC_3
  #H4006_CT_24h_1
  #H4006_CT_24h_2
  #H4006_CT_24h_3
  #H4006_Erlo_DTEC_1
  #H4006_Erlo_DTEC_2
  #H4006_Erlo_DTEC_3
  #H3255_Erlo_21d_1
  #H3255_Erlo_21d_2
  #H3255_Erlo_21d_3
  
  dt.ssgeascore.transposed.final <- dt.ssgeascore.transposed %>% 
    mutate(group = if_else(rownames(dt.ssgeascore.transposed) %in% c("H4006_Erlo_DTC_1","H4006_Erlo_DTC_2","H4006_Erlo_DTC_3","H4006_Erlo_DTEC_1","H4006_Erlo_DTEC_2","H4006_Erlo_DTEC_3"), '0', '1'))
  
  dt.ssgeascore.transposed.final <- mutate_all(dt.ssgeascore.transposed.final, function(x) as.numeric(as.character(x)))

  vector <- c()
  for (i in (1:dim(dt.ssgeascore.transposed.final)[2]) ) {
      
      if(colnames(dt.ssgeascore.transposed.final[i]) =="group"){ break }
      
      str <- paste(c(colnames(dt.ssgeascore.transposed.final[i]),"group"),collapse ="~")
      #print(str)
      rez <- wilcox_test(
        dt.ssgeascore.transposed.final,
        formula (str),
        comparisons = NULL,
        ref.group = NULL,
        p.adjust.method = "holm",
        paired = FALSE,
        exact = TRUE,
        alternative = "two.sided",
        mu = 0,
        conf.level = 0.95,
        detailed = FALSE
      )
      
      if (rez$p < 0.001) {
        vector <- c(vector, colnames(dt.ssgeascore.transposed.final[i]))
      }
      
    }

  subset <- mydataframe[mydataframe$names %in% vector,]

return(subset)
}

filter.byterm<-function(mydataframe,term) {
  
  ind <- grep(term,mydataframe,perl=T)
 
 return(ssGseaScore.clean[ind,]) 
}
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

#   ssGsea


cellline <- "h4"
type     <-"c3"

final.file <- glue("/home/jp/eclipse-workspace/Toulouse/data/results/{cellline}_final.read.tsv")
dataframe.expression <- fread(final.file,data.table=F)

#gmt <- gmtPathways("/home/jp/eclipse-workspace/database/MSigDB/c2.all.v7.2.symbols.gmt")# 50
gmt <- gmtPathways("/home/jp/eclipse-workspace/database/MSigDB/c3.tft.gtrd.v7.2.symbols.gmt")#348
#gmt <- gmtPathways("/home/jp/eclipse-workspace/database/MSigDB/c6.all.v7.2.symbols.gmt") # 189


#gmt <- gmtPathways("/home/jp/eclipse-workspace/database/MSigDB/c5.go.v7.2.symbols.gmt")
#gmt <- gmtPathways("/home/jp/eclipse-workspace/database/MSigDB/h.all.v7.2.symbols.gmt")

# Show the first few pathways, and within those, show only the first few genes. 
#h.All %>% head() %>% lapply(head)

matrix.expression <-matrix.please(dataframe.expression)

ssGseaScore <- gsva(matrix.expression, gmt , method=c( "ssgsea"), kcdf=c( "Poisson"),ssgsea.norm=TRUE,
     abs.ranking=FALSE, min.sz=5, max.sz=250,parallel.sz=4,verbose=TRUE)

names <- rownames(ssGseaScore)
ssGseaScore.clean <- cbind(names,ssGseaScore)

dt.ssgeascore <- as.data.frame (ssGseaScore.clean)

#scores <- filter.byWilcox(dt.ssgeascore)

#term <- "RESISTANCE"
#scores <- filter.byterm(ssGseaScore.clean,term)

# test to clean up, and reduce number...
#write.table(dt.ssgeascore,file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/{type}_{cellline}_ssGSEA.tsv"),quote=F,row.names=F,sep="\t")
write.table(dt.ssgeascore,file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/{type}_{cellline}_ssGSEA.tsv"),quote=F,row.names=F,sep="\t")




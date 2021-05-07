suppressPackageStartupMessages(library(data.table))
library(clusterProfiler)
library(org.Hs.eg.db)
library(optparse)
library(tidyr)
library(dplyr)
library(viridis)
library (glue)
library(ggplot2)

#Over-representation test(Boyle et al. 2004)
# hypergeometric test
#https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/
#http://yulab-smu.top/clusterProfiler-book/chapter3.html
#https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html
#http://yulab-smu.top/clusterProfiler-book/chapter11.html

option_list = list(
  make_option(c("-s", "--subset"), type="character", default=NULL, help="Subset of genes", metavar="character")
); 

parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = 0);
opt <- arguments$options

GO_ALL <-  read.gmt("/home/jp/eclipse-workspace/database/MSigDB/c5.go.v7.2.symbols.gmt")  

#/home/jp/eclipse-workspace/Toulouse/data/results/GO/pc9.genes.1A.2A.txt
#http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/GO/
# ==========   Settings =============================
GO_bp <- read.gmt("/home/jp/eclipse-workspace/database/GO/Human_GO_bp_no_GO_iea_symbol.gmt")

GO_bp <- read.gmt("/home/jp/eclipse-workspace/database/GO/Human_GO_bp_no_GO_iea_symbol.gmt")
GO_mf <- read.gmt("/home/jp/eclipse-workspace/database/GO/Human_GO_mf_no_GO_iea_symbol.gmt")
Reactome <- read.gmt("/home/jp/eclipse-workspace/database/Reactome/ReactomePathways.25022021.gmt")
#file.gmt <- "/home/jp/eclipse-workspace/database/MSigDB/c6.all.v7.2.symbols.gmt"

Reactome <- Reactome[, c("term", "gene")]

GO_bp <- GO_bp %>% separate(term, c("Name","Goterm"), sep = "%GOBP%")
GO_bpterm2gene <- GO_bp[, c("Goterm", "gene")]
GO_bpterm2name <- GO_bp[, c("Goterm", "Name")]

GO_mf <- GO_mf %>% separate(term, c("Name","Goterm"), sep = "%GOMF%")

GO_mfterm2gene <- GO_mf[, c("Goterm", "gene")]
GO_mfterm2name <- GO_mf[, c("Goterm", "Name")]

head(GO_ALL)
GO_ALLterm2gene <- GO_ALL[, c("term", "gene")]


dir <- "genes_2_2"


# ==========================================================================================================
# ================ H4
# ==========================================================================================================

print ("H4")
subset <- opt$subset
genes   <- fread(glue("/home/jp/eclipse-workspace/Toulouse/data/results/{dir}/h4.genes.pattern{subset}.txt"),data.table=F)

# list to character
genes <- unlist(genes, use.names=FALSE)

h4.genes.1A.2A.dataframe <- as.data.frame(genes)
h4.genes.1A.2A.dataframe$geneslist  <- "genes.{subset}"
h4.genes.1A.2A.dataframe$cellline  <- "h4"

genes.enriched.GO_bp <- enricher(genes, TERM2GENE=GO_bpterm2gene, TERM2NAME=GO_bpterm2name,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)
genes.enriched.GO_mf <- enricher(genes, TERM2GENE=GO_mfterm2gene, TERM2NAME=GO_mfterm2name,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)
genes.enriched.Reactome <- enricher(genes, TERM2GENE=Reactome, TERM2NAME=NA,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)
genes.enriched.GO_ALL <- enricher(genes, TERM2GENE=GO_ALLterm2gene, TERM2NAME=NA,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)

if (length(genes.enriched.Reactome) >= 1) {
  
png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/Reactome.h4.genes.{subset}.png"),width=900,height=1+dim(genes.enriched.Reactome)[[1]]*50)
barplot(genes.enriched.Reactome  ,  drop = TRUE,  showCategory = 50,  title = "Reactome", font.size = 12) +scale_fill_viridis(option="viridis",direction = -1)
dev.off()

png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/Reactome.h4.genes.{subset}-heatplot.png"),width=900,height=1+dim(genes.enriched.Reactome)[[1]]*50)
heatplot(genes.enriched.Reactome)+ theme(  axis.text.y = element_text( size = 12))
dev.off()

}

stop()

if (length(genes.enriched.GO_bp) >= 1) {
  
write.csv(genes.enriched.GO_bp   ,file = glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/BP.h4.genes.{subset}.csv"), row.names=FALSE)
png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/BP.h4.genes.{subset}.png"),width=900,height=1+dim(genes.enriched.GO_bp)[[1]]*40)
barplot(genes.enriched.GO_bp    ,  drop = TRUE,  showCategory = 50,  title = "Go-BP", font.size = 12)+scale_fill_viridis(option="viridis",direction = -1)
dev.off()

png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/BP.h4.genes.{subset}-heatplot.png"),width=700,height=1+dim(genes.enriched.GO_bp)[[1]]*30)
heatplot(genes.enriched.GO_bp)+ theme(  axis.text.y = element_text( size = 12))
dev.off()
}

if (length(genes.enriched.GO_mf) >= 1) {
  
write.csv(genes.enriched.GO_mf  ,file = glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/MF.h4.genes.{subset}.csv"), row.names=FALSE)
png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/MF.h4.genes.{subset}.png"),width=700,height=1+dim(genes.enriched.GO_mf)[[1]]*50)
barplot(genes.enriched.GO_mf  ,  drop = TRUE,  showCategory = 50,  title = "Go-MF", font.size = 12)+scale_fill_viridis(option="viridis",direction = -1)
dev.off()

png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/MF.h4.genes.{subset}-heatplot.png"),width=700,height=1+dim(genes.enriched.GO_mf)[[1]]*30)
heatplot(genes.enriched.GO_mf)+ theme(  axis.text.y = element_text( size = 12))
dev.off()
}


if (length(genes.enriched.GO_ALL) >= 1) {
  
  write.csv(genes.enriched.GO_ALL  ,file = glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/ALL.h4.genes.{subset}.csv"), row.names=FALSE)
  png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/ALL.h4.genes.{subset}.png"),width=700,height=1+dim(genes.enriched.GO_mf)[[1]]*50)
  barplot(genes.enriched.GO_ALL  ,  drop = TRUE,  showCategory = 50,  title = "Go-MF", font.size = 12)+scale_fill_viridis(option="viridis",direction = -1)
  dev.off()
  
  png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/ALL.h4.genes.{subset}-heatplot.png"),width=700,height=1+dim(genes.enriched.GO_mf)[[1]]*30)
  heatplot(genes.enriched.GO_ALL)+ theme(  axis.text.y = element_text( size = 12))
  dev.off()
}

# ==========================================================================================================
# ================ PC9
# ==========================================================================================================

print("PC9")
genes   <- fread(glue("/home/jp/eclipse-workspace/Toulouse/data/results/{dir}/pc9.genes.pattern{subset}.txt"),data.table=F)

# list to character
genes <- unlist(genes, use.names=FALSE)

pc9.genes.1A.2A.dataframe <- as.data.frame(genes)
pc9.genes.1A.2A.dataframe$geneslist  <- glue("genes.{subset}")
pc9.genes.1A.2A.dataframe$cellline  <- "pc9"


genes.enriched.GO_bp <- enricher(genes, TERM2GENE=GO_bpterm2gene, TERM2NAME=GO_bpterm2name,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)
genes.enriched.GO_mf <- enricher(genes, TERM2GENE=GO_mfterm2gene, TERM2NAME=GO_mfterm2name,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)
genes.enriched.Reactome <- enricher(genes, TERM2GENE=Reactome, TERM2NAME=NA,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)
genes.enriched.GO_ALL <- enricher(genes, TERM2GENE=GO_ALLterm2gene, TERM2NAME=NA,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)

if (length(genes.enriched.Reactome) >= 1) {
  
write.csv(genes.enriched.Reactome  ,file = glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/Reactome.pc9.genes.{subset}.csv"), row.names=FALSE)
png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/Reactome.pc9.genes.{subset}.png"),width=900,height=1+dim(genes.enriched.Reactome)[[1]]*30)
barplot(genes.enriched.Reactome  ,  drop = TRUE,  showCategory = 50,  title = "Reactome", font.size = 16) +scale_fill_viridis(option="viridis",direction = -1)
dev.off()

png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/Reactome.pc9.genes.{subset}-heatplot.png"),width=900,height=1+dim(genes.enriched.Reactome)[[1]]*50)
heatplot(genes.enriched.Reactome)+ theme(  axis.text.y = element_text( size = 12))
dev.off()

}

if (length(genes.enriched.GO_bp) >= 1) {
  
write.csv(genes.enriched.GO_bp   ,file = glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/BP.pc9.genes.{subset}.csv"), row.names=FALSE)
png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/BP.pc9.genes.{subset}.png"),width=900,height=1+dim(genes.enriched.GO_bp)[[1]]*7)
barplot(genes.enriched.GO_bp    ,  drop = TRUE,  showCategory = 50,  title = "Go-BP", font.size = 14)+scale_fill_viridis(option="viridis",direction = -1)
dev.off()

png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/BP.pc9.genes.{subset}.csv-heatplot.png"),width=900,height=1+dim(genes.enriched.GO_bp)[[1]]*50)
heatplot(genes.enriched.GO_bp) + theme(  axis.text.y = element_text( size = 12))
dev.off()

}

if (length(genes.enriched.GO_mf) >= 1) {
write.csv(genes.enriched.GO_mf  ,file = glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/MF.pc9.genes.{subset}.csv"), row.names=FALSE)
png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/MF.pc9.genes.{subset}png",width=900,height=1+dim(genes.enriched.GO_mf)[[1]]*50)
barplot(genes.enriched.GO_mf  ,  drop = TRUE,  showCategory = 50,  title = "Go-MF", font.size = 12)+scale_fill_viridis(option="viridis",direction = -1)
dev.off()

png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/MF.pc9.genes.{subset}.csv-heatplot.png"),width=900,height=1+dim(genes.enriched.GO_mf)[[1]]*50)
heatplot(genes.enriched.GO_mf)+ theme(  axis.text.y = element_text( size = 12))
dev.off()
}


if (length(genes.enriched.GO_ALL) >= 1) {
  
  write.csv(genes.enriched.GO_ALL  ,file = glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/ALL.pc9.genes.{subset}.csv"), row.names=FALSE)
  png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/ALL.pc9.genes.{subset}.png"),width=700,height=1+dim(genes.enriched.GO_mf)[[1]]*50)
  barplot(genes.enriched.GO_ALL  ,  drop = TRUE,  showCategory = 50,  title = "Go-MF", font.size = 12)+scale_fill_viridis(option="viridis",direction = -1)
  dev.off()
  
  png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/ALL.pc9.genes.{subset}-heatplot.png"),width=700,height=1+dim(genes.enriched.GO_mf)[[1]]*30)
  heatplot(genes.enriched.GO_ALL)+ theme(  axis.text.y = element_text( size = 12))
  dev.off()
}


# ==========================================================================================================
# ================ H3225
# ==========================================================================================================

print("H3")
for(subset in c("UP","DOWN")) {
    
  genes   <- fread(glue("/home/jp/eclipse-workspace/Toulouse/data/results/{dir}/h3.genes.pattern{subset}.txt"),data.table=F)
  
  # list to character
  genes <- unlist(genes, use.names=FALSE)
  
  pc9.genes.1A.2A.dataframe <- as.data.frame(genes)
  pc9.genes.1A.2A.dataframe$geneslist  <- glue("genes.{subset}")
  pc9.genes.1A.2A.dataframe$cellline  <- "h3"
  
  
  genes.enriched.GO_bp <- enricher(genes, TERM2GENE=GO_bpterm2gene, TERM2NAME=GO_bpterm2name,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)
  genes.enriched.GO_mf <- enricher(genes, TERM2GENE=GO_mfterm2gene, TERM2NAME=GO_mfterm2name,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)
  genes.enriched.Reactome <- enricher(genes, TERM2GENE=Reactome, TERM2NAME=NA,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)
  genes.enriched.GO_ALL <- enricher(genes, TERM2GENE=GO_ALLterm2gene, TERM2NAME=NA,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)
  
  if (length(genes.enriched.Reactome) >= 1) {
    
  write.csv(genes.enriched.Reactome  ,file = glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/Reactome.h3.genes.{subset}.csv"), row.names=FALSE)
  png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/Reactome.h3.genes.{subset}.png"),width=900,height=1+dim(genes.enriched.Reactome)[[1]]*30)
  print(barplot(genes.enriched.Reactome  ,  drop = TRUE,  showCategory = 50,  title = "Reactome", font.size = 16) +scale_fill_viridis(option="viridis",direction = -1))
  dev.off()
  
  png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/Reactome.h3.genes.{subset}-heatplot.png"),width=900,height=1+dim(genes.enriched.Reactome)[[1]]*50)
  print(heatplot(genes.enriched.Reactome)+ theme(  axis.text.y = element_text( size = 12)))
  dev.off()
  
  }
  
  if (length(genes.enriched.GO_bp) >= 1) {
    
  write.csv(genes.enriched.GO_bp   ,file = glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/BP.h3.genes.{subset}.csv"), row.names=FALSE)
  png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/BP.h3.genes.{subset}.png"),width=900,height=1+dim(genes.enriched.GO_bp)[[1]]*7)
  print(barplot(genes.enriched.GO_bp    ,  drop = TRUE,  showCategory = 50,  title = "Go-BP", font.size = 14)+scale_fill_viridis(option="viridis",direction = -1))
  dev.off()
  
  png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/BP.h3.genes.{subset}.csv-heatplot.png"),width=900,height=1+dim(genes.enriched.GO_bp)[[1]]*50)
  print(heatplot(genes.enriched.GO_bp) + theme(  axis.text.y = element_text( size = 12)))
  dev.off()
  
  }
  
  if (length(genes.enriched.GO_mf) >= 1) {
  write.csv(genes.enriched.GO_mf  ,file = glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/MF.h3.genes.{subset}.csv"), row.names=FALSE)
  png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/MF.h3.genes.{subset}.png"),width=900,height=1+dim(genes.enriched.GO_mf)[[1]]*50)
  print(barplot(genes.enriched.GO_mf  ,  drop = TRUE,  showCategory = 50,  title = "Go-MF", font.size = 12)+scale_fill_viridis(option="viridis",direction = -1))
  dev.off()
  
  png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/MF.h3.genes.{subset}.csv-heatplot.png"),width=900,height=1+dim(genes.enriched.GO_mf)[[1]]*50)
  print(heatplot(genes.enriched.GO_mf)+ theme(  axis.text.y = element_text( size = 12)))
  dev.off()
  
  }
  
  
  if (length(genes.enriched.GO_ALL) >= 1) {
    
    write.csv(genes.enriched.GO_ALL  ,file = glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/ALL.h3.genes.{subset}.csv"), row.names=FALSE)
    png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/ALL.h3.genes.{subset}.png"),width=700,height=1+dim(genes.enriched.GO_mf)[[1]]*50)
    barplot(genes.enriched.GO_ALL  ,  drop = TRUE,  showCategory = 50,  title = "Go-MF", font.size = 12)+scale_fill_viridis(option="viridis",direction = -1)
    dev.off()
    
    png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/GO/ALL.h3.genes.{subset}-heatplot.png"),width=700,height=1+dim(genes.enriched.GO_mf)[[1]]*30)
    heatplot(genes.enriched.GO_ALL)+ theme(  axis.text.y = element_text( size = 12))
    dev.off()
  }
  
}

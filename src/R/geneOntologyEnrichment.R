library(data.table)
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

#option_list = list(
#  make_option(c("-f", "--file"), type="character", default=NULL, help="Absolute File Input Path", metavar="character"),
#); 

#parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
#arguments = parse_args(parser, positional_arguments = 0);
#opt <- arguments$options
#opt$file
#/home/jp/eclipse-workspace/Toulouse/data/results/GO/pc9.genes.1A.2A.txt

# ==========   Settings =============================

GO_bp <- read.gmt("/home/jp/eclipse-workspace/database/GO/Human_GO_bp_no_GO_iea_symbol.gmt")
GO_mf <- read.gmt("/home/jp/eclipse-workspace/database/GO/Human_GO_mf_no_GO_iea_symbol.gmt")
Reactome <- read.gmt("/home/jp/eclipse-workspace/database/Reactome/ReactomePathways.25022021.gmt")

Reactome <- Reactome[, c("term", "gene")]

GO_bp <- GO_bp %>% separate(term, c("Name","Goterm"), sep = "%GOBP%")
GO_bpterm2gene <- GO_bp[, c("Goterm", "gene")]
GO_bpterm2name <- GO_bp[, c("Goterm", "Name")]

GO_mf <- GO_mf %>% separate(term, c("Name","Goterm"), sep = "%GOMF%")
GO_mfterm2gene <- GO_mf[, c("Goterm", "gene")]
GO_mfterm2name <- GO_mf[, c("Goterm", "Name")]


dir <- "genes_3_2"


# ==========================================================================================================
# ================ H4
# ==========================================================================================================



genes   <- fread(glue("/home/jp/eclipse-workspace/Toulouse/data/results/{dir}/h4.genes.1A.2A.txt"),data.table=F)


# list to character
genes <- unlist(genes, use.names=FALSE)

h4.genes.1A.2A.dataframe <- as.data.frame(genes)
h4.genes.1A.2A.dataframe$geneslist  <- "genes.1A.2A"
h4.genes.1A.2A.dataframe$cellline  <- "h4"

genes.enriched.GO_bp <- enricher(genes, TERM2GENE=GO_bpterm2gene, TERM2NAME=GO_bpterm2name,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)
genes.enriched.GO_mf <- enricher(genes, TERM2GENE=GO_mfterm2gene, TERM2NAME=GO_mfterm2name,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)
genes.enriched.Reactome <- enricher(genes, TERM2GENE=Reactome, TERM2NAME=NA,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)

png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/Reactome.h4.genes.1A.2A.png",width=900,height=1+dim(genes.enriched.Reactome)[[1]]*50)
barplot(genes.enriched.Reactome  ,  drop = TRUE,  showCategory = 50,  title = "Reactome", font.size = 12) +scale_fill_viridis(option="viridis",direction = -1)
dev.off()

png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/Reactome.h4.genes.1A.2A-heatplot.png",width=900,height=1+dim(genes.enriched.Reactome)[[1]]*50)
heatplot(genes.enriched.Reactome)+ theme(  axis.text.y = element_text( size = 12))
dev.off()

write.csv(genes.enriched.GO_bp   ,file = "/home/jp/eclipse-workspace/Toulouse/data/results/GO/BP.h4.genes.1A.2A.csv", row.names=FALSE)
png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/BP.h4.genes.1A.2A.png",width=900,height=1+dim(genes.enriched.GO_bp)[[1]]*40)
barplot(genes.enriched.GO_bp    ,  drop = TRUE,  showCategory = 50,  title = "Go-BP", font.size = 12)+scale_fill_viridis(option="viridis",direction = -1)
dev.off()

png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/BP.h4.genes.1A.2A-heatplot.png",width=700,height=1+dim(genes.enriched.GO_bp)[[1]]*30)
heatplot(genes.enriched.GO_bp)+ theme(  axis.text.y = element_text( size = 12))
dev.off()

write.csv(genes.enriched.GO_mf  ,file = "/home/jp/eclipse-workspace/Toulouse/data/results/GO/MF.h4.genes.1A.2A.csv", row.names=FALSE)
png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/MF.h4.genes.1A.2A.png",width=700,height=1+dim(genes.enriched.GO_mf)[[1]]*50)
barplot(genes.enriched.GO_mf  ,  drop = TRUE,  showCategory = 50,  title = "Go-MF", font.size = 12)+scale_fill_viridis(option="viridis",direction = -1)
dev.off()

png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/MF.h4.genes.1A.2A-heatplot.png",width=700,height=1+dim(genes.enriched.GO_mf)[[1]]*30)
heatplot(genes.enriched.GO_mf)+ theme(  axis.text.y = element_text( size = 12))
dev.off()

# ==========================================================================================================

genes   <- fread(glue("/home/jp/eclipse-workspace/Toulouse/data/results/{dir}/h4.genes.1B.2B.txt"),data.table=F)

# list to character
genes <- unlist(genes, use.names=FALSE)


h4.genes.1B.2B.dataframe <- as.data.frame(genes)
h4.genes.1B.2B.dataframe$geneslist  <- "genes.1B.2B"
h4.genes.1B.2B.dataframe$cellline  <- "h4"


genes.enriched.GO_bp    <- enricher(genes, TERM2GENE=GO_bpterm2gene, TERM2NAME=GO_bpterm2name,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)
genes.enriched.GO_mf    <- enricher(genes, TERM2GENE=GO_mfterm2gene, TERM2NAME=GO_mfterm2name,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)
genes.enriched.Reactome1 <- enricher(genes, TERM2GENE=Reactome, TERM2NAME=NA,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)


png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/Reactome.h4.genes.1B.2B.png",width=1200,height=1+dim(genes.enriched.Reactome)[[1]]*20)
barplot(genes.enriched.Reactome  ,  drop = TRUE,  showCategory = 50,  title = "Reactome", font.size = 14)+scale_fill_viridis(option="viridis",direction = -1)
dev.off()


write.csv(genes.enriched.Reactome  ,file = "/home/jp/eclipse-workspace/Toulouse/data/results/GO/Reactome.h4.genes.1B.2B.csv", row.names=FALSE)
png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/Reactome.h4.genes.1B.2B-heatplot.png",width=700,height=1+dim(genes.enriched.Reactome)[[1]]*30)
heatplot(genes.enriched.Reactome)+ theme(  axis.text.y = element_text( size = 12))
dev.off()

write.csv(genes.enriched.GO_bp   ,file = "/home/jp/eclipse-workspace/Toulouse/data/results/GO/BP.h4.genes.1B.2B.csv", row.names=FALSE)
png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/BP.h4.genes.1B.2B.png",width=900,height=1+dim(genes.enriched.GO_bp)[[1]]*5)
barplot(genes.enriched.GO_bp    ,  drop = TRUE,  showCategory = 50,  title = "Go-BP", font.size = 14)+scale_fill_viridis(option="viridis",direction = -1)
dev.off()

png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/BP.h4.genes.1B.2B-heatplot.png",width=700,height=1+dim(genes.enriched.GO_bp)[[1]]*30)
heatplot(genes.enriched.GO_bp)+ theme(  axis.text.y = element_text( size = 12))
dev.off()

if (length(genes.enriched.GO_mf) >= 1) {
  write.csv(genes.enriched.GO_mf  ,file = "/home/jp/eclipse-workspace/Toulouse/data/results/GO/MF.h4.genes.1B.2B.csv", row.names=FALSE)
  png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/MF.h4.genes.1B.2B.png",width=900,height=1+dim(genes.enriched.GO_mf)[[1]]*50)
  barplot(genes.enriched.GO_mf  ,  drop = TRUE,  showCategory = 50,  title = "Go-MF", font.size = 12)+scale_fill_viridis(option="viridis",direction = -1)
  dev.off()

png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/MF.h4.genes.1B.2B-heatplot.png",width=700,height=1+dim(genes.enriched.GO_mf)[[1]]*30)
heatplot(genes.enriched.GO_mf)+ theme(  axis.text.y = element_text( size = 12))
dev.off()

}



# ==========================================================================================================
# ================ PC9
# ==========================================================================================================


genes   <- fread(glue("/home/jp/eclipse-workspace/Toulouse/data/results/{dir}/pc9.genes.1A.2A.txt"),data.table=F)


# list to character
genes <- unlist(genes, use.names=FALSE)

pc9.genes.1A.2A.dataframe <- as.data.frame(genes)
pc9.genes.1A.2A.dataframe$geneslist  <- "genes.1A.2A"
pc9.genes.1A.2A.dataframe$cellline  <- "pc9"


genes.enriched.GO_bp <- enricher(genes, TERM2GENE=GO_bpterm2gene, TERM2NAME=GO_bpterm2name,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)
genes.enriched.GO_mf <- enricher(genes, TERM2GENE=GO_mfterm2gene, TERM2NAME=GO_mfterm2name,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)
genes.enriched.Reactome <- enricher(genes, TERM2GENE=Reactome, TERM2NAME=NA,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)


write.csv(genes.enriched.Reactome  ,file = "/home/jp/eclipse-workspace/Toulouse/data/results/GO/Reactome.pc9.genes.1A.2A.csv", row.names=FALSE)
png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/Reactome.pc9.genes.1A.2A.png",width=900,height=1+dim(genes.enriched.Reactome)[[1]]*30)
barplot(genes.enriched.Reactome  ,  drop = TRUE,  showCategory = 50,  title = "Reactome", font.size = 16) +scale_fill_viridis(option="viridis",direction = -1)
dev.off()

png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/Reactome.pc9.genes.1A.2A-heatplot.png",width=900,height=1+dim(genes.enriched.Reactome)[[1]]*50)
heatplot(genes.enriched.Reactome)+ theme(  axis.text.y = element_text( size = 12))
dev.off()

write.csv(genes.enriched.GO_bp   ,file = "/home/jp/eclipse-workspace/Toulouse/data/results/GO/BP.pc9.genes.1A.2A.csv", row.names=FALSE)
png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/BP.pc9.genes.1A.2A.png",width=900,height=1+dim(genes.enriched.GO_bp)[[1]]*7)
barplot(genes.enriched.GO_bp    ,  drop = TRUE,  showCategory = 50,  title = "Go-BP", font.size = 14)+scale_fill_viridis(option="viridis",direction = -1)
dev.off()

png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/BP.pc9.genes.1A.2A.csv-heatplot.png",width=900,height=1+dim(genes.enriched.GO_bp)[[1]]*50)
heatplot(genes.enriched.GO_bp) + theme(  axis.text.y = element_text( size = 12))
dev.off()

write.csv(genes.enriched.GO_mf  ,file = "/home/jp/eclipse-workspace/Toulouse/data/results/GO/MF.pc9.genes.1A.2A.csv", row.names=FALSE)
png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/MF.pc9.genes.1A.2A.png",width=900,height=1+dim(genes.enriched.GO_mf)[[1]]*50)
barplot(genes.enriched.GO_mf  ,  drop = TRUE,  showCategory = 50,  title = "Go-MF", font.size = 12)+scale_fill_viridis(option="viridis",direction = -1)
dev.off()

png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/MF.pc9.genes.1A.2A.csv-heatplot.png",width=900,height=1+dim(genes.enriched.GO_mf)[[1]]*50)
heatplot(genes.enriched.GO_mf)+ theme(  axis.text.y = element_text( size = 12))
dev.off()


# ==========================================================================================================

genes   <- fread(glue("/home/jp/eclipse-workspace/Toulouse/data/results/{dir}/pc9.genes.1B.2B.txt"),data.table=F)


# list to character
genes <- unlist(genes, use.names=FALSE)

pc9.genes.1B.2B.dataframe <- as.data.frame(genes)
pc9.genes.1B.2B.dataframe$geneslist  <- "genes.1B.2B"
pc9.genes.1B.2B.dataframe$cellline  <- "pc9"

genes.enriched.GO_bp <- enricher(genes, TERM2GENE=GO_bpterm2gene, TERM2NAME=GO_bpterm2name,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)
genes.enriched.GO_mf <- enricher(genes, TERM2GENE=GO_mfterm2gene, TERM2NAME=GO_mfterm2name,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)
genes.enriched.Reactome <- enricher(genes, TERM2GENE=Reactome, TERM2NAME=NA,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)


if (dim(genes.enriched.Reactome)[[1]] >= 1) {
write.csv(genes.enriched.Reactome  ,file = "/home/jp/eclipse-workspace/Toulouse/data/results/GO/Reactome.pc9.genes.1B.2B.csv", row.names=FALSE)
png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/Reactome.pc9.genes.1B.2B.png",width=700,height=1+dim(genes.enriched.Reactome)[[1]]*40)
barplot(genes.enriched.Reactome  ,  drop = TRUE,  showCategory = 50,  title = "Reactome", font.size = 12)+scale_fill_viridis(option="viridis",direction = -1)
dev.off()

png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/Reactome.pc9.genes.1B.2B.csv-heatplot.png",width=900,height=1+dim(genes.enriched.Reactome)[[1]]*50)
heatplot(genes.enriched.Reactome)+ theme(  axis.text.y = element_text( size = 12))
dev.off()
}
if (dim(genes.enriched.GO_bp)[[1]] >= 1) {
write.csv(genes.enriched.GO_bp   ,file = "/home/jp/eclipse-workspace/Toulouse/data/results/GO/BP.pc9.genes.1B.2B.csv", row.names=FALSE)
png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/BP.pc9.genes.1B.2B.png",width=900,height=1+dim(genes.enriched.GO_bp)[[1]]*8)
barplot(genes.enriched.GO_bp    ,  drop = TRUE,  showCategory = 50,  title = "Go-BP", font.size = 13)+scale_fill_viridis(option="viridis",direction = -1)
dev.off()

png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/BP.pc9.genes.1B.2B-heatplot.png",width=900,height=1+dim(genes.enriched.GO_bp)[[1]]*50)
heatplot(genes.enriched.GO_bp)+ theme(  axis.text.y = element_text( size = 12))
dev.off()
}

if (dim(genes.enriched.GO_mf)[[1]] >= 1) {
write.csv(genes.enriched.GO_mf  ,file = "/home/jp/eclipse-workspace/Toulouse/data/results/GO/MF.pc9.genes.1B.2B.csv", row.names=FALSE)
png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/MF.pc9.genes.1B.2B.png",width=700,height=1+dim(genes.enriched.GO_mf)[[1]]*30)
barplot(genes.enriched.GO_mf  ,  drop = TRUE,  showCategory = 50,  title = "Go-MF", font.size = 12)+scale_fill_viridis(option="viridis",direction = -1)
dev.off()

png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/MF.pc9.genes.1B.2B-heatplot.png",width=900,height=1+dim(genes.enriched.GO_mf)[[1]]*50)
heatplot(genes.enriched.GO_mf)+ theme(  axis.text.y = element_text( size = 12))
dev.off()
}

# ==========================================================================================================
# Dotplot final
# ==========================================================================================================

# Concat everything
pc9.dataframe <- rbind(pc9.genes.1A.2A.dataframe,pc9.genes.1B.2B.dataframe)

h4.dataframe <- rbind(h4.genes.1A.2A.dataframe,h4.genes.1B.2B.dataframe)

all.celllines.dataframe <- rbind(pc9.dataframe,h4.dataframe)

#Bug with points
all.celllines.dataframe$geneslist <- gsub("\\.", "_", all.celllines.dataframe$geneslist)

formula_res <- compareCluster(data=all.celllines.dataframe,genes~geneslist+cellline,  fun =  "enricher",TERM2GENE=Reactome, TERM2NAME=NA,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)

png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/Reactome.genes.doptplot.png",width=1200,height=1000)
dotplot(formula_res, x=~geneslist,showCategory=10)   + ggplot2::facet_grid(~cellline) + theme_gray( )+theme( axis.text.y = element_text( size = 16)) 
dev.off()

formula_res <- compareCluster(data=all.celllines.dataframe,genes~geneslist+cellline,  fun =  "enricher", TERM2GENE=GO_bpterm2gene, TERM2NAME=GO_bpterm2name,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)

png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/Gobp.genes.doptplot.png",width=1200,height=1000)
dotplot(formula_res, x=~geneslist,showCategory=10)   + ggplot2::facet_grid(~cellline) + theme_gray() +theme( axis.text.y = element_text( size = 16))  
dev.off()


formula_res <- compareCluster(data=all.celllines.dataframe,genes~geneslist+cellline,  fun =  "enricher", TERM2GENE=GO_mfterm2gene, TERM2NAME=GO_mfterm2name,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)

png(file="/home/jp/eclipse-workspace/Toulouse/data/results/GO/Gomf.genes.doptplot.png",width=1200,height=1000)
dotplot(formula_res, x=~geneslist,showCategory=10)   + ggplot2::facet_grid(~cellline) +   theme_gray() +theme( axis.text.y = element_text( size = 16)) 
dev.off()


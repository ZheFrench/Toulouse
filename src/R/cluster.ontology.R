#################################################################
#
# date: April 02, 2021
# platform: Ubuntu 10.04
# R.version : 4.0.3
# author: Villemin Jean-Philippe
# team: Bioinformatique et biologie des systèmes du cancer : J. Colinge 
# Institute : IRCM
#
# cluster.ontology.R
# Usage : 
# 
# cluster.ontology.R
# 
# Description : 
#
# Plot dotplots.
#
#
#################################################################

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library (glue))
suppressPackageStartupMessages(library(ggplot2))

# ==========   Settings =============================
base.dir <- "/data/villemin/code/Toulouse-rnaseq/"

GO_bp    <- read.gmt("/data/villemin/annotation/gsea/GO/Human_GO_bp_no_GO_iea_symbol.gmt")
GO_mf    <- read.gmt("/data/villemin/annotation/gsea/GO/Human_GO_mf_no_GO_iea_symbol.gmt")
GO_cc    <- read.gmt("/data/villemin/annotation/gsea/GO/Human_GO_cc_no_GO_iea_symbol.gmt")
Reactome <- read.gmt("/data/villemin/annotation/gsea/Reactome/ReactomePathways.25022021.gmt")
GoAll    <- read.gmt("/data/villemin/annotation/gsea/MSigDB/c5.go.v7.2.symbols.gmt")

Reactome <- Reactome[, c("term", "gene")]
GoAll    <- GoAll[, c("term", "gene")]

GO_bp <- GO_bp %>% separate(term, c("Name","Goterm"), sep = "%GOBP%")
GO_bpterm2gene <- GO_bp[, c("Goterm", "gene")]
GO_bpterm2name <- GO_bp[, c("Goterm", "Name")]

GO_mf <- GO_mf %>% separate(term, c("Name","Goterm"), sep = "%GOMF%")
GO_mfterm2gene <- GO_mf[, c("Goterm", "gene")]
GO_mfterm2name <- GO_mf[, c("Goterm", "Name")]

GO_cc <- GO_cc%>% separate(term, c("Name","Goterm"), sep = "%GOCC%")
GO_ccterm2gene <- GO_cc[, c("Goterm", "gene")]
GO_ccterm2name <- GO_cc[, c("Goterm", "Name")]

dir <- "genes_2_2"

# ==========================================================================================================
# Dotplot final
# ==========================================================================================================
all.celllines= list()

i=1
#for (cellLine in c("pc9","h4") ) {
  
# for (subset in c("1A","1B","2A","2B","3A","3B") ) {
  
    #genes   <- fread(glue("{base.dir}data/results/{dir}/{cellLine}.genes.pattern{subset}.txt"),data.table=F)
    
    # list to character
    #genes <- unlist(genes, use.names=FALSE)
    
    #genes.dataframe <- as.data.frame(genes)
    
    #genes.dataframe$geneslist  <- subset
    #genes.dataframe$cellline   <- cellLine
    ##all.celllines[[i]] <-  genes.dataframe # add it to your list
    #i=i+1
    #print(i)
  #}
#}


for (cellLine in c("pc9.dormance","pc9.expansion","h4.expansion","h4.dormance","h3") ) {
  
  for (subset in c("UP","DOWN") ) {
    
    genes   <- fread(glue("{base.dir}data/results/{dir}/{cellLine}.genes.pattern{subset}.txt"),data.table=F)
    
    # list to character
    genes <- unlist(genes, use.names=FALSE)
    
    genes.dataframe <- as.data.frame(genes)
    
    genes.dataframe$geneslist  <- subset
    genes.dataframe$cellline   <- cellLine
    all.celllines[[i]] <-  genes.dataframe # add it to your list
    i=i+1
    print(i)
  }
}

all.celllines.in.dataframe  = do.call(rbind, all.celllines)

head(all.celllines.in.dataframe )


write.csv(all.celllines.in.dataframe  ,file = glue("{base.dir}data/results/GO/alldatrame.csv"), row.names=FALSE)

#############################################################################
subpath <- file.path(glue("{base.dir}data/results/GO/"), "dotplot")
dir.create(subpath, recursive = TRUE)

#Bug with points
#all.celllines.dataframe$geneslist <- gsub("\\.", "_", all.celllines.in.dataframe$geneslist)

formula_res <- compareCluster(data=all.celllines.in.dataframe ,genes~geneslist+cellline,  fun =  "enricher",TERM2GENE=Reactome, TERM2NAME=NA,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)

png(file=glue("{subpath}/Reactome.genes.doptplot.png"),width=1200,height=1200)
dotplot(formula_res, x=~geneslist,showCategory=10)   + ggplot2::facet_grid(~cellline) + theme_gray( )+theme( axis.text.y = element_text( size = 16)) 
dev.off()

formula_res <- compareCluster(data=all.celllines.in.dataframe,genes~geneslist+cellline,  fun =  "enricher", TERM2GENE=GO_bpterm2gene, TERM2NAME=GO_bpterm2name,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)

png(file=glue("{subpath}/Gobp.genes.doptplot.png"),width=1600,height=1400)
dotplot(formula_res, x=~geneslist,showCategory=10)   + ggplot2::facet_grid(~cellline) + theme_gray() +theme( axis.text.y = element_text( size = 16))  
dev.off()

formula_res <- compareCluster(data=all.celllines.in.dataframe,genes~geneslist+cellline,  fun =  "enricher", TERM2GENE=GO_mfterm2gene, TERM2NAME=GO_mfterm2name,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)

png(file=glue("{subpath}/Gomf.genes.doptplot.png"),width=1600,height=1400)
dotplot(formula_res, x=~geneslist,showCategory=10)   + ggplot2::facet_grid(~cellline) + theme_gray() +theme( axis.text.y = element_text( size = 16))  
dev.off()

formula_res <- compareCluster(data=all.celllines.in.dataframe,genes~geneslist+cellline,  fun =  "enricher", TERM2GENE=GO_ccterm2gene, TERM2NAME=GO_ccterm2name,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)

png(file=glue("{subpath}/Gocc.genes.doptplot.png"),width=1600,height=1400)
dotplot(formula_res, x=~geneslist,showCategory=10)   + ggplot2::facet_grid(~cellline) + theme_gray() +theme( axis.text.y = element_text( size = 16))  
dev.off()

formula_res <- compareCluster(data=all.celllines.in.dataframe,genes~geneslist+cellline,  fun =  "enricher", TERM2GENE=GoAll, TERM2NAME=NA,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)

png(file=glue("{subpath}/GoAll.genes.doptplot.png"),width=1600,height=1400)
dotplot(formula_res, x=~geneslist,showCategory=10)   + ggplot2::facet_grid(~cellline) + theme_gray() +theme( axis.text.y = element_text( size = 16))  
dev.off()


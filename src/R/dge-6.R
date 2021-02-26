library(data.table)
library(edgeR)
#statmod
counts.1 <- fread("../../data/detc-counts.txt",data.table=F)
names(counts.1)[-(1:5)] <- gsub(names(counts.1)[-(1:5)],pattern="_S..?_align",replacement="",perl=T)
names(counts.1)[-(1:5)] <- gsub(names(counts.1)[-(1:5)],pattern="_21j",replacement="_21d",perl=T)
names(counts.1)[-(1:5)] <- gsub(names(counts.1)[-(1:5)],pattern="_75_jours",replacement="_75d",perl=T)

counts.2 <- counts.1[counts.1$type=="protein_coding",]
symbols <- counts.2$symbol
d <- symbols[duplicated(symbols)]
counts <- counts.2[!duplicated(counts.2$symbol),-(1:5)]
rownames(counts) <- symbols[!duplicated(symbols)]

q <- colSums(counts)
pdf("../../data/results/total-reads.pdf",width=5,height=5,pointsize=8,useDingbats=F)
par(mar=c(4,9,2,3)+0.1)
barplot(q,xlab="Read counts",horiz=T,las=2,cex.names=0.6,cex.axis=0.7)
dev.off()

good <- apply(counts,1,function(x) sum(x>5))>=3
sum(good)
counts <- counts[good,]

q <- apply(counts,2,function(x) quantile(x[x>0],prob=0.75))
ncounts <- sweep(counts,2,q/median(q),"/")

condition <- gsub(names(counts.1)[-(1:5)],pattern="[123]$",replacement="",perl=T)
condition <- gsub(condition,pattern="_$",replacement="",perl=T)
condition <- gsub(condition,pattern="_[ABDE]_",replacement="_",perl=T)
condition <- gsub(condition,pattern="-",replacement="_",perl=T)

conditions <- condition

names(conditions) <- names(ncounts)
conditions[conditions%in%c("PC9_3_NT1","PC9_CT_24h")] <- "PC9_CT"

#print (names(conditions))
#Note: it is prefered in R that the first level of a factor be the reference level
#(e.g. control, or untreated samples), so we can relevel the dex factor like so:

dge2by2 <- function  (two.cond,prefix,max.qval=1,min.FC=log2(1),min.count=0,plot=TRUE){
  
  comp <- paste(two.cond,collapse="-")
  fn <- paste0(prefix,"-",comp)
  mat <- cbind(ncounts[,conditions%in%two.cond[1]],ncounts[,conditions%in%two.cond[2]])
  print(head(mat))
  
  dge <- DGEList(mat,genes=rownames(ncounts))
  rownames(dge$counts) <- rownames(ncounts)

  dge$group     <- factor(conditions[names(mat)])
  
  print(conditions[conditions%in%c(two.cond[1])][[1]] )
  
  dge$group <-relevel(dge$group,ref=conditions[conditions%in%c(two.cond[1])][[1]])

  design <- model.matrix(~0+dge$group) # no intercept #x0 = 1,  force model throught the origin
  colnames(design) <- gsub("^dge$group","",colnames(design))
 
  cm <- makeContrasts(contrasts=comp,levels=dge$group)
  
  y     <- estimateDisp(dge,design,robust=T)
  fit.y <- glmFit(y,design)
  lrt   <- glmLRT(fit.y,contrast=cm)
  
  if (plot)
  
  # no filter
  sel.r <- topTags(lrt,p.value=1,adjust.method="BH",n=nrow(y$counts))
  #write.table(sel.r$table,file=paste0(fn,"-fgsea-gene-ratios.txt"),quote=F,row.names=F,sep="\t")
  
  #sel.r <- topTags(lrt,p.value=max.qval,adjust.method="BH",n=nrow(y$counts))
  # Filtering 
  sel.r <- sel.r[sel.r$table$FDR<=max.qval & abs(sel.r$table[[2]])>=min.FC & apply(dge$counts[rownames(sel.r),],1,max)>=min.count,]
  
  sel.rows<- rownames(sel.r)
  sel     <- cbind(sel.r$table,dge$counts[sel.rows,])
  sel.twt <- sel[order(sel$PValue),]
  
  write.table(sel.twt,file=paste0(fn,"-differential-gene-selection.txt"),quote=F,row.names=F,sep="\t")
  #write.table(sel.r$table,file=paste0(fn,"-differential-gene-ratios.txt"),quote=F,row.names=F,sep="\t")
  

} # dge2by2

#"PC9_Erlo_24h-PC9_CT","PC9_Erlo_DTC-PC9_Erlo_24h","PC9_3_DTEC_75d-PC9_Erlo_DTC",
#"H4006_Erlo_24h-H4006_CT_24h","H4006_Erlo_DTC-H4006_Erlo_24h","H4006_Erlo_DTEC-H4006_Erlo_DTC",
#"H3255_Erlo_24h-H3255_NT","H3255_Erlo_21d-H3255_Erlo_24h")
#JP-PC9_Erlo_DTC-PC9_CT-differential-gene-selection.txt

# Manually write each comp
comparisons <- c("PC9_Erlo_DTC-PC9_CT",
"H4006_Erlo_24h-H4006_CT_24h","H3255_Erlo_24h-H3255_NT",
"PC9_Erlo_DTC-PC9_CT","PC9_Erlo_DTC-PC9_Erlo_24h","H4006_Erlo_DTC-H4006_Erlo_24h",
"PC9_3_DTEC_75d-PC9_Erlo_DTC",
"H4006_Erlo_DTC-H4006_CT_24h","H3255_Erlo_21d-H3255_Erlo_24h",
"H4006_Erlo_DTEC-H4006_Erlo_DTC",
"H3255_Erlo_21d-H3255_NT")

for (comp in comparisons){
  two.cond <- strsplit(comp,"-")[[1]][1:2]
  dge2by2(two.cond,"../../data/results/JP")
}

conditions[conditions%in%c("PC9_3_NT1","PC9_CT_24h","H4006_CT_24h")] <- "Primary_CT"
conditions[conditions%in%c("PC9_Erlo_DTC","H4006_Erlo_DTC")] <- "Primary_DTC"

comparisons <- c("Primary_CT-H3255_NT","Primary_DTC-H3255_Erlo_21d")

for (comp in comparisons){
  two.cond <- strsplit(comp,"-")[[1]][1:2]
  dge2by2(two.cond,"../../data/results/JP")
}

#sel.genes <- NULL
#for (comp in comparisons){
#  d <- fread(paste0("../../data/results/JP-",comp,"-differential-gene-selection.txt"),data.table=F)
#  sel.genes <- union(sel.genes,d[["genes"]])
#}
#write.table(sel.genes,file="../../data/results/JP-genes-selection.txt",quote=F,row.names=F,col.names=F)



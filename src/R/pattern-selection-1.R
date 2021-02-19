library(data.table)
library(gplots)
library(ComplexHeatmap)
library(gtools)

pc9.dormance <- fread("../../data/results/JP-PC9_Erlo_DTC-PC9_CT-differential-gene-selection.txt",data.table=F)
pc9.expansion   <- fread("../../data/results/JP-PC9_3_DTEC_75d-PC9_Erlo_DTC-differential-gene-selection.txt",data.table=F)

h4.dormance <- fread("../../data/results/JP-H4006_Erlo_DTC-H4006_CT_24h-differential-gene-selection.txt",data.table=F)
h4.expansion   <- fread("../../data/results/JP-H4006_Erlo_DTEC-H4006_Erlo_DTC-differential-gene-selection.txt",data.table=F)

# Setup Parameters
min.FC <- foldchange2logratio(3, base=2)
print(min.FC)
max.FDR <- 0.01
min.count <- 100
stable.FC <- foldchange2logratio(2, base=2)
print(stable.FC)

#==================================================================================
#================================= Functions ======================================
#==================================================================================

genes_up <- function (list,FC=min.FC,FDR=max.FDR,count=min.count){ 
  
    list.genes.up <- list[
    list$logFC >= FC 
    & list$FDR<=FDR 
    & apply(list[,-(1:6)],1,max)>count
    ,"genes" ]
  length(list.genes.up)
    
  return(list.genes.up)
}

genes_down <- function (list,FC=min.FC,FDR=max.FDR,count=min.count){ 
  
  list.genes.down <- list[
    list$logFC <= - FC 
    & list$FDR<=FDR 
    & apply(list[,-(1:6)],1,max)>count
    ,"genes" ]
  
  length(list.genes.down)
  return(list.genes.down)
}

genes_stable<- function (list,FC=stable.FC,FDR=max.FDR,count=min.count){ 
  
  list.genes.stable <- list[
    list$logFC >= - FC  &  list$logFC <= FC 
    & apply(list[,-(1:6)],1,max)>count
    ,"genes" ]
  
  length(list.genes.stable)
  return(list.genes.stable)
}

# Patterns Definition :
# 1 - A : Up & Down,  B : Down & Up
# 2 - A : Up & Stable,  B : Down & Stable
# 3 - A : Stable & Up,  B : Stable & Down
# 4 - A : Up & Up,  B : Down & Down

#==================================================================================
#==================================================================================

#--------------- Pattern 1 -----------------------
print("===============>Pattern 1")
# PC9 Dormance -> UP and Expansion -> Down
pc9.genes.pattern1A <- intersect(genes_up(pc9.dormance),genes_down(pc9.expansion))
print("intersection pc9.1A")
length(pc9.genes.pattern1A)

# PC9 Dormance -> DOWN and Expansion -> UP
pc9.genes.pattern1B <- intersect(genes_down(pc9.dormance),genes_up(pc9.expansion))
print("intersection pc9.1B ")
length(pc9.genes.pattern1B)

# h4 Dormance -> UP and Expansion -> DOWN
h4.genes.pattern1A <- intersect(genes_up(h4.dormance),genes_down(h4.expansion))
print("intersection h4.1A")
length(h4.genes.pattern1A)

# h4 Dormance -> DOWN and Expansion -> UP
h4.genes.pattern1B <- intersect(genes_down(h4.dormance),genes_up(h4.expansion))
print("intersection h4.1B ")
length(h4.genes.pattern1B)

#--------------- Pattern 2 -----------------------
print("===============>Pattern 2")

# PC9 Dormance -> UP and Expansion -> STABLE
pc9.genes.pattern2A <- intersect(genes_up(pc9.dormance),genes_stable(pc9.expansion))
print("intersection pc9.2A")
length(pc9.genes.pattern2A)

# PC9 Dormance -> DOWN and Expansion -> STABLE
pc9.genes.pattern2B <- intersect(genes_down(pc9.dormance),genes_stable(pc9.expansion))
print("intersection pc9.2B ")
length(pc9.genes.pattern2B)

# h4 Dormance -> UP and Expansion -> STABLE
h4.genes.pattern2A <- intersect(genes_up(h4.dormance),genes_stable(h4.expansion))
print("intersection h4.2A")
length(h4.genes.pattern2A)

# h4 Dormance -> DOWN and Expansion -> STABLE
h4.genes.pattern2B <- intersect(genes_down(h4.dormance),genes_stable(h4.expansion))
print("intersection h4.2B ")
length(h4.genes.pattern2B)


#--------------- Pattern 3 -----------------------
print("===============>Pattern 3")

# PC9 Dormance -> STABLE and Expansion -> UP
pc9.genes.pattern3A <- intersect(genes_stable(pc9.dormance),genes_up(pc9.expansion))
print("intersection pc9.3A")
length(pc9.genes.pattern3A)

# PC9 Dormance -> STABLE and Expansion -> DOWN
pc9.genes.pattern3B <- intersect(genes_stable(pc9.dormance),genes_down(pc9.expansion))
print("intersection pc9.3B ")
length(pc9.genes.pattern3B)

# h4 Dormance -> STABLE and Expansion -> UP
h4.genes.pattern3A <- intersect(genes_stable(h4.dormance),genes_up(h4.expansion))
print("intersection h4.3A")
length(h4.genes.pattern3A)

# h4 Dormance -> STABLE and Expansion -> DOWN
h4.genes.pattern3B <- intersect(genes_stable(h4.dormance),genes_down(h4.expansion))
print("intersection h4.3B ")
length(h4.genes.pattern3B)

#--------------- Pattern 4 -----------------------
print("===============>Pattern 4")

# PC9 Dormance -> UP and Expansion -> UP
pc9.genes.pattern4A <- intersect(genes_up(pc9.dormance),genes_up(pc9.expansion))
print("intersection pc9.4A")
length(pc9.genes.pattern4A)

# PC9 Dormance -> DOWN and Expansion -> DOWN
pc9.genes.pattern4B <- intersect(genes_down(pc9.dormance),genes_down(pc9.expansion))
print("intersection pc9.4B ")
length(pc9.genes.pattern4B)

# h4 Dormance -> UP and Expansion -> UP
h4.genes.pattern4A <- intersect(genes_up(h4.dormance),genes_up(h4.expansion))
print("intersection h4.4A")
length(h4.genes.pattern4A)

# h4 Dormance -> DOWN and Expansion -> DOWN
h4.genes.pattern4B <- intersect(genes_down(h4.dormance),genes_down(h4.expansion))
print("intersection h4.4B ")
length(h4.genes.pattern4B)
head(h4.genes.pattern4B)

#-------------------------    Intersection     ----------------------------
print("===============> Intersection ")

print("intersection pattern1A h4 and pc9 ")
length(intersect(h4.genes.pattern1A,pc9.genes.pattern1A))

print("intersection pattern1B h4 and pc9 ")
length(intersect(h4.genes.pattern1B,pc9.genes.pattern1B))

print("intersection pattern2A h4 and pc9 ")
length(intersect(h4.genes.pattern2A,pc9.genes.pattern2A))

print("intersection pattern2B h4 and pc9 ")
length(intersect(h4.genes.pattern2B,pc9.genes.pattern2B))

print("intersection pattern3A h4 and pc9 ")
length(intersect(h4.genes.pattern3A,pc9.genes.pattern3A))

print("intersection pattern3B h4 and pc9 ")
length(intersect(h4.genes.pattern3B,pc9.genes.pattern3B))

print("intersection pattern4A h4 and pc9 ")
length(intersect(h4.genes.pattern4A,pc9.genes.pattern4A))

print("intersection pattern4B h4 and pc9 ")
length(intersect(h4.genes.pattern4B,pc9.genes.pattern4B))

# Patterns Definition :
# 1 - A : Up & Down,  B : Down & Up
# 2 - A : Up & Stable,  B : Down & Stable
# 3 - A : Stable & Up,  B : Stable & Down
# 4 - A : Up & Up,  B : Down & Down

stop()
 
# ====================================================================================
# control plots
# ====================================================================================

t <- table(c(pc9.24h$genes,pc9.dtc$genes,pc9.dtec$genes,h3.24h$genes,h3.21d$genes))
common.pc9 <- names(t)[t==5]
init.pc9.fc <- setNames(pc9.24h$logFC,pc9.24h$genes)
init.h3.fc <- setNames(h3.24h$logFC,h3.24h$genes)
i21d.h3.fc <- setNames(h3.21d$logFC,h3.21d$genes)
pc9.ratios <- cbind(init.pc9.fc[common.pc9],dtc.pc9.fc[common.pc9],dtec.pc9.fc[common.pc9],init.h3.fc[common.pc9],i21d.h3.fc[common.pc9])
colnames(pc9.ratios) <- c("PC9 24h","PC9 DTC","PC9 DTEC","H3255 24h","H3255 21d")

Heatmap(pc9.ratios[intersect(no.init.opposed.24h.dtec.pc9,common.pc9),],show_column_dend=F,column_order=colnames(pc9.ratios))
Heatmap(pc9.ratios[intersect(no.later.opposed.dtc.dtec.pc9,common.pc9),],show_column_dend=F,column_order=colnames(pc9.ratios))
Heatmap(pc9.ratios[intersect(opposed.init.opposed.24h.dtec.pc9,common.pc9),],show_column_dend=F,column_order=colnames(pc9.ratios))
Heatmap(pc9.ratios[intersect(opposed.later.opposed.dtc.dtec.pc9,common.pc9),],show_column_dend=F,column_order=colnames(pc9.ratios))
Heatmap(pc9.ratios[intersect(opposed.21d.opposed.24h.dtc.pc9,common.pc9),],show_column_dend=F,column_order=colnames(pc9.ratios))



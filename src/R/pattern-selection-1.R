#################################################################
#
# date: April 02, 2021
# platform: Ubuntu 10.04
# R.version : 4.0.3
# author: Villemin Jean-Philippe
# team: Bioinformatique et biologie des systèmes du cancer : J. Colinge 
# Institute : IRCM
#
# pattern-selection-1.R
# Usage : 
# 
# pattern-selection-1.R
# 
# Description : 
#
# Plot dotplots.
#
#################################################################

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(glue))

# ==========   Settings =============================

dir.path        <- "/data/villemin/code/Toulouse-rnaseq/data/results/"
pc9.dormance    <- fread(glue("{dir.path}/pc9/JP-PC9_Erlo_DTC-PC9_CT-differential-gene-selection.txt"),data.table=F)
pc9.expansion   <- fread(glue("{dir.path}/pc9/JP-PC9_3_DTEC_75d-PC9_Erlo_DTC-differential-gene-selection.txt"),data.table=F)

h4.dormance     <- fread(glue("{dir.path}/h4/JP-H4006_Erlo_DTC-H4006_CT_24h-differential-gene-selection.txt"),data.table=F)
h4.expansion    <- fread(glue("{dir.path}/h4/JP-H4006_Erlo_DTEC-H4006_Erlo_DTC-differential-gene-selection.txt"),data.table=F)

h3.control      <- fread(glue("{dir.path}/h4/JP-H3255_Erlo_21d-H3255_NT-differential-gene-selection.txt"),data.table=F)

my.genes        <- fread("/data/villemin/code/Toulouse-rnaseq/data/custom.geneslist.txt",data.table=F)

#print(my.genes)

# Setup Parameters
linear.FC = 2
linear.stable_FC = 2

min.FC <- foldchange2logratio(linear.FC, base=2)
stable.FC <- foldchange2logratio(linear.stable_FC, base=2)

max.FDR <- 0.01
min.count <- 100

# Create dir based of linear FC
dir.create(file.path("{dir.path}/", glue("genes_{linear.FC}_{linear.stable_FC}")), showWarnings = FALSE)

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
    list$logFC > - FC  &  list$logFC <  FC 
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

write.table(sort(genes_up(h3.control)),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/h3.genes.patternUP.txt"),quote=F,row.names=F,sep="\t",col.names=F)
write.table(sort(genes_down(h3.control)),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/h3.genes.patternDOWN.txt"),quote=F,row.names=F,sep="\t",col.names=F)

cat(genes_up(h3.control),genes_down(h3.control),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/h3.control.genes.patternALL.txt"), sep="\n")

write.table(sort(genes_up(h4.dormance)),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/h4.dormance.genes.patternUP.txt"),quote=F,row.names=F,sep="\t",col.names=F)
write.table(sort(genes_down(h4.dormance)),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/h4.dormance.genes.patternDOWN.txt"),quote=F,row.names=F,sep="\t",col.names=F)

cat(genes_up(h4.dormance),genes_down(h4.dormance),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/h4.dormance.genes.patternALL.txt"), sep="\n")

write.table(sort(genes_up(h4.expansion)),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/h4.expansion.genes.patternUP.txt"),quote=F,row.names=F,sep="\t",col.names=F)
write.table(sort(genes_down(h4.expansion)),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/h4.expansion.genes.patternDOWN.txt"),quote=F,row.names=F,sep="\t",col.names=F)

cat(genes_up(h4.expansion),genes_down(h4.expansion),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/h4.expansion.genes.patternALL.txt"), sep="\n")

write.table(sort(genes_up(pc9.expansion)),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/pc9.expansion.genes.patternUP.txt"),quote=F,row.names=F,sep="\t",col.names=F)
write.table(sort(genes_down(pc9.expansion)),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/pc9.expansion.genes.patternDOWN.txt"),quote=F,row.names=F,sep="\t",col.names=F)

cat(genes_up(pc9.expansion),genes_down(pc9.expansion),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/pc9.expansion.genes.patternALL.txt"), sep="\n")

write.table(sort(genes_up(pc9.dormance)),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/pc9.dormance.genes.patternUP.txt"),quote=F,row.names=F,sep="\t",col.names=F)
write.table(sort(genes_down(pc9.dormance)),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/pc9.dormance.genes.patternDOWN.txt"),quote=F,row.names=F,sep="\t",col.names=F)

cat(genes_up(pc9.dormance),genes_down(pc9.dormance),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/pc9.dormance.genes.patternALL.txt"), sep="\n")

stop()
#--------------- Pattern 1 -----------------------
print("===============>Pattern 1")
# PC9 Dormance -> UP and Expansion -> Down
pc9.genes.pattern1A <- intersect(genes_up(pc9.dormance),genes_down(pc9.expansion))
print("intersection pc9.1A")
length(pc9.genes.pattern1A)

write.table(sort(pc9.genes.pattern1A),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/pc9.genes.pattern1A.txt"),quote=F,row.names=F,sep="\t",col.names=F)

# PC9 Dormance -> DOWN and Expansion -> UP
pc9.genes.pattern1B <- intersect(genes_down(pc9.dormance),genes_up(pc9.expansion))
print("intersection pc9.1B ")
length(pc9.genes.pattern1B)
write.table(sort(pc9.genes.pattern1B),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/pc9.genes.pattern1B.txt"),quote=F,row.names=F,sep="\t",col.names=F)

# h4 Dormance -> UP and Expansion -> DOWN
h4.genes.pattern1A <- intersect(genes_up(h4.dormance),genes_down(h4.expansion))
print("intersection h4.1A")
length(h4.genes.pattern1A)
write.table(sort(h4.genes.pattern1A),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/h4.genes.pattern1A.txt"),quote=F,row.names=F,sep="\t",col.names=F)

# h4 Dormance -> DOWN and Expansion -> UP
h4.genes.pattern1B <- intersect(genes_down(h4.dormance),genes_up(h4.expansion))
print("intersection h4.1B ")
length(h4.genes.pattern1B)
write.table(sort(h4.genes.pattern1B),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/h4.genes.pattern1B.txt"),quote=F,row.names=F,sep="\t",col.names=F)

#--------------- Pattern 2 -----------------------
print("===============>Pattern 2")

# PC9 Dormance -> UP and Expansion -> STABLE
pc9.genes.pattern2A <- intersect(genes_up(pc9.dormance),genes_stable(pc9.expansion))
print("intersection pc9.2A")
length(pc9.genes.pattern2A)
write.table(sort(pc9.genes.pattern2A),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/pc9.genes.pattern2A.txt"),quote=F,row.names=F,sep="\t",col.names=F)

# PC9 Dormance -> DOWN and Expansion -> STABLE
pc9.genes.pattern2B <- intersect(genes_down(pc9.dormance),genes_stable(pc9.expansion))
print("intersection pc9.2B ")
length(pc9.genes.pattern2B)
write.table(sort(pc9.genes.pattern2B),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/pc9.genes.pattern2B.txt"),quote=F,row.names=F,sep="\t",col.names=F)

# h4 Dormance -> UP and Expansion -> STABLE
h4.genes.pattern2A <- intersect(genes_up(h4.dormance),genes_stable(h4.expansion))
print("intersection h4.2A")
length(h4.genes.pattern2A)
write.table(sort(h4.genes.pattern2A),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/h4.genes.pattern2A.txt"),quote=F,row.names=F,sep="\t",col.names=F)

# h4 Dormance -> DOWN and Expansion -> STABLE
h4.genes.pattern2B <- intersect(genes_down(h4.dormance),genes_stable(h4.expansion))
print("intersection h4.2B ")
length(h4.genes.pattern2B)
write.table(sort(h4.genes.pattern2B),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/h4.genes.pattern2B.txt"),quote=F,row.names=F,sep="\t",col.names=F)


#--------------- Pattern 3 -----------------------
print("===============>Pattern 3")

# PC9 Dormance -> STABLE and Expansion -> UP
pc9.genes.pattern3A <- intersect(genes_stable(pc9.dormance),genes_up(pc9.expansion))
print("intersection pc9.3A")
length(pc9.genes.pattern3A)
write.table(sort(pc9.genes.pattern3A),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/pc9.genes.pattern3A.txt"),quote=F,row.names=F,sep="\t",col.names=F)

# PC9 Dormance -> STABLE and Expansion -> DOWN
pc9.genes.pattern3B <- intersect(genes_stable(pc9.dormance),genes_down(pc9.expansion))
print("intersection pc9.3B ")
length(pc9.genes.pattern3B)
write.table(sort(pc9.genes.pattern3B),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/pc9.genes.pattern3B.txt"),quote=F,row.names=F,sep="\t",col.names=F)

# h4 Dormance -> STABLE and Expansion -> UP
h4.genes.pattern3A <- intersect(genes_stable(h4.dormance),genes_up(h4.expansion))
print("intersection h4.3A")
length(h4.genes.pattern3A)
write.table(sort(h4.genes.pattern3A),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/h4.genes.pattern3A.txt"),quote=F,row.names=F,sep="\t",col.names=F)


# h4 Dormance -> STABLE and Expansion -> DOWN
h4.genes.pattern3B <- intersect(genes_stable(h4.dormance),genes_down(h4.expansion))
print("intersection h4.3B ")
length(h4.genes.pattern3B)
write.table(sort(h4.genes.pattern3B),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/h4.genes.pattern3B.txt"),quote=F,row.names=F,sep="\t",col.names=F)

#--------------- Pattern 4 -----------------------
print("===============>Pattern 4")

# PC9 Dormance -> UP and Expansion -> UP
pc9.genes.pattern4A <- intersect(genes_up(pc9.dormance),genes_up(pc9.expansion))
print("intersection pc9.4A")
length(pc9.genes.pattern4A)
write.table(sort(pc9.genes.pattern4A),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/pc9.genes.pattern4A.txt"),quote=F,row.names=F,sep="\t",col.names=F)

# PC9 Dormance -> DOWN and Expansion -> DOWN
pc9.genes.pattern4B <- intersect(genes_down(pc9.dormance),genes_down(pc9.expansion))
print("intersection pc9.4B ")
length(pc9.genes.pattern4B)
write.table(sort(pc9.genes.pattern4B),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/pc9.genes.pattern4B.txt"),quote=F,row.names=F,sep="\t",col.names=F)

# h4 Dormance -> UP and Expansion -> UP
h4.genes.pattern4A <- intersect(genes_up(h4.dormance),genes_up(h4.expansion))
print("intersection h4.4A")
length(h4.genes.pattern4A)
write.table(sort(h4.genes.pattern4A),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/h4.genes.pattern4A.txt"),quote=F,row.names=F,sep="\t",col.names=F)

# h4 Dormance -> DOWN and Expansion -> DOWN
h4.genes.pattern4B <- intersect(genes_down(h4.dormance),genes_down(h4.expansion))
print("intersection h4.4B ")
length(h4.genes.pattern4B)
head(h4.genes.pattern4B)
write.table(sort(h4.genes.pattern4B),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/h4.genes.pattern4B.txt"),quote=F,row.names=F,sep="\t",col.names=F)

#-------------------------    Intersection     ----------------------------
print("===============> Intersection ")

print("intersection pattern1A h4 and pc9 ")
length(intersect(h4.genes.pattern1A,pc9.genes.pattern1A))
intersect.pattern1A = intersect(h4.genes.pattern1A,pc9.genes.pattern1A)
print(intersect.pattern1A)
write.table(sort(intersect.pattern1A),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/intersect.pattern1A.txt"),quote=F,row.names=F,sep="\t",col.names=F)

print("intersection pattern1B h4 and pc9 ")
length(intersect(h4.genes.pattern1B,pc9.genes.pattern1B))
intersect.pattern1B <- intersect(h4.genes.pattern1B,pc9.genes.pattern1B)
print(intersect.pattern1B)
write.table(sort(intersect.pattern1B),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/intersect.pattern1B.txt"),quote=F,row.names=F,sep="\t",col.names=F)

print("intersection pattern2A h4 and pc9 ")
length(intersect(h4.genes.pattern2A,pc9.genes.pattern2A))
intersect.pattern2A <- intersect(h4.genes.pattern2A,pc9.genes.pattern2A)
print(intersect.pattern2A)
write.table(sort(intersect.pattern2A),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/intersect.pattern2A.txt"),quote=F,row.names=F,sep="\t",col.names=F)

print("intersection pattern2B h4 and pc9 ")
length(intersect(h4.genes.pattern2B,pc9.genes.pattern2B))
intersect.pattern2B <- intersect(h4.genes.pattern2B,pc9.genes.pattern2B)
print(intersect.pattern2B)
write.table(sort(intersect.pattern2B),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/intersect.pattern2B.txt"),quote=F,row.names=F,sep="\t",col.names=F)

print("intersection pattern3A h4 and pc9 ")
length(intersect(h4.genes.pattern3A,pc9.genes.pattern3A))
intersect.pattern3A<- intersect(h4.genes.pattern3A,pc9.genes.pattern3A)
write.table(sort(intersect.pattern3A),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/intersect.pattern3A.txt"),quote=F,row.names=F,sep="\t",col.names=F)

print("intersection pattern3B h4 and pc9 ")
length(intersect(h4.genes.pattern3B,pc9.genes.pattern3B))
intersect.pattern3B<- intersect(h4.genes.pattern3B,pc9.genes.pattern3B)
write.table(sort(intersect.pattern3B),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/intersect.pattern3B.txt"),quote=F,row.names=F,sep="\t",col.names=F)

print("intersection pattern4A h4 and pc9 ")
length(intersect(h4.genes.pattern4A,pc9.genes.pattern4A))
intersect.pattern4A<- intersect(h4.genes.pattern4A,pc9.genes.pattern4A)
write.table(sort(intersect.pattern4A),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/intersect.pattern4A.txt"),quote=F,row.names=F,sep="\t",col.names=F)

print("intersection pattern4B h4 and pc9 ")
length(intersect(h4.genes.pattern4B,pc9.genes.pattern4B))
intersect.pattern4B<- intersect(h4.genes.pattern4B,pc9.genes.pattern4B)
write.table(sort(intersect.pattern4B),file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/intersect.pattern4B.txt"),quote=F,row.names=F,sep="\t",col.names=F)

# Patterns Definition :
# 1 - A : Up & Down,  B : Down & Up
# 2 - A : Up & Stable,  B : Down & Stable
# 3 - A : Stable & Up,  B : Stable & Down
# 4 - A : Up & Up,  B : Down & Down

data.pc9 <- inner_join(select(pc9.dormance, -c(3:last_col())), select(pc9.expansion, -c(3:last_col())), by = c("genes"),keep=FALSE)
colnames(data.pc9) <- gsub('logFC.x', 'pc9.dormance', colnames(data.pc9), fixed=TRUE)
colnames(data.pc9) <- gsub('logFC.y', 'pc9.expansion', colnames(data.pc9), fixed=TRUE)

#print(filter(data.pc9,genes %in%  my.genes[[1]]))

data.h4 <- inner_join(select(h4.dormance, -c(3:last_col())), select(h4.expansion, -c(3:last_col())), by =c("genes"),keep=FALSE) 
colnames(data.h4) <- gsub('logFC.x', 'h4.dormance', colnames(data.h4), fixed=TRUE)
colnames(data.h4) <- gsub('logFC.y', 'h4.expansion', colnames(data.h4), fixed=TRUE)

data.h4.pc9 <- inner_join(data.pc9,data.h4, by =c("genes"),keep=FALSE)

data.h4.pc9.h3 <- inner_join(select(h3.control, -c(3:last_col())) , data.h4.pc9 , by =c("genes"),keep=FALSE) 
colnames(data.h4.pc9.h3) <- gsub('logFC', 'h3.control', colnames(data.h4.pc9.h3), fixed=TRUE)

#write.table(filter(data.h4.pc9.h3, genes%in%  my.genes),file=glue("{dir.path}/custom_matrice_heatmap_all_patterns.txt"),quote=F,row.names=F,sep="\t")

data.h4.pc9.h3$pattern <- 0
data.h4.pc9.h3 <- data.h4.pc9.h3 %>% select(pattern, everything())

data.h4.pc9.h3$pattern[data.h4.pc9.h3$genes %in% intersect.pattern1A] <- "1A"
data.h4.pc9.h3$pattern[data.h4.pc9.h3$genes %in% intersect.pattern1B] <- "1B"
data.h4.pc9.h3$pattern[data.h4.pc9.h3$genes %in% intersect.pattern2A] <- "2A"
data.h4.pc9.h3$pattern[data.h4.pc9.h3$genes %in% intersect.pattern2B] <- "2B"
data.h4.pc9.h3$pattern[data.h4.pc9.h3$genes %in% intersect.pattern3A] <- "3A"
data.h4.pc9.h3$pattern[data.h4.pc9.h3$genes %in% intersect.pattern3B] <- "3B"
#data.h4.pc9.h3<- data.h4.pc9.h3[data.h4.pc9.h3$pattern!=0,]

write.table(data.h4.pc9.h3[data.h4.pc9.h3$pattern!=0,],file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/matrice_heatmap_intersect_patterns.txt"),quote=F,row.names=F,sep="\t")


data.h4.pc9.h3$patternpc9 <- 0
data.h4.pc9.h3 <- data.h4.pc9.h3 %>% select(patternpc9, everything())

data.h4.pc9.h3$patternpc9[data.h4.pc9.h3$genes %in% pc9.genes.pattern1A] <- "1A"
data.h4.pc9.h3$patternpc9[data.h4.pc9.h3$genes %in% pc9.genes.pattern1B] <- "1B"
data.h4.pc9.h3$patternpc9[data.h4.pc9.h3$genes %in% pc9.genes.pattern2A] <- "2A"
data.h4.pc9.h3$patternpc9[data.h4.pc9.h3$genes %in% pc9.genes.pattern2B] <- "2B"
data.h4.pc9.h3$patternpc9[data.h4.pc9.h3$genes %in% pc9.genes.pattern3A] <- "3A"
data.h4.pc9.h3$patternpc9[data.h4.pc9.h3$genes %in% pc9.genes.pattern3B] <- "3B"
data.h4.pc9.h3$patternpc9[data.h4.pc9.h3$genes %in% pc9.genes.pattern4A] <- "4A"
data.h4.pc9.h3$patternpc9[data.h4.pc9.h3$genes %in% pc9.genes.pattern4B] <- "4B"
#data.h4.pc9.h3<- data.h4.pc9.h3[data.h4.pc9.h3$pattern!=0,]


#write.table(data.h4.pc9.h3,file=glue("{dir.path}/matrice_heatmap_pc9_patterns.txt"),quote=F,row.names=F,sep="\t")

data.h4.pc9.h3$patternh4 <- 0
data.h4.pc9.h3 <- data.h4.pc9.h3 %>% select(patternh4, everything())

data.h4.pc9.h3$patternh4[data.h4.pc9.h3$genes %in% h4.genes.pattern1A] <- "1A"
data.h4.pc9.h3$patternh4[data.h4.pc9.h3$genes %in% h4.genes.pattern1B] <- "1B"
data.h4.pc9.h3$patternh4[data.h4.pc9.h3$genes %in% h4.genes.pattern2A] <- "2A"
data.h4.pc9.h3$patternh4[data.h4.pc9.h3$genes %in% h4.genes.pattern2B] <- "2B"
data.h4.pc9.h3$patternh4[data.h4.pc9.h3$genes %in% h4.genes.pattern3A] <- "3A"
data.h4.pc9.h3$patternh4[data.h4.pc9.h3$genes %in% h4.genes.pattern3B] <- "3B"
data.h4.pc9.h3$patternh4[data.h4.pc9.h3$genes %in% h4.genes.pattern4A] <- "4A"
data.h4.pc9.h3$patternh4[data.h4.pc9.h3$genes %in% h4.genes.pattern4B] <- "4B"
#data.h4.pc9.h3<- data.h4.pc9.h3[data.h4.pc9.h3$pattern!=0,]

#write.table(data.h4.pc9.h3,file=glue("{dir.path}/matrice_heatmap_h4_patterns.txt"),quote=F,row.names=F,sep="\t")

data.h4.pc9.h3$contrast <- "h3.undef" 
data.h4.pc9.h3 <- data.h4.pc9.h3 %>% select(contrast, everything())

data.h4.pc9.h3$contrast[data.h4.pc9.h3$genes %in% genes_up(h3.control)] <- "h3.up"
data.h4.pc9.h3$contrast[data.h4.pc9.h3$genes %in% genes_down(h3.control)] <- "h3.down"
data.h4.pc9.h3$contrast[data.h4.pc9.h3$genes %in% genes_stable(h3.control)] <- "h3.stable"


data.h4.pc9.h3.pc9Contrast <- data.h4.pc9.h3[data.h4.pc9.h3$patternpc9!=0 & data.h4.pc9.h3$contrast !="h3.undef"  ,]
data.h4.pc9.h3.pc9Contrast <- select(data.h4.pc9.h3.pc9Contrast , -c(	h4.dormance,	h4.expansion,patternh4))

write.table(data.h4.pc9.h3.pc9Contrast,file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/matrice_heatmap_pc9_contrasted_patterns.txt"),quote=F,row.names=F,sep="\t")

data.h4.pc9.h3.h4Contrast <- data.h4.pc9.h3[data.h4.pc9.h3$patternh4!=0 & data.h4.pc9.h3$contrast !="h3.undef" ,]

data.h4.pc9.h3.h4Contrast <- select(data.h4.pc9.h3.h4Contrast , -c(	pc9.dormance,	pc9.expansion,patternpc9)  ) 

countspc9 <- rename(count(data.h4.pc9.h3.pc9Contrast,patternpc9, contrast), Freq = n)
countsh4  <- rename(count(data.h4.pc9.h3.h4Contrast,patternh4, contrast), Freq = n)
write.table(countspc9,file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/freq_pc9_contrasted_patterns.txt"),quote=F,row.names=F,sep="\t")
write.table(countsh4,file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/freq_h4_contrasted_patterns.txt"),quote=F,row.names=F,sep="\t")


write.table(data.h4.pc9.h3.h4Contrast,file=glue("{dir.path}/genes_{linear.FC}_{linear.stable_FC}/matrice_heatmap_h4_contrasted_patterns.txt"),quote=F,row.names=F,sep="\t")


#==================================================================================
#=====================         trash                 ==============================
#==================================================================================

#data.h4.pc9.h3 <-  unite(data.h4.pc9.h3,col="pc9_h4_pattern" ,patternpc9:patternh4, sep = "_", remove = FALSE, na.rm = FALSE)
#data.h4.pc9.h3 <-  unite(data.h4.pc9.h3,col="pc9_contrast" ,patternpc9:contrast, sep = "_", remove = FALSE, na.rm = FALSE)
#data.h4.pc9.h3 <-  unite(data.h4.pc9.h3,col="h4_contrast" ,patternh4:contrast, sep = "_", remove = FALSE, na.rm = FALSE)


#primary.vs.meta.untreated <- fread("{dir.path}/JP-Primary_CT-H3255_NT-differential-gene-selection.txt",data.table=F)
#primary.vs.meta.treated <- fread("{dir.path}/JP-Primary_DTC-H3255_Erlo_21d-differential-gene-selection.txt",data.table=F)

# Patterns Definition :
# 1 - Up & Up,  
# 2 - Up & Down, 
# 3 -  Down & Up, 
# 4 -Down & Down 

#primary.vs.meta.genes.pattern1 <- intersect(genes_up(primary.vs.meta.untreated ),genes_up(primary.vs.meta.treated))
#primary.vs.meta.genes.pattern2 <- intersect(genes_up(primary.vs.meta.untreated ),genes_down(primary.vs.meta.treated))
#primary.vs.meta.genes.pattern3 <- intersect(genes_down(primary.vs.meta.untreated ),genes_up(primary.vs.meta.treated))
#primary.vs.meta.genes.pattern4 <- intersect(genes_down(primary.vs.meta.untreated ),genes_down(primary.vs.meta.treated))

#print("primary.vs.meta.genes.pattern1 ")
#length(primary.vs.meta.genes.pattern1)
#print(primary.vs.meta.genes.pattern1)
#write.table(primary.vs.meta.genes.pattern1,file=glue("{dir.path}/primary.vs.meta.pattern1.txt"),quote=F,row.names=F,sep="\t")

#print("primary.vs.meta.genes.pattern2 ")
#length(primary.vs.meta.genes.pattern2)
#print(primary.vs.meta.genes.pattern2)

#print("primary.vs.meta.genes.pattern3 ")
#length(primary.vs.meta.genes.pattern3)
#print(primary.vs.meta.genes.pattern3)

#print("primary.vs.meta.genes.pattern4 ")
#length(primary.vs.meta.genes.pattern4)
#print(primary.vs.meta.genes.pattern4)
#write.table(primary.vs.meta.genes.pattern4,file=glue("{dir.path}/primary.vs.meta.pattern4.txt"),quote=F,row.names=F,sep="\t")

#data.primary.vs.meta <- inner_join(select(primary.vs.meta.untreated , -c(3:last_col())), select(primary.vs.meta.treated, -c(3:last_col())), by = c("genes"),keep=FALSE)

# clean up the colnames (.x and .y were added during the merge)
#colnames(data.primary.vs.meta) <- gsub('logFC.x', 'Prim.vs.Meta.Untreated', colnames(data.primary.vs.meta), fixed=TRUE)
#colnames(data.primary.vs.meta) <- gsub('logFC.y', 'Prim.vs.Meta.Treated', colnames(data.primary.vs.meta), fixed=TRUE)

#print(my.genes)
#filter(primary.vs.meta.untreated,genes %in%  my.genes)
#filter(primary.vs.meta.untreated,genes %in% "ESRP1" )

#write.table(filter(data.primary.vs.meta, genes %in%  my.genes),file=glue("{dir.path}/custom_matrice_heatmap_primary_vs_meta.txt"),quote=F,row.names=F,sep="\t")

#data.primary.vs.meta$pattern <- 0
#data.primary.vs.meta <- data.primary.vs.meta %>% select(pattern, everything())

#data.primary.vs.meta$pattern[data.primary.vs.meta$genes %in% primary.vs.meta.genes.pattern1] <- "1"
#data.primary.vs.meta$pattern[data.primary.vs.meta$genes %in% primary.vs.meta.genes.pattern2] <- "2"
#data.primary.vs.meta$pattern[data.primary.vs.meta$genes %in% primary.vs.meta.genes.pattern3] <- "3"
#data.primary.vs.meta$pattern[data.primary.vs.meta$genes %in% primary.vs.meta.genes.pattern4] <- "4"
#data.primary.vs.meta<- data.primary.vs.meta[data.primary.vs.meta$pattern!=0,]

#write.table(data.primary.vs.meta,file=glue("{dir.path}/matrice_heatmap_primary_vs_meta.txt"),quote=F,row.names=F,sep="\t")

# ====================================================================================
# Trash
# ====================================================================================

# remove duplicate columns
#data.h4.pc9 <- data.h4.pc9%>% select(-contains(".y"))
#clean up the colnames (.x and .y were added during the merge)
#colnames(data.h4.pc9) <- gsub('.x', '', colnames(data.h4.pc9), fixed=TRUE)

#cat(headerLine1,file="{dir.path}/matrice_heatmap.txt",sep="\t")
#cat("\n",file="{dir.path}/matrice_heatmap.txt",sep="\t",append=TRUE)
#cat(headerLine2,file="{dir.path}/matrice_heatmap.txt",sep="\t",append=TRUE)
#cat("\n",file="{dir.path}/matrice_heatmap.txt",sep="\t",append=TRUE)

# Renaming clearly
#colnames(data.h4.pc9.h3) <- gsub(colnames(data.h4.pc9.h3),pattern="^PC9_Erlo_DTC.*",replacement="PC9_Erlo_DTC",perl=T)
#colnames(data.h4.pc9.h3) <- gsub(colnames(data.h4.pc9.h3),pattern="^PC9-3_NT1.*",replacement="PC9-3_NT1",perl=T)
#colnames(data.h4.pc9.h3) <- gsub(colnames(data.h4.pc9.h3),pattern="^PC9-3_DTEC.*",replacement="PC9-3_DTEC",perl=T)
#colnames(data.h4.pc9.h3) <- gsub(colnames(data.h4.pc9.h3),pattern="^PC9_CT_24h.*",replacement="PC9_CT_24h",perl=T)
#colnames(data.h4.pc9.h3) <- gsub(colnames(data.h4.pc9.h3),pattern="^H4006_CT_24h.*",replacement="H4006_CT_24h",perl=T)
#colnames(data.h4.pc9.h3) <- gsub(colnames(data.h4.pc9.h3),pattern="^H4006_Erlo_DTEC.*",replacement="H4006_Erlo_DTEC",perl=T)
#colnames(data.h4.pc9.h3) <- gsub(colnames(data.h4.pc9.h3),pattern="^H4006_Erlo_DTC.*",replacement="H4006_Erlo_DTC",perl=T)
#colnames(data.h4.pc9.h3) <- gsub(colnames(data.h4.pc9.h3),pattern="^H3255_NT.*",replacement="H3255_NT",perl=T)
#colnames(data.h4.pc9.h3) <- gsub(colnames(data.h4.pc9.h3),pattern="^H3255_Erlo_21d.*",replacement="H3255_Erlo_21d",perl=T)

#headerLine1 <- gsub(colnames(data.h4.pc9.h3),pattern="^PC9_Erlo_DTC.*",replacement="T2",perl=T)
#headerLine1 <- gsub(headerLine1 ,pattern="^PC9-3_NT1.*",replacement="T1",perl=T)
#headerLine1 <- gsub(headerLine1 ,pattern="^PC9-3_DTEC.*",replacement="T3",perl=T)
#headerLine1 <- gsub(headerLine1 ,pattern="^PC9_CT_24h.*",replacement="T1",perl=T)
#headerLine1 <- gsub(headerLine1 ,pattern="^H4006_CT_24h.*",replacement="T1",perl=T)
#headerLine1 <- gsub(headerLine1 ,pattern="^H4006_Erlo_DTEC.*",replacement="T3",perl=T)
#headerLine1 <- gsub(headerLine1 ,pattern="^H4006_Erlo_DTC.*",replacement="T2",perl=T)
#headerLine1 <- gsub(headerLine1 ,pattern="^H3255_NT.*",replacement="T1",perl=T)
#headerLine1 <- gsub(headerLine1 ,pattern="^H3255_Erlo_21d.*",replacement="T2",perl=T)

#headerLine2 <- gsub(colnames(data.h4.pc9.h3),pattern="^PC9_Erlo_DTC.*",replacement="PC9",perl=T)
#headerLine2 <- gsub(headerLine2,pattern="^PC9-3_NT1.*",replacement="PC9",perl=T)
#headerLine2 <- gsub(headerLine2,pattern="^PC9-3_DTEC.*",replacement="PC9",perl=T)
#headerLine2 <- gsub(headerLine2,pattern="^PC9_CT_24h.*",replacement="PC9",perl=T)
#headerLine2 <- gsub(headerLine2,pattern="^H4006_CT_24h.*",replacement="H4006",perl=T)
#headerLine2 <- gsub(headerLine2,pattern="^H4006_Erlo_DTEC.*",replacement="H4006",perl=T)
#headerLine2 <- gsub(headerLine2,pattern="^H4006_Erlo_DTC.*",replacement="H4006",perl=T)
#headerLine2 <- gsub(headerLine2,pattern="^H3255_NT.*",replacement="H3255",perl=T)
#headerLine2 <- gsub(headerLine2,pattern="^H3255_Erlo_21d.*",replacement="H3255",perl=T)


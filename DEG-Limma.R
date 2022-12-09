library("edgeR")
library("RColorBrewer")
library("genefilter")
library("ggrepel")
library("ggfortify")
library("cluster")
library("factoextra")
library("ggplot2")
library("M3C")
library("AnnotationDbi")
library("xlsx")
library("clusterProfiler")
library("png")
library("curl")
library("corrplot")
library("tidyr")
library("dplyr")
library("WGCNA")

# FUNCTIONS
# Filter
filter<-function(filter="standard",data,min_group=3){
  if(filter == "standard"){
    keep<-rowSums(cpm(data)>1) >= min_group
    data<-data[keep,]
    data$samples$lib.size <-colSums(data$counts)
  }
  else if(filter == "bin"){
    keep<-rowSums(data$counts)>0
    data<-data[keep,]
    data$samples$lib.size <-colSums(data$counts)
  }
  else{
    stop("At the moment only bin/standard are supported")
  }
  return(data)
}

# DATA
# Reading counts
rawdata<-read.delim(file="/mnt/beegfs/sbonilla/Projects/Utero/Analisis/miARma/20221018_strandNo-Readcount_results/str-ReadCount.tab",header=TRUE,row.names=1,check.names = F)
colnames(rawdata)<-gsub("_nat","",colnames(rawdata))
colnames(rawdata)<-gsub("clean_reads_","",colnames(rawdata))
rawdata<-rawdata[,sort(colnames(rawdata))]

# Reading pheno
targets<-read.table(file = "/mnt/beegfs/sbonilla/Projects/Utero/Analisis/miARma/sampleinfo.txt",header = T,stringsAsFactors=F)
targets$Muestras_Utero<-substr(targets$Muestras_Utero, start = 1, stop = 9)
targets$Muestras_Utero<-gsub("_E","",targets$Muestras_Utero)

# Delete fail sample (Up_G008)
targets <- targets [-c(which(targets$Muestras_Utero == "Up_G008")),]
rownames(targets) = c(1:nrow(targets))

# Reorder rawdata and targets by group
targets <- targets[order(targets$Diagnosis_estandarizada), ]
rownames(targets) = c(1:nrow(targets))
idx = match(targets$Muestras_Utero,colnames(rawdata))
rawdata=rawdata[,idx]
colnames(rawdata)

# Check data: match and same order
all(colnames(rawdata) %in% targets$Muestras_Utero)
all(colnames(rawdata) == targets$Muestras_Utero)

# Group
group=targets$Diagnosis_estandarizada

# DELETE OUTLIERS SAMPLES
# Up_G027 and Up_G043
rawdata <- select(rawdata, -Up_G027, -Up_G043)
which(targets$Muestras_Utero == "Up_G043")
which(targets$Muestras_Utero == "Up_G027")
targets = targets[-c(which(targets$Muestras_Utero == "Up_G027"), which(targets$Muestras_Utero == "Up_G043")),]
rownames(targets) = c(1:nrow(targets))

# Check data: match and same order
all(colnames(rawdata) %in% targets$Muestras_Utero)
all(colnames(rawdata) == targets$Muestras_Utero)

# Group
group=targets$Diagnosis_estandarizada

# Reading gene length
gene.length<-read.table(file="/mnt/beegfs/sbonilla/Projects/Utero/Analisis/miARma/20221018_strandNo-Readcount_results/str-Size.tab",header=T)

# Reorder gene length like rawdata in result_counts
idx<-match(rownames(rawdata),gene.length$Gene)
results_counts<-gene.length[idx,]
results_counts[is.na(results_counts$Length),"Length"]<-0
nrow(results_counts)
# 60591
all(rownames(rawdata) == results_counts$Gene)

# ANALYSIS
# Design Matrix
design <- model.matrix(~0+group)

# DGE Norm and Filter by exp
dge<- DGEList(counts=rawdata, group=group, genes=results_counts) 

# How many samples in each group?
# table(targets$Diagnosis_estandarizada)
# Endometriosis   Factor_masculino                RIF Sin_causa_aparente        Utero_doble 
# 11                  8                 13                 10                  2 

# Filter DGEList with the corresponding min group:
#1 Factor_masculino  vs Endometriosis        --> min_group = 8
#2 Factor_masculino  vs RIF                  --> min_group = 8
#3 Factor_masculino  vs Sin_causa_aparente   --> min_group = 8
#4 Factor_masculino  vs Utero_doble          --> min_group = 2
#5 Endometriosis     vs RIF                  --> min_group = 11
#6 Endometriosis     vs Sin_causa_aparente   --> min_group = 10
#7 Endometriosis     vs Utero_doble          --> min_group = 2
#8 RIF               vs Sin_causa_aparente   --> min_group = 10
#9 RIF               vs Utero_doble          --> min_group = 2
#10 Sin_causa_aparente vs Utero_doble        --> min_group = 2


# Filter data by filter "standard": keep<-rowSums(cpm(data)>1) >= min_group

dge8 <-filter(filter="standard",dge,min_group = 8)
dge2 <-filter(filter="standard",dge,min_group = 2)
dge11 <-filter(filter="standard",dge,min_group = 11)
dge10 <-filter(filter="standard",dge,min_group = 10)

# n genes before filter: 
nrow(dge)
#60591
# n genes after filter:
nrow(dge8)
#25667
nrow(dge2)
#28313
nrow(dge11)
#25024
nrow(dge10)
#25217

# Library size normalization: Calculate scaling factors to convert raw library sizes into effective library sizes.
dgenorm8 <- calcNormFactors(dge8)
dgenorm2 <- calcNormFactors(dge2)
dgenorm11 <- calcNormFactors(dge11)
dgenorm10 <- calcNormFactors(dge10)

# RPKMs: Compute reads per kilobase per million
rpkm8 <-rpkm(dgenorm8,normalized.lib.sizes=TRUE)
rpkm2 <-rpkm(dgenorm2,normalized.lib.sizes=TRUE)
rpkm11<-rpkm(dgenorm11,normalized.lib.sizes=TRUE)
rpkm10<-rpkm(dgenorm10,normalized.lib.sizes=TRUE)

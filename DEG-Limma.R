# Reading counts
rawdata<-read.delim(file="/mnt/beegfs/sbonilla/Projects/Utero/Analisis/miARma/20221018_strandNo-Readcount_results/str-ReadCount.tab",header=TRUE,row.names=1,check.names = F)
colnames(rawdata)<-gsub("_nat","",colnames(rawdata))
colnames(rawdata)<-gsub("clean_reads_","",colnames(rawdata))
rawdata<-rawdata[,sort(colnames(rawdata))]

# Reading pheno
targets<-read.table(file = "/mnt/beegfs/sbonilla/Projects/Utero/Analisis/miARma/sampleinfo.txt",header = T,stringsAsFactors=F)
targets$Muestras_Utero<-substr(targets$Muestras_Utero, start = 1, stop = 9)
targets$Muestras_Utero<-gsub("_E","",targets$Muestras_Utero)

# Reorder targets by Sample Name
# targets <- targets[order(targets$Muestras_Utero), ]
# rownames(targets) = c(1:nrow(targets))

# Delete fail sample (Up_G008)
targets <- targets [-c(which(targets$Muestras_Utero == "Up_G008")),]
rownames(targets) = c(1:nrow(targets))

# Reorder rawdata and targets by group
targets <- targets[order(targets$Diagnosis_estandarizada), ]
rownames(targets) = c(1:nrow(targets))
idx = match(targets$Muestras_Utero,colnames(rawdata))
rawdata=rawdata[,idx]
colnames(rawdata)

#Check data: match and same order
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

#reorder gene length like rawdata in result_counts
idx<-match(rownames(rawdata),gene.length$Gene)
results_counts<-gene.length[idx,]
results_counts[is.na(results_counts$Length),"Length"]<-0
nrow(results_counts)
# 60591
all(rownames(rawdata) == results_counts$Gene)


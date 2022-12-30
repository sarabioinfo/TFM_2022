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
Sara_max <-function(data,n=30){
  data_df <- as.data.frame(data)
  data_df_n <- sample_n(data_df, n, replace = TRUE)
  index_data_df_n <- as.data.frame(apply(data_df_n, 1, which.max))
  max = vector()
  for (i in 1:n){
    colnames(data_df_n[index_data_df_n[i,]])
    max = sort(c(max,colnames(data_df_n[index_data_df_n[i,]]) ))
  }
  data_sort <- sort(table(max), decreasing = TRUE)
  cat(paste("How many times is the most expressed in",n,"random elements: "))
  return (data_sort)
  
}
Sara_QC<-function(label,dge,data){
  x<-dge
  group<-x$samples$group
  log_Edu<-log(rawdata,2)
  cpm <- cpm(x)
  lcpm <- cpm(x, log=TRUE)
  L <- mean(x$samples$lib.size) * 1e-6
  M <- median(x$samples$lib.size) * 1e-6
  c(L, M)
  
  summary(lcpm)
  table(rowSums(x$counts==0)==6)
  
  pdf(paste(path,label,"_QC_test.pdf",sep=""),paper="A4")
  lcpm.cutoff <- log2(10/M + 2/L)
  
  nsamples <- ncol(x)
  col <- rainbow(ncol(x))
  
  plot(density(log_Edu[,1]), col=col[1], lwd=2, las=2, main="", xlab="")
  title(main="A. Raw data", xlab="Log-cpm")
  abline(v=lcpm.cutoff, lty=3)
  for (i in 2:nsamples){
    den <- density(log_Edu[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright",legend=targets$Muestras_Utero, text.col=col, bty = "n", cex = 0.4)
  
  
  lcpm <- cpm(x, log=TRUE)
  
  plot(density(lcpm[,1]), col=col[1], lwd=2, las=2, main="", xlab="")
  title(main="B. Filtered data", xlab="Log-cpm")
  abline(v=lcpm.cutoff, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", legend=targets$Muestras_Utero, text.col=col, bty="n", cex = 0.4)
  
  #normalizing
  x <- calcNormFactors(x)
  x$samples
  x <- estimateCommonDisp(x, robust=TRUE)
  x <- estimateTagwiseDisp(x)
  
  col.group <- as.factor(group)
  if (nlevels(col.group) >=3 ) {
    levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
  }
  else 
  {levels(col.group) <- c("red", "darkgreen")}
  col.group <- as.character(col.group)
  
  lcpm <- cpm(dge, log=TRUE)
  boxplot(lcpm, las=2, col=col.group, main="", names=targets$Muestras_Utero, cex.axis=0.4)
  title(main="Unnormalised data",ylab="Log-cpm")
  
  lcpm <- cpm(x, log=TRUE)
  boxplot(lcpm, las=2, col=col.group, main="", names=targets$Muestras_Utero, cex.axis=0.4)
  title(main="Normalised data",ylab="Log-cpm")
  
  #Libray size
  barplot(x$samples$lib.size,names.arg = targets$Muestras_Utero,las=2, main="Library Size",col=col.group, ylim=range(pretty(c(0, x$samples$lib.size))), cex.names = 0.4)
  
  lcpm <- cpm(x, log= FALSE) # cpm de los datos normalizados (sin log)
  
  corrplot(cor(lcpm,method="spearman"), method='number',type = 'upper', tl.cex = 0.4, number.cex = 0.3)
  corrplot(cor(lcpm,method="spearman"), method='color',type = 'upper', tl.cex = 0.4) # method=color porque si no, no se ve bien
  corrplot(cor(lcpm,method="spearman"), order='AOE',tl.cex = 0.4)
  
  col.group <- as.factor(group)
  levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
  col.group <- as.character(col.group)
  levels(col.group)
  lcpm <- cpm(x, log=TRUE)
  z<-plotMDS(lcpm, labels=targets$Muestras_Utero, col=col.group, gene.selection = "pairwise", plot=F)
  edge<-sd(z$x)
  #par(mar=c(5, 5, 4, 11), xpd=TRUE)
  #par(xpd=TRUE)
  plotMDS(lcpm, labels=targets$Muestras_Utero, col=col.group, gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))

  title(main="A. MDS-PCoA Sample Names")
  
  # with legend
  plotMDS(lcpm, labels=targets$Muestras_Utero, col=col.group, gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
  title(main="A. MDS-PCoA Sample Names")
  # legend:
  group <- as.factor(group)
  col.group <- as.factor(group)
  col.group <-  brewer.pal(nlevels(group), "Set1")
  col.group
  legend("topleft", legend =  levels(group), text.col = col.group)
  
  
  #PCA por tipos
  # as.matrix(x) = matriz de conteos
  data_pca<-as.matrix(x)
  # t transpone los datos: las muestras son las filas y los genes las columnas
  data_pca<-as.data.frame(t(data_pca))
  
  #rownames(data_pca)<-targets$Muestras_Utero
  
  #prcomp  = Performs a principal components analysis on the given data matrix and returns the results as an object of class prcomp.
  # prcomp transforma los resultados para que tengan de media 0. Si queremos que la desv estandard sea 1, añadir scale = TRUE
  
  data_pca.PC = prcomp(data_pca)
  data_pca$Name<-targets$Muestras_Utero
  data_pca$Age<-targets$Edad
  data_pca$BMI<-targets$BMI
  data_pca$Diagnosis_estandarizada<-targets$Diagnosis_estandarizada
  
  
  #plot(autoplot(data_pca.PC,label=T,data=data_pca,colour='Name',xlim = c(-10,10))) + theme(legend.position='none')
  plot(autoplot(data_pca.PC,label=T,data=data_pca,colour='Name',xlim = c(-0.25,0.8)))  + theme(legend.position='none')
  plot(autoplot(data_pca.PC,label=T,data=data_pca,colour='Age',xlim = c(-0.25,0.8)))
  plot(autoplot(data_pca.PC,label=T,data=data_pca,colour='BMI',xlim = c(-0.25,0.8)))
  plot(autoplot(data_pca.PC,label=T,data=data_pca,colour='Diagnosis_estandarizada',xlim = c(-0.25,0.8)))
  
  # Con desv estandard 1 scale=true para que no destaquen genes por encima de otros:
  data_pca<-as.matrix(x)
  data_pca<-as.data.frame(t(data_pca))
  data_pca.PC = prcomp(data_pca, scale = TRUE)
  data_pca$Name<-targets$Muestras_Utero
  data_pca$Age<-targets$Edad
  data_pca$BMI<-targets$BMI
  data_pca$Diagnosis_estandarizada<-targets$Diagnosis_estandarizada
  #plot(autoplot(data_pca.PC,label=T,data=data_pca,colour='Name',xlim = c(-10,10)))
  plot(autoplot(data_pca.PC,label=T,data=data_pca,colour='Name',xlim = c(-0.25,0.8)))   + theme(legend.position='none')
  plot(autoplot(data_pca.PC,label=T,data=data_pca,colour='Age',xlim = c(-0.25,0.8)))
  plot(autoplot(data_pca.PC,label=T,data=data_pca,colour='BMI',xlim = c(-0.25,0.8)))
  plot(autoplot(data_pca.PC,label=T,data=data_pca,colour='Diagnosis_estandarizada',xlim = c(-0.25,0.8)))
  
  #### HEATMAP ####
  # desv estandard por fila de los valores de conteo normalizados de todas las muestras:
  rsd <- rowSds(as.matrix(x)) 
  
  # Ordenar la desviación estandard de mayor a menor, y seleccionar los 250 con más desviación
  sel <- order(rsd, decreasing=TRUE)[1:250]
  samplenames<-as.character(targets$Muestras_Utero)
  heatmap(na.omit(as.matrix(x[sel,])),margins=c(10,8),main="Heatmap 250 most diff entities",cexRow=0.01,cexCol=0.5,labCol=samplenames)
  
  #cluster
  par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
  #pr.hc.c<- hclust(na.omit(dist(t(data))))
  #plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of ", label, sep=""), labels=targets$Filemane, cex=0.5)
  #Normalized clustering analysis plot
  #pr.hc.c<- hclust(na.omit(dist(t(dgenorm$counts))))
  pr.hc.c<- hclust(na.omit(dist(t(cpm(x$counts,log=T)),method = "euclidean")))
  plot(pr.hc.c, xlab="Sample Distance",main="Hierarchical Clustering of Normalized samples", labels=targets$Muestras_Utero, cex=0.5)
  
  #tSNE
  #a<-tsne(x$counts,seed=100,labels=as.factor(targets$Type), perplex=perplex, legendtitle="Types",text=targets$Type ,dotsize=3, legendtextsize = 8) + ggtitle("Tsne") + theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
  #plot(a)  
  dev.off()
  
  # jpeg(paste0(path,"PCoA_QC.jpg"))
  # plotMDS(lcpm, labels=targets$Type, col=col.group, gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
  # dev.off()
  # 
  # jpeg(paste0(path,"Samples_HTree_QC.jpg"))
  # plot(pr.hc.c, xlab="Sample Distance", labels=targets$Name, cex=0.5)
  # dev.off()
  
}
Sara_DGE_limma <- function(top, pvalue = 0.05 , logFC = 1){
  cat(paste("Number of genes before filter: ", nrow(top),sep=""),"\n")
  data <- top[top$P.Value <= pvalue & abs(top$logFC) >= logFC,]
  cat(paste("Number of genes after filter: ", nrow(data),sep=""),"\n")
  return (data)
  
}
Sara_enrichment<-function(de_data,FA_label,pvalue=0.05,cutoff=0.05,organism="human",logFC){
  dbName<-NULL
  keggDB<-NULL
  if(organism == "human"){
    library("org.Hs.eg.db")
    dbName<-org.Hs.eg.db
    keggDB<-"hsa"
  }else if (organism == "mouse"){
    library("org.Mm.eg.db")
    dbName<-org.Mm.eg.db
    keggDB<-"mmu"
  }else{
    stop("At the moment only human and mouse are supported")
  }
  
  de_data_filter <- Sara_DGE_limma(de_data, pvalue = pvalue, logFC = logFC)
  gene <- de_data_filter$Gene
  
  eg = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=dbName)
  all_eg = bitr(de_data$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=dbName)
  
  new_path=paste(path,"Enrichment/",sep="")
  system(paste("mkdir -p ", new_path,sep=""))
  ego_BP <- enrichGO(gene          = as.vector(eg$ENTREZID),
                     universe      = as.vector(all_eg$ENTREZID),
                     OrgDb         = dbName,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = cutoff,
                     qvalueCutoff  = cutoff,
                     readable      = TRUE)
  ego_MF <- enrichGO(gene          = as.vector(eg$ENTREZID),
                     universe      = as.vector(all_eg$ENTREZID),
                     OrgDb         = dbName,
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = cutoff,
                     qvalueCutoff  = cutoff,
                     readable      = TRUE)
  
  ego_CC <- enrichGO(gene          = as.vector(eg$ENTREZID),
                     universe      = as.vector(all_eg$ENTREZID),
                     OrgDb         = dbName,
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = cutoff,
                     qvalueCutoff  = cutoff,
                     readable      = TRUE)
  
  if(nrow(as.data.frame(ego_BP))>0){
    write.table(format_enrichment(ego_BP,de_data),file = paste(new_path,"Enrichment_GO_BP_",FA_label,".xls",sep=""),sep="\t",row.names = F)
    pdf(file = paste(new_path,"Enrichment_GO_BP_",FA_label,".pdf",sep=""),paper="a4")
    y<-barplot(ego_BP)
    y$labels$title<-"GO terms from BP\n"
    plot(y)
    
    x<-dotplot(ego_BP)
    x$labels$title<-"GO terms from BP\n"
    plot(x)
    
    #nose<-ego_BP
    #class(nose)<-NULL
    #z<-cnetplot(ego_BP,categorySize = "pvalue",main="GO terms from BP\n",vertex.label.cex=0.5)
    #plot(z)
    
    plotGOgraph(ego_BP, firstSigNodes = 5,  sigForAll = TRUE, useFullNames = TRUE)
    dev.off()
    
    # png(file = paste(new_path,"GO_BP_dotplot",FA_label,".png",sep=""))
    # plot(x)
    # dev.off()
  }
  if(nrow(as.data.frame(ego_MF))>0){
    
    write.table(format_enrichment(ego_MF,de_data),file = paste(new_path,"Enrichment_GO_MF_",FA_label,".xls",sep=""),sep="\t",row.names = F)
    pdf(file = paste(new_path,"Enrichment_GO_MF_",FA_label,".pdf",sep=""),paper="a4")
    y<-barplot(ego_MF)
    y$labels$title<-"GO terms from MF\n"
    plot(y)
    
    x<-dotplot(ego_MF)
    x$labels$title<-"GO terms from MF\n"
    plot(x)
    
    #z<-cnetplot(ego_MF,categorySize = "pvalue",main="GO terms from MF\n",vertex.label.cex=0.5)
    #plot(z)
    plotGOgraph(ego_MF, firstSigNodes = 5,  sigForAll = TRUE, useFullNames = TRUE)
    dev.off()
    
    # png(file = paste(new_path,"GO_MF_dotplot",FA_label,".png",sep=""))
    # plot(x)
    # dev.off()
  }
  if(nrow(as.data.frame(ego_CC))>0){
    write.table(format_enrichment(ego_CC,de_data),file = paste(new_path,"Enrichment_GO_CC_",FA_label,".xls",sep=""),sep="\t",row.names = F)
    pdf(file = paste(new_path,"Enrichment_GO_CC_",FA_label,".pdf",sep=""),paper="a4")
    y<-barplot(ego_CC)
    y$labels$title<-"GO terms from CC\n"
    plot(y)
    
    x<-dotplot(ego_CC)
    x$labels$title<-"GO terms from CC\n"
    plot(x)
    
    #z<-cnetplot(ego_CC,categorySize = "pvalue",main="GO terms from CC\n",vertex.label.cex=0.3,node.color="red",edge.width=0.2)
    #plot(z)
    plotGOgraph(ego_CC, firstSigNodes = 5,  sigForAll = TRUE, useFullNames = TRUE)
    dev.off()
    
    # png(file = paste(new_path,"GO_CC_dotplot",FA_label,".png",sep=""))
    # plot(x)
    # dev.off()
  }
  
  KEGG_enrich <- enrichKEGG(
    gene= as.vector(eg$ENTREZID),
    keyType = "ncbi-geneid",
    organism     = keggDB,
    pvalueCutoff = cutoff,
    use_internal_data = F
  )
  if(is.null(KEGG_enrich) == F){
    if(nrow(KEGG_enrich)>0){
      write.table(KEGG_enrich,file = paste(new_path,"Enrichment_KEGG_",FA_label,".xls",sep=""),sep="\t",row.names = F)
      
      result<-NA
      for(cont in 1:length(KEGG_enrich$ID)){
        cat(paste("For pathway named: ", KEGG_enrich$ID[cont]))
        result[cont]<-browseKEGG(KEGG_enrich, KEGG_enrich$ID[cont])
        #print(result[cont])
        Sys.sleep(2)
        file<-readLines(result[cont])
        #cat(file)
        pattern='<.+ src=\"(.+)\" srcset='
        #cat(pattern)
        path<-gsub(pattern,"\\1", grep("pathwayimage",file,value=TRUE))
        #cat(path)
        result[cont]<-strsplit(path," ")[[1]][1]
        #cat(result[cont])
        result[cont]<-gsub("\\\"","",result[cont])
        result[cont]<-gsub("\\t","",result[cont])
        result[cont]<-gsub("\\.png.*",".png",result[cont])
        #cat(result[cont])
        my_url<-paste("http://www.kegg.jp",result[cont],sep="")
        #cat(my_url)
        print(paste("Downloading ", my_url))
        curl_download(my_url,paste(KEGG_enrich$ID[cont],".png",sep=""),mode="wb")
      }
      
      files <- list.files(path=".", pattern="*.png", all.files=T, full.names=T)
      
      pdf(file = paste(new_path,"Enrichment_KEGG_",FA_label,".pdf",sep=""),paper = "a4")
      for(i in 1:length(files)) {
        img <- readPNG(files[i])
        plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n",axes=F,xlab="", ylab="")
        rasterImage(img,0,0,1,1)
      }
      dev.off()
      a<-apply(as.data.frame(files),1,function(x) unlink(x))
    }  
  } 
}
Sara_Volcano_limma <- function (data,  title = "GroupA_GroupB", logFC = 1, selected = 20){
  data$Expression<-"Unchanged"
  data[data$P.Value <= 0.05 & data$logFC >= logFC ,"Expression"]<-"Up-regulated"
  data[data$P.Value <= 0.05 & data$logFC <= -logFC ,"Expression"]<-"Down-regulated"  
  
  real_DE<-data[data$P.Value<=0.05,]
  selected_FC<-head(real_DE[order(abs(real_DE$logFC),decreasing = T),], selected)
  max_value = max(abs(data$logFC))
  g = ggplot(data=data, aes(x=logFC, y=-log10(P.Value), colour=Expression)) + coord_cartesian(xlim = c(-max_value, max_value )) +
    geom_point(alpha=0.4, size=1.75) +
    xlab("log2 Fold Change") + ylab("-log10 P Value") +
    geom_text_repel(data=selected_FC, aes(label=Gene),colour="black",size=3) + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept = c(-logFC, logFC),lty=4, size = 0.1) +
    geom_hline(yintercept = -log10(0.05),lty=4, size = 0.1) +
    scale_color_manual(values=c("#56B4E9","#999999","#E69F00"))
  
  # jpeg(paste0(path,"Volcano_plot_",title,".jpg"))
  # plot(g)
  # dev.off()
  pdf(paste(path,"Volcano_plot",title,".pdf",sep=""),paper="a4")
  plot(g)
  dev.off()
}
Sara_myDGE<-function(filter,min_group=3){
  #creating DGE object
  dge<-DGEList(counts=rawdata, group=group, genes=results_counts) 
  dge$sample
  #                        group lib.size norm.factors
  # Up_G022        Endometriosis 37625120            1
  # Up_G023        Endometriosis 42215195            1
  # Up_G040        Endometriosis 39820221            1
  # Up_G041        Endometriosis 47078941            1
  # Up_G044        Endometriosis 38507072            1
  
  #filtering low expressed genes
  table(dge$samples$group)
  # Endometriosis   Factor_masculino        RIF   Sin_causa_aparente        Utero_doble 
  #            11                  8        13                    10                  2
  
  cat(paste("Number of genes before filter [ ",filter," ]: ", nrow(dge),sep=""),"\n")
  dge<-filter(filter=filter,dge,min_group = min_group)
  # Number of genes before filter [ standard ]: 60591 
  
  cat(paste("Number of genes after filter [ ",filter," ]: ", nrow(dge),sep=""),"\n")
  nrow(dge)
  #25667
  
  #normalizing
  dgenorm <- calcNormFactors(dge)
  # 3 elements: $counts, $samples (sample, lib.size and norm.factors) and $genes
  dgenorm <- estimateCommonDisp(dgenorm, robust=TRUE)
  # 7 elements: previous 3 + $common.dispersion (0.1389135), $pseudo.counts, $pseudo.lib.size (40182342) and $AveLogCPM
  
  dgenorm <- estimateTagwiseDisp(dgenorm)
  # 11 elements: previous 7 + $prior.df (10), $prior.n (0.2564103), $tagwise.dispersion and $span (0.3)
  
  #Quality
  Sara_QC(label=label,dge=dge,data=rawdata)
  
  return(dgenorm)
}
filter<-function(filter="standard",data,min_group=3){
  if(filter == "standard"){
    # Mantiene los genes donde haya al menos min_groups muestras con cpm > 1
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
project<-function(){
  path=paste(getwd(),"/",label,"/",sep="")
  system(paste("mkdir -p ", path,sep=""))
  cat(paste("All files will be saved under folder:", path),"\n")
  return(path)
}
createExcel<-function(data,file,organism="human"){
  dbName<-NULL
  if(organism == "human"){
    library("org.Hs.eg.db")
    dbName<-org.Hs.eg.db
  }else if (organism == "mouse"){
    library("org.Mm.eg.db")
    dbName<-org.Mm.eg.db
    
  }else{
    stop("At the moment only human and mouse are supported")
  }
  options(java.parameters = "-Xmx2048m")
  
  wb<-createWorkbook(type="xlsx")
  
  CellStyle(wb, dataFormat=NULL, alignment=NULL, border=NULL, fill=NULL, font=NULL)
  
  TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) +
    Alignment(wrapText=FALSE, horizontal="ALIGN_CENTER") +
    Border(color="black", position=c("TOP", "BOTTOM", "LEFT", "RIGHT"), 
           pen=c("BORDER_THIN", "BORDER_THIN", "BORDER_THIN", "BORDER_THIN"))
  
  TextFormat = CellStyle(wb, dataFormat=DataFormat("@"))
  FormatList = list('1'=TextFormat, '2'=TextFormat,'3' = TextFormat,'4' = TextFormat)
  
  sheet1 <- xlsx::createSheet(wb, sheetName='DE Genes')
  
  x<-as.data.frame(mapIds(
    dbName, 
    keys = data$Gene, 
    "GENENAME", 
    "SYMBOL",
    fuzzy = TRUE,
    multiVals = "first"))
  
  data$Description<-x[,1]
  
  data<-data[,c(1,7,2:6)]
  xlsx::addDataFrame(data, sheet1, col.names=TRUE, colStyle=FormatList, 
                     row.names=FALSE, colnamesStyle = TABLE_COLNAMES_STYLE)
  
  addAutoFilter(sheet1, paste(LETTERS[1],LETTERS[ncol(data)],sep="-"))
  for (column in 1:ncol(data)) {
    #print(paste("LA columna ", column, " tiene un total de ",  max(nchar(na.omit(data[,column])))))
    size=max(nchar(na.omit(data[,column])))
    if(size <10) size<-15
    setColumnWidth(sheet1, column, size)
  }
  
  saveWorkbook(wb, file = paste(path,file,sep=""))
  write.table(data, paste0(path, gsub("DEG_","",gsub(".xlsx",".tsv",file))), sep = "\t", row.names = F)
}
format_enrichment<-function(ego,de_data){
  #de_data=Endometriosis_FactorMasculino
  data<-data.frame(ego)
  colnames(data)<-gsub("Count","Total",colnames(data))
  genes_down<-de_data[de_data$P.Value<=0.05 & de_data$logFC< (-0.58),"Gene"] 
  genes_up<-de_data[de_data$P.Value<=0.05 & de_data$logFC> 0.58,"Gene"] 
  data$GenesID_up<-NA
  data$GenesID_dn<-NA
  data$Total_up<-NA
  data$Total_dn<-NA
  
  for(n in 1:nrow(data)){
    up<-genes_up[genes_up %in% unlist(strsplit(data[n,"geneID"],"/"))]
    dn<-genes_down[genes_down %in% unlist(strsplit(data[n,"geneID"],"/"))]
    data[n,"GenesID_up"]<-paste(up,collapse = "/")
    data[n,"GenesID_dn"]<-paste(dn,collapse = "/")
    data[n,"Total_up"]<-length(up)
    data[n,"Total_dn"]<-length(dn)
  }
  return(data[,c(1:8,10,11,9,12,13)])
}

#--------------------------------------------------------------------------------------------------#
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
targets = targets[-c(which(targets$Muestras_Utero == "Up_G027"), which(targets$Muestras_Utero == "Up_G043")),]
rownames(targets) = c(1:nrow(targets))

# DELETE DOUBLE UTERUS SAMPLES
# Up_G051_D and Up_G051_I
rawdata <- select(rawdata, -Up_G051_D, -Up_G051_I)
targets = targets[-c(which(targets$Muestras_Utero == "Up_G051_D"), which(targets$Muestras_Utero == "Up_G051_I")),]
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
# Path
label<-("Utero")
path<-project()

# Design Matrix
design <- model.matrix(~0+group)

# DGE Norm and Filter by exp
dge<- DGEList(counts=rawdata, group=group, genes=results_counts) 

# How many samples in each group?
# table(targets$Diagnosis_estandarizada)
# Endometriosis   Factor_masculino                RIF Sin_causa_aparente 
# 11                  8                 13                 10        

# Filter DGEList with the corresponding min group:
#1 Factor_masculino  vs Endometriosis        --> min_group = 8
#2 Factor_masculino  vs RIF                  --> min_group = 8
#3 Factor_masculino  vs Sin_causa_aparente   --> min_group = 8
#4 Endometriosis     vs RIF                  --> min_group = 11
#5 Endometriosis     vs Sin_causa_aparente   --> min_group = 10
#6 RIF               vs Sin_causa_aparente   --> min_group = 10


# Filter data by filter "standard": keep<-rowSums(cpm(data)>1) >= min_group

dge8 <-filter(filter="standard",dge,min_group = 8)
dge11 <-filter(filter="standard",dge,min_group = 11)
dge10 <-filter(filter="standard",dge,min_group = 10)

# n genes before filter: 
nrow(dge)
#60591
# n genes after filter:
nrow(dge8)
#25574
nrow(dge11)
#24943
nrow(dge10)
#25128

# Library size normalization: Calculate scaling factors to convert raw library sizes into effective library sizes.
dgenorm8 <- calcNormFactors(dge8)
dgenorm11 <- calcNormFactors(dge11)
dgenorm10 <- calcNormFactors(dge10)

# RPKMs: Compute reads per kilobase per million
rpkm8 <-rpkm(dgenorm8,normalized.lib.sizes=TRUE)
rpkm11<-rpkm(dgenorm11,normalized.lib.sizes=TRUE)
rpkm10<-rpkm(dgenorm10,normalized.lib.sizes=TRUE)

# Transform count data to log2CPM estimate the mean-variance relationship and use this to compute appropriate observation-level weights. 
# The data are then ready for linear modelling (normalize="quantile")
v8 <- voom(dgenorm8, design, plot = TRUE)
v11 <- voom(dgenorm11, design, plot = TRUE)
v10 <- voom(dgenorm10, design, plot = TRUE)

# lmFit (mean counts matrix per condition). Fit linear model for each gene given a series of arrays
fit8 <- lmFit(v8, design)
fit11 <- lmFit(v11, design)
fit10 <- lmFit(v10, design)

# contrast = all groups vs all groups
contrast <- makeContrasts(
  groupEndometriosis      - groupFactor_masculino,
  groupRIF                - groupFactor_masculino,
  groupSin_causa_aparente - groupFactor_masculino,
  groupRIF                - groupEndometriosis,
  groupSin_causa_aparente - groupEndometriosis,
  groupSin_causa_aparente - groupRIF,
  levels=design)

# contrasts.fit: Compute Contrasts from Linear Model Fit
# Given a linear model fit to microarray data, compute estimated coefficients and standard errors for a given set of contrasts.
contrast_fit8 <- contrasts.fit(fit8, contrast)
contrast_fit11 <- contrasts.fit(fit11, contrast)
contrast_fit10 <- contrasts.fit(fit10, contrast)

# eBayes: Empirical Bayes Statistics for Differential Expression
# Given a linear model fit from lmFit, compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression by empirical Bayes moderation of the standard errors towards a global value.

contrast_result8 <- eBayes(contrast_fit8)
contrast_result11 <- eBayes(contrast_fit11)
contrast_result10 <- eBayes(contrast_fit10)

# Tables by comparisons
Endometriosis_FactorMasculino    <- topTable(contrast_result8, coef = 1, sort.by = "P", number = Inf)
RIF_FactorMasculino              <- topTable(contrast_result8, coef = 2, sort.by = "P", number = Inf)
SinCausa_FactorMasculino         <- topTable(contrast_result8, coef = 3, sort.by = "P", number = Inf)
RIF_Endometriosis                <- topTable(contrast_result11, coef = 4, sort.by = "P", number = Inf)
SinCausa_Endometriosis           <- topTable(contrast_result10, coef = 5, sort.by = "P", number = Inf)
SinCausa_RIF                     <- topTable(contrast_result10, coef = 6, sort.by = "P", number = Inf)

# Create excels
createExcel(Endometriosis_FactorMasculino,"Endometriosis_FactorMasculino.xlsx",organism = "human")
createExcel(RIF_FactorMasculino,"RIF_FactorMasculino.xlsx",organism = "human")
createExcel(SinCausa_FactorMasculino,"SinCausa_FactorMasculino.xlsx",organism = "human")
createExcel(RIF_Endometriosis,"RIF_Endometriosis.xlsx",organism = "human")
createExcel(SinCausa_Endometriosis,"SinCausa_Endometriosis.xlsx",organism = "human")
createExcel(SinCausa_RIF,"SinCausa_RIF.xlsx",organism = "human")

# Volcano plots
Sara_Volcano_limma (Endometriosis_FactorMasculino,  title = "Endometriosis vs Control", logFC = 0.58, selected = 10)
Sara_Volcano_limma (RIF_FactorMasculino,  title = "RIF vs Control", logFC = 0.58, selected = 10)
Sara_Volcano_limma (SinCausa_FactorMasculino,  title = "Unexplained Infertility  vs Control", logFC = 0.58, selected = 10)
Sara_Volcano_limma (RIF_Endometriosis,  title = "RIF vs Endometriosis", logFC = 0.58, selected = 10)
Sara_Volcano_limma (SinCausa_Endometriosis,  title = "Unexplained Infertility vs Endometriosis", logFC = 0.58, selected = 10)
Sara_Volcano_limma (SinCausa_RIF,  title = "Unexplained Infertility vs RIF", logFC = 0.58, selected = 10)

# Filter by pvalue =< 0.05 and logFC >= |0.58|
dge_Endometriosis_FactorMasculino_logFC0.58 <- Sara_DGE_limma(Endometriosis_FactorMasculino, pvalue = 0.05, logFC = 0.58)
# Number of genes before filter: 25574  
# Number of genes after filter (logFC = 0.58): 743 
dge_RIF_FactorMasculino_logFC0.58 <- Sara_DGE_limma(RIF_FactorMasculino, pvalue = 0.05, logFC = 0.58)
# Number of genes before filter:  25574 
# Number of genes after filter (logFC = 0.58):  587 
dge_SinCausa_FactorMasculino_logFC0.58 <- Sara_DGE_limma(SinCausa_FactorMasculino, pvalue = 0.05, logFC = 0.58)
# Number of genes before filter:  25574 
# Number of genes after filter (logFC = 0.58):  140 
dge_RIF_Endometriosis_logFC0.58 <- Sara_DGE_limma(RIF_Endometriosis, pvalue = 0.05, logFC = 0.58)
# Number of genes before filter:  24943 
# Number of genes after filter (logFC = 0.58):  245 
dge_SinCausa_Endometriosis_logFC0.58 <- Sara_DGE_limma(SinCausa_Endometriosis, pvalue = 0.05, logFC = 0.58)
# Number of genes before filter:  25128 
# Number of genes after filter (logFC = 0.58):  709 
dge_SinCausa_RIF_logFC0.58 <- Sara_DGE_limma(SinCausa_RIF, pvalue = 0.05, logFC = 0.58) 
# Number of genes before filter:  25128 
# Number of genes after filter (logFC = 0.58):  288 

# Enrichment 
# Endometriosis_FactorMasculino
Sara_enrichment(Endometriosis_FactorMasculino,"Endometriosis_FactorMasculino_logFC0.58",pvalue=0.05,organism="human", logFC = 0.58)
# RIF_FactorMasculino
Sara_enrichment(RIF_FactorMasculino,"RIF_FactorMasculino_logFC0.58",pvalue=0.05,organism="human", logFC = 0.58)
# SinCausa_FactorMasculino
Sara_enrichment(SinCausa_FactorMasculino,"SinCausa_FactorMasculino_logFC0.58",pvalue=0.05,organism="human", logFC = 0.58)
Sara_enrichment(SinCausa_FactorMasculino,"SinCausa_FactorMasculino_logFC0.58",pvalue=0.05,organism="human", logFC = 1)
# RIF_Endometriosis
Sara_enrichment(RIF_Endometriosis,"RIF_Endometriosis_logFC0.58",pvalue=0.05,organism="human", logFC = 0.58)
# SinCausa_Endometriosis
Sara_enrichment(SinCausa_Endometriosis,"SinCausa_Endometriosis_logFC0.58",pvalue=0.05,organism="human", logFC = 0.58)
# SinCausa_RIF
Sara_enrichment(SinCausa_RIF,"SinCausa_RIF_logFC0.58",pvalue=0.05,organism="human", logFC = 0.58)
               
# Venn diagram
library(ggven)
x <- list("Endometriosis vs Ctrl" = dge_Endometriosis_FactorMasculino_logFC0.58$Gene, 
          "RIF vs Ctrl" = dge_RIF_FactorMasculino_logFC0.58$Gene, 
          "Unexplained Infertility vs Ctrl" = dge_SinCausa_FactorMasculino_logFC0.58$Gene)
ggvenn(x, 
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
       stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE)


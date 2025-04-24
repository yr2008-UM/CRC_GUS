# ######################################################################
#
## Correlation between gmGUS and Metabolites/KOs (Figure 5)
## Author: JR C
## Date: 01.01.2025
#
# ######################################################################
library(dplyr)
library(ggplot2)

##### Functions #####
enrichMBandKO <- function(edgeOUTV2,mark){
  if(mark == 'KO'){
    ##### enrichment analysis #####
    inputKOs <- edgeOUTV2$from %>% as.vector() %>% unique
    mSet <- enricher(inputKOs,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 1,pAdjustMethod = "BH",qvalueCutoff = 1)
    mSetEnrich <- NULL
    mSetClass <- NULL
    mSetMap <- NULL
  }else{
    if(mark == 'MB'){
      ##### enrichment analysis #####
      cmpd.vec <- edgeOUTV2$from %>% unique()
      cmpd.vec <- gsub("_.*","",cmpd.vec)
      cmpd.vec <- gsub("\"","",cmpd.vec)
      mSet <- InitDataObjects("conc", "msetora", FALSE)
      mSet<-Setup.MapData(mSet, cmpd.vec);
      mSet<-CrossReferencing(mSet, "kegg");
      mSet<-CreateMappingResultTable(mSet)
      mSet<-SetMetabolomeFilter(mSet, F);
      mSet<-SetCurrentMsetLib(mSet, "main_class", 2);
      mSet<-CalculateHyperScore(mSet)
      mSetEnrich <- mSet$analSet$ora.mat %>% as.data.frame() %>% dplyr::filter(FDR < 0.05) %>% rownames()
      mSetClass <- mSet$analSet$ora.hits %>% unlist %>% as.data.frame() %>% dplyr::rename(Name=1) %>% tibble::rownames_to_column(var="Class") %>% dplyr::mutate(Class=gsub("[1234567890].*$","",Class)) %>% dplyr::filter(Class %in% mSetEnrich) %>% dplyr::arrange(Class)
      mSetMap <- mSet$dataSet$map.table %>% as.data.frame() %>% dplyr::filter(Match %in% mSetClass$Name) %>% tibble::column_to_rownames(var="Match")
      mSetClass$Des <- paste0(mSetMap[as.vector(mSetClass$Name),"Query"],"_",as.vector(mSetClass$Name))
      mSetClass <- mSetClass %>% distinct(Des, .keep_all = TRUE)  %>% tibble::column_to_rownames(var = "Des") %>% dplyr::select(Class)
    }else{
      mSet <- NULL
      mSetEnrich <- NULL
      mSetClass <- NULL
      mSetMap <- NULL
    }
  }
  return(list(mSet=mSet,mSetEnrich=mSetEnrich,mSetClass=mSetClass,mSetMap=mSetMap))
}
exactandEnrich <- function(type){
  result <- list()
  for(phase in c("MP","S0","S12","S34")){
    if(phase == 'MP'){inGUS <- signTaxG0;}
    if(phase == 'S0'){inGUS <- signTaxG1;}
    if(phase == 'S12'){inGUS <- signTaxG2;}
    if(phase == 'S34'){inGUS <- signTaxG3;}
    
    if(type == 'KO'){
      mid_KO_Corr <- Total_KO_Corr$edge %>% dplyr::filter(to %in% inGUS)
      mid_KO_Corr$Mark <- as.vector(signKOs[as.vector(mid_KO_Corr$from),"V2"])
      mid_KO_Corr$gusMark <- gus36Mark[as.vector(mid_KO_Corr$to)] %>% as.vector()
    }
    if(type == 'MB'){
      mid_KO_Corr <- Total_MB_Corr$edge %>% dplyr::filter(to %in% inGUS)
      mid_KO_Corr$Mark <- as.vector(signMBs[as.vector(mid_KO_Corr$from),"V2"])
      mid_KO_Corr$gusMark <- gus36Mark[as.vector(mid_KO_Corr$to)] %>% as.vector()
    }
    if(type == 'Spe'){
      mid_KO_Corr <- Total_spe_Corr$edge %>% dplyr::filter(to %in% inGUS)
      mid_KO_Corr$Mark <- as.vector(signSps[as.vector(mid_KO_Corr$from),"V2"])
      mid_KO_Corr$gusMark <- gus36Mark[as.vector(mid_KO_Corr$to)] %>% as.vector()
    }
    
    if(phase == 'MP'){
      mid_KO_Corr1 <- mid_KO_Corr[grep("\\+...",as.vector(mid_KO_Corr$Mark)),]
      mid_KO_Corr2 <- mid_KO_Corr[grep("-...",as.vector(mid_KO_Corr$Mark)),]
      mid_KO_Corr <- rbind(mid_KO_Corr1,mid_KO_Corr2)
      midResult <- c()
      for (lineN in 1:nrow(mid_KO_Corr)) {
        mark1 <- substr(mid_KO_Corr[lineN,"Mark"],1,1)
        mark2 <- substr(mid_KO_Corr[lineN,"gusMark"],1,1)
        if(mark1 ==  mark2 & mid_KO_Corr[lineN,"r"] > 0){
          midResult <- c(midResult,lineN)
        }
        if(mark1 !=  mark2 & mid_KO_Corr[lineN,"r"] < 0){
          midResult <- c(midResult,lineN)
        }
      }
      mid_KO_Corr <- mid_KO_Corr[midResult,]
      mid_KO_Corr1 <- mid_KO_Corr[grep("\\+...",as.vector(mid_KO_Corr$Mark)),]
      mid_KO_Corr2 <- mid_KO_Corr[grep("-...",as.vector(mid_KO_Corr$Mark)),]
    }
    if(phase == 'S0'){
      mid_KO_Corr1 <- mid_KO_Corr[grep(".\\+..",as.vector(mid_KO_Corr$Mark)),]
      mid_KO_Corr2 <- mid_KO_Corr[grep(".-..",as.vector(mid_KO_Corr$Mark)),]
      mid_KO_Corr <- rbind(mid_KO_Corr1,mid_KO_Corr2)
      midResult <- c()
      for (lineN in 1:nrow(mid_KO_Corr)) {
        mark1 <- substr(mid_KO_Corr[lineN,"Mark"],2,2)
        mark2 <- substr(mid_KO_Corr[lineN,"gusMark"],2,2)
        if(mark1 ==  mark2 & mid_KO_Corr[lineN,"r"] > 0){
          midResult <- c(midResult,lineN)
        }
        if(mark1 !=  mark2 & mid_KO_Corr[lineN,"r"] < 0){
          midResult <- c(midResult,lineN)
        }
      }
      mid_KO_Corr <- mid_KO_Corr[midResult,]
      mid_KO_Corr1 <- mid_KO_Corr[grep(".\\+..",as.vector(mid_KO_Corr$Mark)),]
      mid_KO_Corr2 <- mid_KO_Corr[grep(".-..",as.vector(mid_KO_Corr$Mark)),]
    }
    if(phase == 'S12'){
      mid_KO_Corr1 <- mid_KO_Corr[grep("..\\+.",as.vector(mid_KO_Corr$Mark)),]
      mid_KO_Corr2 <- mid_KO_Corr[grep("..-.",as.vector(mid_KO_Corr$Mark)),]
      mid_KO_Corr <- rbind(mid_KO_Corr1,mid_KO_Corr2)
      midResult <- c()
      for (lineN in 1:nrow(mid_KO_Corr)) {
        mark1 <- substr(mid_KO_Corr[lineN,"Mark"],3,3)
        mark2 <- substr(mid_KO_Corr[lineN,"gusMark"],3,3)
        if(mark1 ==  mark2 & mid_KO_Corr[lineN,"r"] > 0){
          midResult <- c(midResult,lineN)
        }
        if(mark1 !=  mark2 & mid_KO_Corr[lineN,"r"] < 0){
          midResult <- c(midResult,lineN)
        }
      }
      mid_KO_Corr <- mid_KO_Corr[midResult,]
      mid_KO_Corr1 <- mid_KO_Corr[grep("..\\+.",as.vector(mid_KO_Corr$Mark)),]
      mid_KO_Corr2 <- mid_KO_Corr[grep("..-.",as.vector(mid_KO_Corr$Mark)),]
    }
    if(phase == 'S34'){
      mid_KO_Corr1 <- mid_KO_Corr[grep("...\\+",as.vector(mid_KO_Corr$Mark)),]
      mid_KO_Corr2 <- mid_KO_Corr[grep("...-",as.vector(mid_KO_Corr$Mark)),]
      mid_KO_Corr <- rbind(mid_KO_Corr1,mid_KO_Corr2)
      midResult <- c()
      for (lineN in 1:nrow(mid_KO_Corr)) {
        mark1 <- substr(mid_KO_Corr[lineN,"Mark"],4,4)
        mark2 <- substr(mid_KO_Corr[lineN,"gusMark"],4,4)
        if(mark1 ==  mark2 & mid_KO_Corr[lineN,"r"] > 0){
          midResult <- c(midResult,lineN)
        }
        if(mark1 !=  mark2 & mid_KO_Corr[lineN,"r"] < 0){
          midResult <- c(midResult,lineN)
        }
      }
      mid_KO_Corr <- mid_KO_Corr[midResult,]
      mid_KO_Corr1 <- mid_KO_Corr[grep("...\\+",as.vector(mid_KO_Corr$Mark)),]
      mid_KO_Corr2 <- mid_KO_Corr[grep("...-",as.vector(mid_KO_Corr$Mark)),]
    }
    
    if(nrow(mid_KO_Corr1) > 0){mid_KO_Corr1$Type <- 'Up'}
    if(nrow(mid_KO_Corr2)){mid_KO_Corr2$Type <- 'Down'}
    mid_KO_Corr <- rbind(mid_KO_Corr1,mid_KO_Corr2)
    if(nrow(mid_KO_Corr)){mid_KO_Corr$Stage <- phase}
    if(type == 'KO'){
      midUp <- enrichMBandKO(mid_KO_Corr1,"KO")
      midDown <- enrichMBandKO(mid_KO_Corr2,"KO")
      midTotal <- enrichMBandKO(mid_KO_Corr,"KO")
      result[[phase]] = list(corr=mid_KO_Corr,upEnrich=midUp,downEnrich=midDown,totalEnrich=midTotal)
    }else{
      result[[phase]] = list(corr=mid_KO_Corr)
    }
  }
  total_corr <- rbind(result$MP$corr,result$S0$corr,result$S12$corr,result$S34$corr)
  if(type == 'KO'){
    Total_Enrich <- enrichMBandKO(total_corr,"KO")
  }else{
    if(type == 'MB'){
      Total_Enrich <- enrichMBandKO(total_corr,"MB")
    }else{
      Total_Enrich <- NA
    }
  }
  return(list(total=list(corr=total_corr,enrich=Total_Enrich),stage=result))
}
CorrGUSs2KOs <- function(signSpecies,dataGUS,phase,cutoff=0.2){
  signSps <- signKOs %>% dplyr::mutate(TT='test')
  if(phase == 'MP'){
    whatCut <- meanKOs %>% dplyr::filter(Healthy_Mean > 1e-5 | MP_Mean > 1e-5) %>% dplyr::filter(Healthy_Median > 0 | MP_Median > 0) %>% rownames()
    inSpeciesUP <- intersect(signSps[grep("\\+...",as.vector(signSps$V2)),] %>% rownames(),whatCut)
    inSpeciesDown <- intersect(signSps[grep("-...",as.vector(signSps$V2)),] %>% rownames(),whatCut)
    inSamples <- group %>% dplyr::filter(Stage == 'MP') %>% rownames()
    inGUSsUP <- signSpecies$Up$MP
    inGUSsDown <- signSpecies$Down$MP
  }
  if(phase == 'S0'){
    whatCut <- meanKOs %>% dplyr::filter(Healthy_Mean > 1e-5 | S0_Mean > 1e-5) %>% dplyr::filter(Healthy_Median > 0 | S0_Median > 0) %>% rownames()
    inSpeciesUP <- intersect(signSps[grep(".\\+..",as.vector(signSps$V2)),] %>% rownames(),whatCut)
    inSpeciesDown <- intersect(signSps[grep(".-..",as.vector(signSps$V2)),] %>% rownames(),whatCut)
    inSamples <- group %>% dplyr::filter(Stage == 'S0') %>% rownames()
    inGUSsUP <- signSpecies$Up$S0
    inGUSsDown <- signSpecies$Down$S0
  }
  if(phase == 'S12'){
    whatCut <- meanKOs %>% dplyr::filter(Healthy_Mean > 1e-5 | SI_II_Mean > 1e-5) %>% dplyr::filter(Healthy_Median > 0 | SI_II_Median > 0) %>% rownames()
    inSpeciesUP <- intersect(signSps[grep("..\\+.",as.vector(signSps$V2)),] %>% rownames(),whatCut)
    inSpeciesDown <- intersect(signSps[grep("..-.",as.vector(signSps$V2)),] %>% rownames(),whatCut)
    inSamples <- group %>% dplyr::filter(Stage == 'SI_II') %>% rownames()
    inGUSsUP <- signSpecies$Up$S12
    inGUSsDown <- signSpecies$Down$S12
  }
  if(phase == 'S34'){
    whatCut <- meanKOs %>% dplyr::filter(Healthy_Mean > 1e-5 | SIII_IV_Mean > 1e-5) %>% dplyr::filter(Healthy_Median > 0 | SIII_IV_Median > 0) %>% rownames()
    inSpeciesUP <- intersect(signSps[grep("...\\+",as.vector(signSps$V2)),] %>% rownames(),whatCut)
    inSpeciesDown <- intersect(signSps[grep("...-",as.vector(signSps$V2)),] %>% rownames(),whatCut)
    inSamples <- group %>% dplyr::filter(Stage == 'SIII_IV') %>% rownames()
    inGUSsUP <- signSpecies$Up$S34
    inGUSsDown <- signSpecies$Down$S34
  }
  inGUSs <- c(inGUSsUP,inGUSsDown)
  inSpecies <- c(inSpeciesUP,inSpeciesDown)
  inSamples <- c(inSamples,group %>% dplyr::filter(Group == 'Healthy') %>% rownames())
  
  dataMKS <- KOAbun[inSpecies,inSamples]
  #cutoff: 10%
  mid <- dataMKS
  mid[mid > 0] <- 1
  mark <- rowSums(mid) >= ncol(dataMKS) * 0.1
  dataMKS <- dataMKS[mark,]
  #cutoff: mean > 1e-5
  #meanMKS <- apply(dataMKS,1,mean)
  #dataMKS <- dataMKS[meanMKS[meanMKS > 1e-5] %>% names(),]
  print(dim(dataMKS))
  
  
  dataGUS <- dataGUS[inGUSs,inSamples]
  CTR.mat <- corr.test(t(dataMKS),t(dataGUS),method="spearman",adjust = "fdr")
  edgeOUTV2 <- data.frame(from=NA,to=NA,r=NA,p=NA)
  i <- 1
  for (mid.tax in colnames(CTR.mat$r)) {
    for (mid.rna in rownames(CTR.mat$r)) {
      if(CTR.mat$p.adj[mid.rna,mid.tax] < 0.05 & !is.na(CTR.mat$p.adj[mid.rna,mid.tax]) & !is.na(CTR.mat$r[mid.rna,mid.tax])){
        edgeOUTV2[i,] <- c(mid.rna,mid.tax,CTR.mat$r[mid.rna,mid.tax],CTR.mat$p.adj[mid.rna,mid.tax])
        i <- i+1
      }
    }
  }
  edgeOUTV2$r <- as.numeric(edgeOUTV2$r)
  edgeOUTV2$p <- as.numeric(edgeOUTV2$p)
  edgeOUTV2 <- edgeOUTV2 %>% dplyr::filter(abs(r) > cutoff)
  
  if(nrow(edgeOUTV2) > 0){
    edgeOUTV2$Mark <- as.vector(signSps[as.vector(edgeOUTV2$from),'V2'])
    edgeOUTV2$Name <- as.vector(KOinfo[as.vector(edgeOUTV2$from),'Name'])
    edgeOUTV2$Description <- as.vector(KOinfo[as.vector(edgeOUTV2$from),'Description'])
    midResult <- rbind(
      edgeOUTV2 %>% dplyr::filter(from %in% inSpeciesUP & to %in% inGUSsUP) %>% dplyr::filter(r > 0) %>% dplyr::mutate(Type = 'Up'),
      edgeOUTV2 %>% dplyr::filter(from %in% inSpeciesDown & to %in% inGUSsDown) %>% dplyr::filter(r > 0) %>% dplyr::mutate(Type = 'Down'),
      edgeOUTV2 %>% dplyr::filter(from %in% inSpeciesUP & to %in% inGUSsDown) %>% dplyr::filter(r < 0) %>% dplyr::mutate(Type = 'Up'),
      edgeOUTV2 %>% dplyr::filter(from %in% inSpeciesDown & to %in% inGUSsUP) %>% dplyr::filter(r < 0) %>% dplyr::mutate(Type = 'Down')
    ) %>% dplyr::mutate(Stage = phase)
    midUp <- enrichMBandKO(midResult %>% dplyr::filter(Type == 'Up'),"KO")
    midDown <- enrichMBandKO(midResult %>% dplyr::filter(Type == 'Down'),"KO")
    midTotal <- enrichMBandKO(midResult,"KO")
  }else{
    midResult <- NA
    midUp <- NA
    midDown <- NA
    midTotal <- NA
  }
  return(list(mat=CTR.mat,corr=midResult,enrich=midTotal,enrichUp=midUp,enrichDown=midDown))
}
CorrGUSs2MBs <- function(signSpecies,dataGUS,phase,cutoff=0.2){
  signSps <- signMBs %>% dplyr::mutate(TT='test')
  if(phase == 'MP'){
    mid_KO_Corr1 <- signSps[grep("\\+...",as.vector(signSps$V2)),] 
    mid_KO_Corr2 <- signSps[grep("-...",as.vector(signSps$V2)),]
    inSpecies <- rbind(mid_KO_Corr1,mid_KO_Corr2)
    inSamples <- group %>% dplyr::filter(Stage == 'MP') %>% rownames()
    inGUSsUP <- signSpecies$Up$MP
    inGUSsDown <- signSpecies$Down$MP
  }
  if(phase == 'S0'){
    mid_KO_Corr1 <- signSps[grep(".\\+..",as.vector(signSps$V2)),]
    mid_KO_Corr2 <- signSps[grep(".-..",as.vector(signSps$V2)),]
    inSpecies <- rbind(mid_KO_Corr1,mid_KO_Corr2)
    inSamples <- group %>% dplyr::filter(Stage == 'S0') %>% rownames()
    inGUSsUP <- signSpecies$Up$S0
    inGUSsDown <- signSpecies$Down$S0
  }
  if(phase == 'S12'){
    mid_KO_Corr1 <- signSps[grep("..\\+.",as.vector(signSps$V2)),]
    mid_KO_Corr2 <- signSps[grep("..-.",as.vector(signSps$V2)),]
    inSpecies <- rbind(mid_KO_Corr1,mid_KO_Corr2)
    inSamples <- group %>% dplyr::filter(Stage == 'SI_II') %>% rownames()
    inGUSsUP <- signSpecies$Up$S12
    inGUSsDown <- signSpecies$Down$S12
  }
  if(phase == 'S34'){
    mid_KO_Corr1 <- signSps[grep("...\\+",as.vector(signSps$V2)),]
    mid_KO_Corr2 <- signSps[grep("...-",as.vector(signSps$V2)),]
    inSpecies <- rbind(mid_KO_Corr1,mid_KO_Corr2)
    inSamples <- group %>% dplyr::filter(Stage == 'SIII_IV') %>% rownames()
    inGUSsUP <- signSpecies$Up$S34
    inGUSsDown <- signSpecies$Down$S34
  }
  inGUSs <- c(inGUSsUP,inGUSsDown)
  inSpeciesUP <- rownames(mid_KO_Corr1)
  inSpeciesDown <- rownames(mid_KO_Corr2)
  inSamples <- intersect(c(inSamples,group %>% dplyr::filter(Group == 'Healthy') %>% rownames()),colnames(metaAbun))
  
  dataMKS <- metaAbun[rownames(inSpecies),inSamples]
  dataGUS <- dataGUS[inGUSs,inSamples]
  CTR.mat <- corr.test(t(dataMKS),t(dataGUS),method="spearman",adjust = "fdr")
  edgeOUTV2 <- data.frame(from=NA,to=NA,r=NA,p=NA)
  i <- 1
  for (mid.tax in colnames(CTR.mat$r)) {
    for (mid.rna in rownames(CTR.mat$r)) {
      if(CTR.mat$p.adj[mid.rna,mid.tax] < 0.05 & !is.na(CTR.mat$p.adj[mid.rna,mid.tax]) & !is.na(CTR.mat$r[mid.rna,mid.tax])){
        edgeOUTV2[i,] <- c(mid.rna,mid.tax,CTR.mat$r[mid.rna,mid.tax],CTR.mat$p.adj[mid.rna,mid.tax])
        i <- i+1
      }
    }
  }
  edgeOUTV2$r <- as.numeric(edgeOUTV2$r)
  edgeOUTV2$p <- as.numeric(edgeOUTV2$p)
  edgeOUTV2 <- edgeOUTV2 %>% dplyr::filter(abs(r) > cutoff)
  
  if(nrow(edgeOUTV2) > 0){
    edgeOUTV2$Mark <- as.vector(signSps[as.vector(edgeOUTV2$from),'V2'])
    #edgeOUTV2$Name <- as.vector(KOinfo[as.vector(edgeOUTV2$from),'Name'])
    #edgeOUTV2$Description <- as.vector(KOinfo[as.vector(edgeOUTV2$from),'Description'])
    midResult <- rbind(
      edgeOUTV2 %>% dplyr::filter(from %in% inSpeciesUP & to %in% inGUSsUP) %>% dplyr::filter(r > 0) %>% dplyr::mutate(Type = 'Up'),
      edgeOUTV2 %>% dplyr::filter(from %in% inSpeciesDown & to %in% inGUSsDown) %>% dplyr::filter(r > 0) %>% dplyr::mutate(Type = 'Down'),
      edgeOUTV2 %>% dplyr::filter(from %in% inSpeciesUP & to %in% inGUSsDown) %>% dplyr::filter(r < 0) %>% dplyr::mutate(Type = 'Up'),
      edgeOUTV2 %>% dplyr::filter(from %in% inSpeciesDown & to %in% inGUSsUP) %>% dplyr::filter(r < 0) %>% dplyr::mutate(Type = 'Down')
    ) %>% dplyr::mutate(Stage = phase)
  }else{
    midResult <- NA
  }
  return(list(mat=CTR.mat,corr=midResult))
}
group_mean <- function(group.subgroup,data,colInfo="Group"){
  tmpg=unique(group.subgroup[,colInfo])
  tmpn=length(tmpg)
  rows=nrow(data)
  result <- matrix(nrow=rows,ncol=tmpn*3)
  name <- NULL
  for(i in tmpg){ name <- c(name,paste(i,'Median',sep='_'),paste(i,'Mean',sep='_'),paste(i,'SD',sep='_'))    }
  colnames(result)=name
  for(i in tmpg){
    tmpdata=data[,rownames(group.subgroup[which(group.subgroup[,colInfo]==i),,drop=F])]
    name <- c(paste(i,'Median',sep='_'),paste(i,'Mean',sep='_'),paste(i,'SD',sep='_'))
    for (j in 1:rows){
      tmp.median=median(as.numeric(tmpdata[j,]))
      tmp.mean=mean(as.numeric(tmpdata[j,]))
      tmp.SD=sd(as.numeric(tmpdata[j,]))
      result[j,name] <- c(tmp.median,tmp.mean,tmp.SD)
    }
  }
  rownames(result) <- rownames(data)
  result <- data.frame(result)
  return(result)
}

##### Get Data #####
GUSabun_TPM <- read.csv("00.rawdata/GUSabun_TPM.csv",header = T,row.names = 1)
group <- read.csv("00.rawdata/group.csv",header = T,row.names = 1)
groupV2 <- group %>% dplyr::mutate(ID=Subject_ID) %>% tibble::remove_rownames() %>% tibble::column_to_rownames('ID')

signGUSs <- readRDS("00.rawdata/signGUSs.RDS")
signKOs <- read.table("00.rawdata/supp/supp.signKOs.mark.txt",sep="\t",header = F,row.names = 1,check.names = F) %>% dplyr::filter(V2 != "____")
signMBs <- read.table("00.rawdata/supp/supp.signMBs.mark.txt",sep="\t",header = F,row.names = 1,check.names = F) %>% dplyr::filter(V2 != "____")

metaAbun <- read.table("00.rawdata/supp/Supp.metabolites.txt",sep = "\t",header = T,quote = "",row.names = 1,check.names = F)
rownames(metaAbun) <- gsub("'","_",gsub("\"$","",gsub("^\"","",rownames(metaAbun))))
metaAbun <- metaAbun[rownames(signMBs),]
metaAbun <- metaAbun[,as.character(group[intersect(groupV2[colnames(metaAbun),"Sample_ID"],rownames(group)),"Subject_ID"])]
colnames(metaAbun) <- groupV2[colnames(metaAbun),"Sample_ID"]

KOAbun <- read.table("00.rawdata/supp/Supp.ko.txt",sep = "\t",header = T,row.names = 1,check.names = F,quote = "")[rownames(signKOs),as.character(group$Subject_ID)]
colnames(KOAbun) <- groupV2[colnames(KOAbun),"Sample_ID"]

KOinfo <- read.table("00.rawdata/KOinfo/DB.koInfo.txt",sep = "\t",header = F,row.names = 1,quote = "")
colnames(KOinfo) <- c("Name","Description")

KO2names <- read.table("00.rawdata/KOinfo/DB.koInfo.txt",sep="\t",header = F,quote = "",row.names = 1)
pathway2category <- read.table("00.rawdata/KOinfo/pathway2category.txt",sep = "\t",header = T,row.names = 1)

term2name <- read.csv("00.rawdata/KOinfo/formatDB.pathwayInfo.txt",sep="\t",header = F)
colnames(term2name) <- c("PathwayID","PathwayName")
term2name2 <- term2name %>% tibble::column_to_rownames(var='PathwayID')

term2gene <- read.table("00.rawdata/KOinfo/formatDB.ko2pathway.txt",sep="\t",header = F)[,c("V2","V1")]
colnames(term2gene) <- c("PathwayID","KOID")
term2gene <- term2gene %>% dplyr::filter(PathwayID %in% setdiff(unique(term2gene$PathwayID),c("ko01100","ko01110","ko01120","ko01200","ko01210","ko01212","ko01230","ko01232","ko01250","ko01240","ko01220")))
term2gene$KOinfo <- paste0(as.vector(KOinfo[as.vector(term2gene$KOID),"Name"]),": ",as.vector(KOinfo[as.vector(term2gene$KOID),"Description"]))
term2gene$Pathinfo <- as.vector(term2name2[as.vector(term2gene$PathwayID),"PathwayName"])

meanKOs <- group_mean(group,KOAbun,'Stage')
meanMBs <- group_mean(group[colnames(metaAbun),],metaAbun,'Stage')

##### Calculate of correlation ######################
library(psych)
library(clusterProfiler)
Total_KO_Enrich <- list()
Total_MB_Enrich <- list()
for (midStage in c('MP','S0','S12','S34')) {
  gus2KO <- CorrGUSs2KOs(signGUSs$signSpecies,GUSabun_TPM,midStage)
  gus2MB <- CorrGUSs2MBs(signGUSs$signSpecies,GUSabun_TPM,midStage)
  Total_KO_Enrich$stage[[midStage]] <- gus2KO
  Total_MB_Enrich$stage[[midStage]] <- gus2MB
}

##### recalculate edge for metabolites network #####
calEdge <- function(mid_spe_Corr,output=NA){
  mid_spe_Corr$gusSpecies <- gsub(".GUS([1234567890])*$","",as.vector(mid_spe_Corr$to))
  spe2gusSpe <- table(mid_spe_Corr$from,mid_spe_Corr$gusSpecies)
  recalEdge <- data.frame(gusSpecies=NA,from=NA,mean=NA,n=NA,percent=NA)
  for (mid.spe in rownames(spe2gusSpe)) {
    for(mid.gus in colnames(spe2gusSpe)){
      mid.result <- mid_spe_Corr %>% dplyr::filter(from == mid.spe & gusSpecies == mid.gus) %>% dplyr::group_by(gusSpecies,from) %>% dplyr::summarise(mean=mean(r),n=n()) %>% as.data.frame()
      mid.result$percent <- mid.result[,"n"]/as.vector(table(gsub(".GUS([1234567890])*$","",unique(unlist(signGUSs$signSpecies))))[mid.gus])
      recalEdge <- rbind(recalEdge,mid.result) 
    }
  }
  recalEdge <- recalEdge[-1,]
  recalEdge <- recalEdge %>% dplyr::mutate(rr=abs(mean),CorrType=ifelse(mean>0,'up','down'))
  return(recalEdge)
}
calEdge_S34 <- calEdge(Total_MB_Enrich$stage$S34$corr) #used as input for cytoscape

##### metaboAnalyst analysis ######################
# performed online with the correlated metabolites: unique(calEdge_S34$from)
# then visualized using the following codes
p1 <- read.table("msea_ora_result_disease.csv",sep = ",",header = T) %>% dplyr::filter(FDR < 0.05) %>% dplyr::mutate(ratio=hits/total,FDRCol=ifelse(FDR<0.05,FDR,NA),data='disease') %>% dplyr::arrange(-ratio) 
p2 <- read.table("msea_ora_result_KEGG.csv",sep = ",",header = T) %>% dplyr::arrange(-hits) %>% dplyr::filter(hits > 1) %>% dplyr::mutate(ratio=hits/total,FDRCol=ifelse(FDR<0.05,FDR,NA),data='KEGG') %>% dplyr::arrange(-ratio) 
rbind(p1,p2) %>% ggplot(aes(-log10(FDR),reorder(X,-FDR),color=FDRCol,size=ratio)) + geom_point(alpha=0.8) + ylab("") + scale_color_gradient(low = '#e06663',high="#327eba",na.value = 'grey70') + theme_bw()  + ggforce::facet_row(vars(data), scales = 'free', space = 'free')

##### enrichment based on all correlated KOs #####
inKOs <- unique((rbind(Total_KO_Enrich$stage$MP$corr,Total_KO_Enrich$stage$S0$corr,Total_KO_Enrich$stage$S12$corr,Total_KO_Enrich$stage$S34$corr))$from)
inKOsEnrich <- enricher(inKOs,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 1,pAdjustMethod = "BH",qvalueCutoff = 1)
midEnrich <- as.data.frame(inKOsEnrich) %>% dplyr::arrange(-Count) %>% dplyr::filter(Count > 1)
midResult <- NULL
for (midNum in 1:nrow(midEnrich)) {
  midData <- KOinfo[unlist(strsplit(midEnrich[midNum,'geneID'],'\\/')),]
  midData$Mark <- signKOs[unlist(strsplit(midEnrich[midNum,'geneID'],'\\/')),]
  midData$Path <- midEnrich[midNum,'Description']
  midResult <- rbind(midResult,midData)
}
test <- enricher(inKOs,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.05,pAdjustMethod = "fdr",qvalueCutoff = 0.05) 
enrichplot::dotplot(test)

##### venn analysis of correlated KOs #####
library(VennDiagram)
library(grid)
test <- venn.diagram(list(
  MP=unique(Total_KO_Enrich$stage$MP$corr$from),
  S0=unique(Total_KO_Enrich$stage$S0$corr$from),
  S12=unique(Total_KO_Enrich$stage$S12$corr$from),
  S34=unique(Total_KO_Enrich$stage$S34$corr$from)
),filename = NULL,fill = c('yellow', 'purple',"red","blue"))
grid.draw(test)

##### visualization based on self-organized information #####
selfKOs <- read.table("00.rawdata/Liter/Summary of KO.txt",sep = "\t",header = T,quote = "")
selfKOs %>% dplyr::mutate(MP=substr(Mark,1,1),S0=substr(Mark,2,2),S12=substr(Mark,3,3),S34=substr(Mark,4,4))
orderPath <- unique(selfKOs$Path)
orderCategory <- selfKOs %>% dplyr::select(Path,Category) %>% dplyr::distinct(Path,.keep_all = T) %>% tibble::column_to_rownames('Path')
selfKOs <- selfKOs %>% dplyr::mutate(MP=substr(Mark,1,1),S0=substr(Mark,2,2),S12=substr(Mark,3,3),S34=substr(Mark,4,4)) %>% dplyr::select(Path,MP,S0,S12,S34) %>% reshape2::melt("Path") %>% dplyr::mutate(up=ifelse(value == '+',1,0),down=ifelse(value == '-',1,0),value = ifelse(value == '_',0,1)) %>% dplyr::group_by(Path,variable) %>% dplyr::summarise(signNum=sum(value),Up=sum(up),Down=sum(down)) %>% as.data.frame() %>% dplyr::mutate(Color=ifelse(Up > 0 & Down > 0,'Mix',ifelse(Up > 0 & Down == 0,'Up',ifelse(Down > 0 & Up == 0,'Down','Nosign'))))
selfKOs$Path <- factor(as.vector(selfKOs$Path),levels=orderPath)
selfKOs$Category <- factor(as.vector(orderCategory[as.vector(selfKOs$Path),"Category"]),levels=unique(orderCategory$Category))
selfKOs$variable <- factor(as.vector(selfKOs$variable),levels = rev(c("MB","MP","S0","S12","S34")))
selfMBs <- read.table("00.rawdata/Liter/Summary of MB.txt",sep = "\t",header = T,quote = "") %>% dplyr::mutate(up=ifelse(Type == 'Up',1,0),down=ifelse(Type == 'Down',1,0),value = ifelse(Type == 'Up' | Type == 'Down',1,0)) %>% dplyr::group_by(Pathway) %>% dplyr::summarise(signNum=sum(value),Up=sum(up),Down=sum(down)) %>% as.data.frame() %>% dplyr::mutate(Color=ifelse(Up > 0 & Down > 0,'Mix',ifelse(Up > 0 & Down == 0,'Up',ifelse(Down > 0 & Up == 0,'Down','Nosign')))) %>% dplyr::mutate(Path=Pathway,variable='MB') %>% dplyr::select(-Pathway)
selfMBs$Category <- as.vector(orderCategory[as.vector(selfMBs$Path),"Category"])
selfMBs <- selfMBs[,colnames(selfKOs)] %>% dplyr::filter(Category!='NA')
selfInfo <- rbind(selfMBs,selfKOs) 
selfInfo$Path <- factor(as.vector(selfInfo$Path),levels=orderPath)
selfInfo$Category <- factor(as.vector(selfInfo$Category),levels=unique(orderCategory$Category))
selfInfo %>% ggplot(aes(Path,variable,color=Color,size=signNum)) + geom_point(alpha=0.8) + theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = 'bottom') + xlab('') + ylab('') + scale_color_manual(values = c('Down'='navy','Mix'='grey70','Up'='firebrick3')) + ggforce::facet_row(vars(Category), scales = 'free', space = 'free') + ylim(rev(c("MB","MP","S0","S12","S34")))

# ######################################################################
#
## Analysis for gmGUS-species correlation  (Figure 4)
## Author: JR C
## Date: 01.01.2025
#
# ######################################################################
library(dplyr)
library(ggplot2)

##### Function #####
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

signSps <- read.table("00.rawdata/supp/supp.signSpecies.mark.txt",sep="\t",header = F,row.names = 1) %>% dplyr::filter(V2 != "____")

SpeciesAbun <- read.table("00.rawdata/supp/Supp.Silva.txt",sep = "\t",header = T,row.names = 1,check.names = F,quote = "")[rownames(signSps),]
SpeciesAbun <- SpeciesAbun[,as.character(group[intersect(groupV2[colnames(SpeciesAbun),"Sample_ID"],rownames(group)),"Subject_ID"])]
colnames(SpeciesAbun) <- groupV2[colnames(SpeciesAbun),"Sample_ID"]

meanSpecies <- group_mean(group[colnames(SpeciesAbun),],SpeciesAbun,'Stage')

GlcAList <- data.frame(Tax=unique(gsub("^s__","",as.vector((read.table("00.rawdata/GlcA.LCA.output.m8.tax.xls",sep = "\t",header = T,row.names = 1))$S))),Mark='GlcA') %>% tibble::column_to_rownames('Tax')
LCA196List <- read.table("00.rawdata/Liter/summary of species.txt",sep = "\t",header = T,row.names = 1)#[corrSpecies,]
GUSlist <- read.csv("00.rawdata/GUSsStat.csv",header = T,row.names = 1)
GUSlist <- data.frame(Tax=unique(as.vector(GUSlist$species)),Mark='GUS') %>% tibble::column_to_rownames(var='Tax')

signGUSs <- readRDS("00.rawdata/signGUSs.RDS")
ggtreeOrderGUS <- c(
  "Bac_cellulosilyticus.GUS2",
  "Bac_cellulosilyticus.GUS18",
  "Bac_cellulosilyticus.GUS5",
  "Bac_cellulosilyticus.GUS3",
  "Bac_cellulosilyticus.GUS15",
  "Bac_cellulosilyticus.GUS1",
  "Bac_nordii.GUS3",
  "Bac_nordii.GUS6",
  "Bac_nordii.GUS5",
  "Bac_nordii.GUS4",
  "Bac_nordii.GUS2",
  "Bac_nordii.GUS1",
  "Bac_ovatus.GUS4",
  "Bac_ovatus.GUS3",
  "Bac_helcogenes.GUS2",
  "Bac_faecium.GUS7",
  "Bac_thetaiotaomicron.GUS4",
  "Par_goldsteinii.GUS6",
  "Par_goldsteinii.GUS3",
  "Par_goldsteinii.GUS2",
  "Par_chongii.GUS7",
  "Par_chongii.GUS1",
  "Ech_strongylocentroti.GUS",
  "Dor_longicatena.GUS2",
  "Dor_longicatena.GUS1",
  "Lac_eligens.GUS2",
  "Rum_torques.GUS",
  "Ent_clostridioformis.GUS2",
  "Hun_hathewayi.GUS1",
  "Osc_bacterium.GUS7",
  "Osc_bacterium.GUS6",
  "Fae_prausnitzii.GUS3",
  "Fae_prausnitzii.GUS7",
  "Fae_prausnitzii.GUS10",
  "Bif_bifidum.GUS2",
  "Bif_bifidum.GUS1",
  "Cel_cellulans.GUS",
  "Stu_stutzeri.GUS"
)

##### calculation of correlations between gmGUS and species #####
library(psych)
CorrGUSs2Species <- function(signSpecies,signSps,dataGUS,dataMKS,phase,cutoff=0.2){
  signSps <- signSps %>% dplyr::mutate(TT='test')
  if(phase == 'MP'){
    whatCut <- meanSpecies %>% dplyr::filter(Healthy_Mean > 1e-5 | MP_Mean > 1e-5) %>% dplyr::filter(Healthy_Median > 0 | MP_Median > 0) %>% rownames()
    inSpeciesUP <- intersect(signSps[grep("\\+...",as.vector(signSps$V2)),] %>% rownames(),whatCut)
    inSpeciesDown <- intersect(signSps[grep("-...",as.vector(signSps$V2)),] %>% rownames(),whatCut)
    inSamples <- group %>% dplyr::filter(Stage == 'MP') %>% rownames()
    inGUSsUP <- signSpecies$Up$MP
    inGUSsDown <- signSpecies$Down$MP
    inGUSs <- c(inGUSsUP,inGUSsDown)
  }
  if(phase == 'S0'){
    whatCut <- meanSpecies %>% dplyr::filter(Healthy_Mean > 1e-5 | S0_Mean > 1e-5) %>% dplyr::filter(Healthy_Median > 0 | S0_Median > 0) %>% rownames()
    inSpeciesUP <- intersect(signSps[grep(".\\+..",as.vector(signSps$V2)),] %>% rownames(),whatCut)
    inSpeciesDown <- intersect(signSps[grep(".-..",as.vector(signSps$V2)),] %>% rownames(),whatCut)
    inSamples <- group %>% dplyr::filter(Stage == 'S0') %>% rownames()
    inGUSsUP <- signSpecies$Up$S0
    inGUSsDown <- signSpecies$Down$S0
    inGUSs <- c(inGUSsUP,inGUSsDown)
  }
  if(phase == 'S12'){
    whatCut <- meanSpecies %>% dplyr::filter(Healthy_Mean > 1e-5 | SI_II_Mean > 1e-5) %>% dplyr::filter(Healthy_Median > 0 | SI_II_Median > 0) %>% rownames()
    inSpeciesUP <- intersect(signSps[grep("..\\+.",as.vector(signSps$V2)),] %>% rownames(),whatCut)
    inSpeciesDown <- intersect(signSps[grep("..-.",as.vector(signSps$V2)),] %>% rownames(),whatCut)
    inSamples <- group %>% dplyr::filter(Stage == 'SI_II') %>% rownames()
    inGUSsUP <- signSpecies$Up$S12
    inGUSsDown <- signSpecies$Down$S12
    inGUSs <- c(inGUSsUP,inGUSsDown)
  }
  if(phase == 'S34'){
    whatCut <- meanSpecies %>% dplyr::filter(Healthy_Mean > 1e-5 | SIII_IV_Mean > 1e-5) %>% dplyr::filter(Healthy_Median > 0 | SIII_IV_Median > 0) %>% rownames()
    inSpeciesUP <- intersect(signSps[grep("...\\+",as.vector(signSps$V2)),] %>% rownames(),whatCut)
    inSpeciesDown <- intersect(signSps[grep("...-",as.vector(signSps$V2)),] %>% rownames(),whatCut)
    inSamples <- group %>% dplyr::filter(Stage == 'SIII_IV') %>% rownames()
    inGUSsUP <- signSpecies$Up$S34
    inGUSsDown <- signSpecies$Down$S34
    inGUSs <- c(inGUSsUP,inGUSsDown)
  }
  inSpecies <- c(inSpeciesUP,inSpeciesDown)
  inSamples <- c(inSamples,group %>% dplyr::filter(Group == 'Healthy') %>% rownames())
  
  dataMKS <- SpeciesAbun[inSpecies,inSamples]
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
  edgeOUTV2$Mark <- as.vector(signSps[as.vector(edgeOUTV2$from),'V2'])
  edgeOUTV2 <- edgeOUTV2 %>% dplyr::filter(abs(r) > cutoff)
  
  midResult <- rbind(
    edgeOUTV2 %>% dplyr::filter(from %in% inSpeciesUP & to %in% inGUSsUP) %>% dplyr::filter(r > 0) %>% dplyr::mutate(Type = 'Up'),
    edgeOUTV2 %>% dplyr::filter(from %in% inSpeciesDown & to %in% inGUSsDown) %>% dplyr::filter(r > 0) %>% dplyr::mutate(Type = 'Down'),
    edgeOUTV2 %>% dplyr::filter(from %in% inSpeciesUP & to %in% inGUSsDown) %>% dplyr::filter(r < 0) %>% dplyr::mutate(Type = 'Up'),
    edgeOUTV2 %>% dplyr::filter(from %in% inSpeciesDown & to %in% inGUSsUP) %>% dplyr::filter(r < 0) %>% dplyr::mutate(Type = 'Down')
  ) %>% dplyr::mutate(Stage = phase)
  return(list(mat=CTR.mat,corr=midResult))
}
Total_Spe_Enrich <- list()
for (midStage in c('MP','S0','S12','S34')) {
  gus2S <- CorrGUSs2Species(signGUSs$signSpecies,signSps,GUSabun_TPM,SpeciesAbun,midStage)
  Total_Spe_Enrich$stage[[midStage]] <- gus2S
}

##### species number (normalized by differential number) across groups ####
rbind(
  table((Total_Spe_Enrich$stage$MP$corr %>% dplyr::select(from,Type) %>% dplyr::distinct(from,.keep_all = T))$Type) %>% as.data.frame() %>% dplyr::mutate(Stage='MP',Total=length(signSps[grep("[-+]...",as.vector(signSps$V2)),])),
  table((Total_Spe_Enrich$stage$S0$corr %>% dplyr::select(from,Type) %>% dplyr::distinct(from,.keep_all = T))$Type) %>% as.data.frame() %>% dplyr::mutate(Stage='S0',Total=length(signSps[grep(".[-+]..",as.vector(signSps$V2)),])),
  table((Total_Spe_Enrich$stage$S12$corr %>% dplyr::select(from,Type) %>% dplyr::distinct(from,.keep_all = T))$Type) %>% as.data.frame() %>% dplyr::mutate(Stage='S12',Total=length(signSps[grep("..[-+].",as.vector(signSps$V2)),])),
  table((Total_Spe_Enrich$stage$S34$corr %>% dplyr::select(from,Type) %>% dplyr::distinct(from,.keep_all = T))$Type) %>% as.data.frame() %>% dplyr::mutate(Stage='S34',Total=length(signSps[grep("...[-+]",as.vector(signSps$V2)),]))
) %>% dplyr::mutate(Perc=Freq/Total) %>% ggplot(aes(Stage,weight=Perc,fill=Var1)) + geom_bar() + theme_classic() + xlab('') + ylab('gus-correlated species') + scale_fill_manual(values = c('navy','firebrick3'))

##### network input: edge.txt #####
calEdge <- function(mid_spe_Corr,output=NA){
  mid_spe_Corr$gusSpecies <- gsub(".GUS([1234567890])*$","",as.vector(mid_spe_Corr$to))
  spe2gusSpe <- table(mid_spe_Corr$from,mid_spe_Corr$gusSpecies)
  recalEdge <- data.frame(gusSpecies=NA,from=NA,mean=NA,n=NA,percent=NA)
  for (mid.spe in rownames(spe2gusSpe)) {
    for(mid.gus in colnames(spe2gusSpe)){
      mid.result <- mid_spe_Corr %>% dplyr::filter(from == mid.spe & gusSpecies == mid.gus) %>% dplyr::group_by(gusSpecies,from) %>% dplyr::summarise(mean=mean(r),n=n()) %>% as.data.frame()
      mid.result$percent <- mid.result[,"n"]/as.vector(table(gsub(".GUS([1234567890])*$","",ggtreeOrderGUS))[mid.gus])
      recalEdge <- rbind(recalEdge,mid.result) 
    }
  }
  recalEdge <- recalEdge[-1,]
  recalEdge <- recalEdge %>% dplyr::mutate(rr=abs(mean),CorrType=ifelse(mean>0,'up','down'))
  return(recalEdge)
}
calEdge_MP <- calEdge(Total_Spe_Enrich$stage$MP$corr)
calEdge_S0 <- calEdge(Total_Spe_Enrich$stage$S0$corr)
calEdge_S12 <- calEdge(Total_Spe_Enrich$stage$S12$corr)
calEdge_S34 <- calEdge(Total_Spe_Enrich$stage$S34$corr)
calEdge_Total <- rbind(
  calEdge_MP %>% dplyr::mutate(from=paste0("MP_",from),gusSpecies=paste0("MP_",gusSpecies)),
  calEdge_S0 %>% dplyr::mutate(from=paste0("S0_",from),gusSpecies=paste0("S0_",gusSpecies)),
  calEdge_S12 %>% dplyr::mutate(from=paste0("S12_",from),gusSpecies=paste0("S12_",gusSpecies)),
  calEdge_S34 %>% dplyr::mutate(from=paste0("S34_",from),gusSpecies=paste0("S34_",gusSpecies))
) # used as input for cytoscape

##### network input: node.txt #####
nodeInfo <- rbind(
  Total_Spe_Enrich$stage$MP$corr %>% dplyr::select(from,Type) %>% dplyr::distinct(from,.keep_all=T) %>% dplyr::mutate(keyF=paste0("MP_",from)),
  Total_Spe_Enrich$stage$S0$corr %>% dplyr::select(from,Type) %>% dplyr::distinct(from,.keep_all=T) %>% dplyr::mutate(keyF=paste0("S0_",from)),
  Total_Spe_Enrich$stage$S12$corr %>% dplyr::select(from,Type) %>% dplyr::distinct(from,.keep_all=T) %>% dplyr::mutate(keyF=paste0("S12_",from)),
  Total_Spe_Enrich$stage$S34$corr %>% dplyr::select(from,Type) %>% dplyr::distinct(from,.keep_all=T) %>% dplyr::mutate(keyF=paste0("S34_",from))
)
nodeInfo <- cbind(
  nodeInfo,
  LCA196List[as.vector(nodeInfo$from),] %>% tibble::remove_rownames()
)
nodeInfo$GlcA <- as.vector(GlcAList[as.vector(nodeInfo$from),])
nodeInfo$GUS <- as.vector(GUSlist[as.vector(nodeInfo$from),]) # used as node information for cytoscape
nodeInfo <- nodeInfo %>% dplyr::mutate(Group = gsub("_.*","",keyF))

##### GlcA and GUS number #####
nodeInfoMid <- nodeInfo %>% dplyr::filter(GlcA == 'GlcA')
nodeInfoMid2 <- nodeInfo %>% dplyr::filter(GUS == 'GUS')
rbind(
  table(as.vector(nodeInfoMid$Group),as.vector(nodeInfoMid$Type)) %>% as.data.frame() %>% reshape2::dcast(Var1~Var2) %>% tibble::column_to_rownames('Var1') %>% as.data.frame() %>% tibble::rownames_to_column('Group') %>% reshape2::melt(by='Group') %>% dplyr::mutate(Mark='GlcA'),
  table(as.vector(nodeInfoMid2$Group),as.vector(nodeInfoMid2$Type)) %>% as.data.frame() %>% reshape2::dcast(Var1~Var2) %>% tibble::column_to_rownames('Var1') %>% as.data.frame() %>% tibble::rownames_to_column('Group') %>% reshape2::melt(by='Group') %>% dplyr::mutate(Mark='GUS')
) %>% ggplot(aes(Group,weight=value,fill=variable)) + geom_bar() + xlab("") + ylab("") + theme_classic() + scale_fill_manual(values = c('navy','firebrick3')) + facet_wrap('Mark',nrow = 2,ncol = 1,strip.position = 'right') + coord_flip()

##### association number across groups ####
rbind(
  table(calEdge_MP$CorrType) %>% as.data.frame() %>% dplyr::mutate(Stage='MP'),
  table(calEdge_S0$CorrType) %>% as.data.frame() %>% dplyr::mutate(Stage='S0'),
  table(calEdge_S12$CorrType) %>% as.data.frame() %>% dplyr::mutate(Stage='S12'),
  table(calEdge_S34$CorrType) %>% as.data.frame() %>% dplyr::mutate(Stage='S34')
)  %>% ggplot(aes(Stage,weight=Freq,fill=Var1)) + geom_bar() + theme_classic() + xlab('') + ylab('Number of associations') + scale_fill_manual(values = c('navy','firebrick3'))

##### fisher's exact test, figure5DE #####
nodeInfoFisher <- as.data.frame(matrix(nrow = 0,ncol = 10))
for (midType in c('Total')) {
  nodeInfoMid <- nodeInfo #nodeInfo %>% dplyr::filter(Type == midType)
  midN <- nodeInfoMid$from %>% unique() %>% length()
  for (midVar in colnames(nodeInfoMid)[c(4,5,6,8:13)]) {
    for (midkey in unique(as.vector(unlist(nodeInfoMid[midVar])))) {
      if(midkey != "NA" & midkey != ""){
        midM <- (nodeInfoMid %>% dplyr::filter(get(midVar) == midkey))$from %>% unique() %>% length()
        for (midG in unique((as.vector(unlist(nodeInfoMid$Group))))) {
          midn <- nodeInfoMid %>% dplyr::filter(Group == midG) %>% nrow()
          midm <- nodeInfoMid %>% dplyr::filter(get(midVar) == midkey & Group == midG) %>% nrow()
          if(midm > 0){
            midp <- phyper(midm-1,midM, midN-midM, midn, lower.tail=FALSE)
            d <- data.frame(gene.not.interest=c(midM-midm, midN-midM-midn+midm), gene.in.interest=c(midm, midn-midm))
            row.names(d) <- c("In_category", "not_in_category")
            midpV2 <- (fisher.test(d,alternative = 'two.sided'))$p.value
            nodeInfoFisher <- rbind(nodeInfoFisher,c(midType,midVar,midkey,midG,midp,midpV2,midN,midM,midn,midm))
          }else{
            nodeInfoFisher <- rbind(nodeInfoFisher,c(midType,midVar,midkey,midG,NA,NA,midN,midM,midn,midm))
          }
        }
      }else{
        nodeInfoFisher <- rbind(nodeInfoFisher,c(midType,midVar,midkey,midG,NA,NA,midN,midM,midn,midm))
      }
    }
  }
}
colnames(nodeInfoFisher) <- c('Type','Feature',"Level","Group","pValue","pValueV2","N","M","n","m")
nodeInfoFisher$GeneRatio <- as.numeric(nodeInfoFisher$m)/as.numeric(nodeInfoFisher$n)
nodeInfoFisher$BgRatio <- as.numeric(nodeInfoFisher$M)/as.numeric(nodeInfoFisher$N)
nodeInfoFisher$m <- as.numeric(nodeInfoFisher$m)
nodeInfoFisher$pValue <- as.numeric(nodeInfoFisher$pValue)
nodeInfoFisher %>% dplyr::filter(Feature == 'Genus' & Level %in% as.vector((nodeInfoFisher %>% dplyr::filter(Feature == 'Genus' & pValue < 0.05))$Level)) %>% dplyr::mutate(pValue=ifelse(pValue > 0.05,NA,pValue)) %>% ggplot(aes(Group,Level,color=pValue,size=GeneRatio)) + geom_point() + xlab("") + ylab("") + scale_colour_gradient2(low = 'firebrick3',mid='blue',high = 'blue',midpoint = 0.05,na.value = 'grey70') + ylim(rev(c('Lachnospira','Phocaeicola','Prevotella','Parabacteroides','Bacteroides','Bifidobacterium','Ruminiclostridium','Alistipes','Acetivibrio','Fusobacterium'))) + theme_classic()
nodeInfoFisher %>% dplyr::filter(Feature == 'Level.Last' & Level %in% c('Primary degrader','Pathogen','Butyrate producers','Lactic acid bacterium','Acetogen','Mucin degrader','PBA converting bacterium','Equol-producing bacterium','Sulphate reducer','Putrescine-fermenting bacterium','Flavonoid-degrading bacterium')) %>% ggplot(aes(Group,reorder(Level,m),color=log2(m),size=GeneRatio))  + geom_point() + xlab("") + ylab("") + scale_colour_gradient(low = 'blue',high = 'firebrick3',na.value = 'grey70') + theme_classic()

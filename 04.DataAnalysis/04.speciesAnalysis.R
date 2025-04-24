# ######################################################################
#
## Species-level analysis (Figure 2FGH)
## Author: JR C
## Date: 01.01.2025
#
# ######################################################################
library(dplyr)
library(ggplot2)
library(ggpubr)

##### Get Data #####
GUSsStat <- read.csv("00.rawdata/GUSsStat.csv",header = T,row.names = 1)
GUSabun_TPM <- read.csv("00.rawdata/GUSabun_TPM.csv",header = T,row.names = 1)
group <- read.csv("00.rawdata/group.csv",header = T,row.names = 1)
analysis.list.info <- list(
  compaired = list(
    "Stage"=list(c("Healthy" ,"MP"),c("Healthy" ,"S0"),c("Healthy" ,"SI_II"),c("Healthy" ,"SIII_IV"))
  ),
  levels = list(
    "Stage"=c("Healthy" ,"MP","S0","SI_II","SIII_IV")
  ),
  colors = list(
    "Stage"=c("Healthy" ="#cecccb","MP"="#55B7E6","S0"="#193E8F","SI_II"="#F09739","SIII_IV"="#E53528")
  )
)

##### Functions #####
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

##### get species abundance (cumulative GUS abundance) #####
sel.taxes <- (table(GUSsStat[,'species']) %>% as.data.frame() %>% dplyr::filter(Freq > 1) %>% dplyr::arrange(-Freq))$Var1 %>% as.vector()
species.abun <- matrix(nrow=length(sel.taxes),ncol = ncol(GUSabun_TPM))
rownames(species.abun) <-sel.taxes
colnames(species.abun) <- colnames(GUSabun_TPM)
for (mid.sp in sel.taxes) {
  species.abun[mid.sp,] <- GUSabun_TPM[GUSsStat %>% dplyr::filter(species == mid.sp) %>% rownames(),] %>% colSums() %>% as.vector()
}

##### difference of abundance #####
wilcoxon.FDR.TAX <- function(data,group){
      
      data=data[,rownames(group)]
      tmpg=unique(group$Stage)
      species.abun <- data

      #statistics
      result <- group_mean(group,data,'Stage')
      
      #wilcoxon
      for (i in c('Healthy')) {
        for (j in c('MP','S0',"SI_II","SIII_IV")){
          #wilcoxon test
          sample1=rownames(group[which(group[,'Stage']==i),,drop=F])
          sample2=rownames(group[which(group[,'Stage']==j),,drop=F])
          wilcox.01.p <- apply(data,1,function(x) {if(sum(x[c(sample1,sample2)]) == 0){return('NA')}else{test <- wilcox.test(x[sample1],x[sample2],conf.int = T);p <- test$p.value;return(p);}}) %>% as.numeric()
          wilcox.01.q <- rep("NA",length(wilcox.01.p))
          wilcox.01.q[which(wilcox.01.p < 0.05)] <- p.adjust(wilcox.01.p[which(wilcox.01.p < 0.05)], method = "fdr")
          wilname.p <- data.frame(p=wilcox.01.p,q=wilcox.01.q)
          colnames(wilname.p) <- paste0(c("wilcox.test.p","wilcox.test.q"),paste0("(",i,' VS ',j,")"))
          result <- cbind(result,wilname.p)
        }
      }
      
      
      # sign species
      cutoff <- 61.6 #total number: 616
      cutFC <- log2(2)
      signSpecies = list()
      
      signTaxMP <- result %>% dplyr::mutate(log2FC=ifelse(MP_Mean > Healthy_Mean,log2(MP_Mean/Healthy_Mean),-log2(Healthy_Mean/MP_Mean))) %>% dplyr::mutate(log2FC=ifelse(!is.finite(log2FC),ifelse(MP_Mean>Healthy_Mean,10,-10),log2FC))
      signTax <- signTaxMP %>% dplyr::filter(`wilcox.test.p(Healthy VS MP)` < 0.05) 
      data <- species.abun[signTax %>% dplyr::arrange(-log2FC) %>% rownames(),]
      data[data>0] <- 1
      signTaxG0 <- intersect(signTax %>% dplyr::filter(abs(log2FC) >= cutFC) %>% rownames(),rownames(data)[rowSums(data) >= cutoff])
      signSpecies[['Up']][['MP']] <- intersect(signTax %>% dplyr::filter(log2FC >= cutFC) %>% rownames(),rownames(data)[rowSums(data) >= cutoff])
      signSpecies[['Down']][['MP']] <- intersect(signTax %>% dplyr::filter(log2FC <= -cutFC) %>% rownames(),rownames(data)[rowSums(data) >= cutoff])
      
      signTaxS0 <- result %>% dplyr::mutate(log2FC=ifelse(S0_Mean > Healthy_Mean,log2(S0_Mean/Healthy_Mean),-log2(Healthy_Mean/S0_Mean))) %>% dplyr::mutate(log2FC=ifelse(!is.finite(log2FC),ifelse(S0_Mean>Healthy_Mean,10,-10),log2FC))
      signTax <- signTaxS0 %>% dplyr::filter(`wilcox.test.q(Healthy VS S0)` < 0.05) 
      data <- species.abun[signTax %>% dplyr::arrange(-log2FC) %>% rownames(),]
      data[data>0] <- 1
      signTaxG1 <- intersect(signTax %>% dplyr::filter(abs(log2FC) >= cutFC) %>% rownames(),rownames(data)[rowSums(data) >= cutoff])
      signSpecies[['Up']][['S0']] <- intersect(signTax %>% dplyr::filter(log2FC >= cutFC) %>% rownames(),rownames(data)[rowSums(data) >= cutoff])
      signSpecies[['Down']][['S0']] <- intersect(signTax %>% dplyr::filter(log2FC <= -cutFC) %>% rownames(),rownames(data)[rowSums(data) >= cutoff])

      signTaxS12 <- result %>% dplyr::mutate(log2FC=ifelse(SI_II_Mean > Healthy_Mean,log2(SI_II_Mean/Healthy_Mean),-log2(Healthy_Mean/SI_II_Mean))) %>% dplyr::mutate(log2FC=ifelse(!is.finite(log2FC),ifelse(SI_II_Mean>Healthy_Mean,10,-10),log2FC))
      signTax <- signTaxS12 %>% dplyr::filter(`wilcox.test.q(Healthy VS SI_II)` < 0.05) 
      data <- species.abun[signTax %>% dplyr::arrange(-log2FC) %>% rownames(),]
      data[data>0] <- 1
      signTaxG2 <- intersect(signTax %>% dplyr::filter(abs(log2FC) >= cutFC) %>% rownames(),rownames(data)[rowSums(data) >= cutoff])
      signSpecies[['Up']][['S12']] <- intersect(signTax %>% dplyr::filter(log2FC >= cutFC) %>% rownames(),rownames(data)[rowSums(data) >= cutoff])
      signSpecies[['Down']][['S12']] <- intersect(signTax %>% dplyr::filter(log2FC <= -cutFC) %>% rownames(),rownames(data)[rowSums(data) >= cutoff])

      signTaxS34 <- result %>% dplyr::mutate(log2FC=ifelse(SIII_IV_Mean > Healthy_Mean,log2(SIII_IV_Mean/Healthy_Mean),-log2(Healthy_Mean/SIII_IV_Mean))) %>% dplyr::mutate(log2FC=ifelse(!is.finite(log2FC),ifelse(SIII_IV_Mean>Healthy_Mean,10,-10),log2FC))
      signTax <- signTaxS34 %>% dplyr::filter(`wilcox.test.q(Healthy VS SIII_IV)` < 0.05) 
      data <- species.abun[signTax %>% dplyr::arrange(-log2FC) %>% rownames(),]
      data[data>0] <- 1
      signTaxG3 <- intersect(signTax %>% dplyr::filter(abs(log2FC) >= cutFC) %>% rownames(),rownames(data)[rowSums(data) >= cutoff])
      signSpecies[['Up']][['S34']] <- intersect(signTax %>% dplyr::filter(log2FC >= cutFC) %>% rownames(),rownames(data)[rowSums(data) >= cutoff])
      signSpecies[['Down']][['S34']] <- intersect(signTax %>% dplyr::filter(log2FC <= -cutFC) %>% rownames(),rownames(data)[rowSums(data) >= cutoff])

      allSign <- unique(c(signTaxG0,signTaxG1,signTaxG2,signTaxG3))
      signFeatureStageFC <- data.frame(Tax=allSign,MP=signTaxMP[allSign,'log2FC'],S0=signTaxS0[allSign,'log2FC'],S12=signTaxS12[allSign,'log2FC'],S34=signTaxS34[allSign,'log2FC']) %>% tibble::column_to_rownames('Tax')
      signFeatureStageMean <- result[allSign,c("Healthy_Mean","MP_Mean","S0_Mean",'SI_II_Mean',"SIII_IV_Mean")]
      
      return(list(testResult=result,signFC=signFeatureStageFC,signMean=signFeatureStageMean,signSpecies=signSpecies))
}
wilcoxonResult <- wilcoxon.FDR.TAX(species.abun,group)

##### Visualization #####
heat.in <- wilcoxonResult$signMean
colnames(heat.in) <- c("HC","MP","S0","S12","S34")
heat.in <- t(scale(t(heat.in),center=T,scale=T))
heat.in %>% as.data.frame() %>% tibble::rownames_to_column('Tax') %>% reshape2::melt(by = 'Tax') %>% dplyr::mutate(Zscore=value) %>% ggplot(aes(variable,Tax,fill=Zscore)) + geom_tile() + xlab("")  + theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "top") + ylab("") + scale_fill_gradient2(low='navy',mid='white',high='firebrick3') + ylim(rownames(wilcoxonResult$signMean))

wilcoxonResult$signFC %>% tibble::rownames_to_column('Tax') %>% reshape2::melt(by='Tax') %>% ggplot(aes(variable,Tax,fill=value)) + geom_tile() + xlab("") + theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "top",axis.text.y = element_blank()) + ylab("") + scale_fill_gradient2(low='navy',mid='white',high='firebrick3')  + ylim(rownames(wilcoxonResult$signMean))

signFeatureStage <- rbind(
  wilcoxonResult$signSpecies$Up %>% unlist() %>% as.data.frame() %>% dplyr::rename(Tax = 1) %>% tibble::rownames_to_column('Group') %>% dplyr::mutate(Group=gsub("S34.*","S34",Group),Type='Up'),
  wilcoxonResult$signSpecies$Down %>% unlist() %>% as.data.frame() %>% dplyr::rename(Tax = 1) %>% tibble::rownames_to_column('Group') %>% dplyr::mutate(Group=gsub("S12.*","S12",Group),Type='Down')
)
signFeatureStage %>% ggplot(aes(Group,Tax,shape=Type,fill=Type)) + geom_point(size=2) + scale_shape_manual(values = c("Down"=25,"Up"=24)) + scale_fill_manual(values = c("Down"="blue","Up"="red")) + xlab("") + ylab("") + ylim(rownames(wilcoxonResult$signMean)) + theme_classic() + theme(legend.position = "top",axis.text.x = element_text(angle = 45,hjust = 1))

##### difference of cumulative GUS number #####
sel.taxes <- (table(GUSsStat[,'species']) %>% as.data.frame() %>% dplyr::filter(Freq > 1) %>% dplyr::arrange(-Freq))$Var1 %>% as.vector()
data <- GUSabun_TPM
data[data > 0] <- 1
species.abun <- matrix(nrow=length(sel.taxes),ncol = ncol(GUSabun_TPM))
rownames(species.abun) <-sel.taxes
colnames(species.abun) <- colnames(GUSabun_TPM)
for (mid.sp in sel.taxes) {
  species.abun[mid.sp,] <- data[GUSsStat %>% dplyr::filter(species == mid.sp) %>% rownames(),] %>% colSums() %>% as.vector()
}
species.abun[rownames(wilcoxonResult$signMean),] %>% as.data.frame() %>% tibble::rownames_to_column('ID') %>% reshape2::melt(by='ID') %>% dplyr::mutate(Group=group[variable,'Stage']) %>% ggplot(aes(Group,value,fill=Group)) + geom_boxplot() + xlab("")  + theme_pubclean() + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "top") + ylab("") + facet_wrap('ID',nrow = 1)  + geom_signif(comparisons = analysis.list.info$compaired$Stage, color="black",step_increase = 0.08, map_signif_level=function(p) if(p < 0.01){'+'}else{if(p < 0.05){'*'}else{sprintf("%.3f",p)}}, test = wilcox.test)
species.abun[c("Bacteroides cellulosilyticus","Bacteroides faecium","Bacteroides nordii"),] %>% as.data.frame() %>% tibble::rownames_to_column('ID') %>% reshape2::melt(by='ID') %>% dplyr::mutate(Group=group[variable,'Stage']) %>% ggplot(aes(Group,value,fill=Group)) + geom_violin() + geom_boxplot(col='black',size=1,width =0.1,outlier.colour = NA) + xlab("")  + theme_pubclean() + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "top") + ylab("") + facet_wrap('ID',nrow = 1)  + geom_signif(comparisons = analysis.list.info$compaired$Stage, color="black",step_increase = 0.08, map_signif_level=function(p) if(p < 0.01){sprintf("%.3f",p)}else{if(p < 0.05){sprintf("%.3f",p)}else{sprintf("%.3f",p)}}, test = wilcox.test) + scale_fill_manual(values=analysis.list.info$colors$Stage)

##### analysis based on mOTU and metaphlan #####
wilcoxon.FDR.mOTU4 <- function(data,group){
  data=data[,rownames(group)]
  tmpg=unique(group$Stage)
  
  #statistics
  result <- group_mean(group,data,'Stage')
  
  #wilcoxon
  for (i in c('Healthy')) {
    for (j in c('MP','S0',"SI_II","SIII_IV")){
      #wilcoxon test
      sample1=rownames(group[which(group[,'Stage']==i),,drop=F])
      sample2=rownames(group[which(group[,'Stage']==j),,drop=F])
      wilcox.01.p <- apply(data,1,function(x) {if(sum(x[c(sample1,sample2)]) == 0){return('NA')}else{test <- wilcox.test(x[sample1],x[sample2],conf.int = T);p <- test$p.value;return(p);}}) %>% as.numeric()
      wilname.p <- data.frame(p=wilcox.01.p)
      colnames(wilname.p) <- paste0(c("wilcox.test.p"),paste0("(",i,' VS ',j,")"))
      result <- cbind(result,wilname.p)
    }
  }
  return(result)
}
#based on mOTU4
dbInfo <- read.table("00.rawdata/mOTUs4/mOTUsv4.0.gtdb.taxonomy.80mv.tsv",sep = "\t",header = T,row.names = 2)
data <- read.table("00.rawdata/mOTUs4/merged.mOTUs.table.txt",sep = "\t",header = T,row.names = 1)
data$Species <- as.vector(dbInfo[rownames(data),'species'])
data <- data %>% as.data.frame() %>% reshape2::melt('Species') %>% dplyr::group_by(Species,variable) %>% dplyr::summarise(meanS=sum(value)) %>% as.data.frame() %>% reshape2::dcast(Species ~ variable)
data <- data %>% as.data.frame() %>% dplyr::filter(Species != '') %>% tibble::column_to_rownames('Species')
selSpecies <- rownames(data)[grep("Bacteroides.faecium|helcogenes|nordii|Subdoligranulum.variabile|bifidum|cellulosilyticus|muciniphila|longicatena|eligens",rownames(data))]
speciesAbun <- data[selSpecies,rownames(group)]
testResult <- wilcoxon.FDR.mOTU4(speciesAbun,group)

#based on metaphlan
data <- read.table("00.rawdata/supp/Supp.metaphlan.txt",sep = "\t",header = T,row.names = 1)
selSpecies <- rownames(data)[grep("Bacteroides.faecium|helcogenes|nordii|Subdoligranulum.variabile|bifidum|cellulosilyticus|muciniphila|longicatena|eligens",rownames(data))]

groupV2 <- group %>% dplyr::mutate(ID=Subject_ID) %>% tibble::remove_rownames() %>% tibble::column_to_rownames('ID')
speciesAbun <- data[selSpecies,]
speciesAbun <- speciesAbun[,paste0("X",group[intersect(groupV2[gsub("^X","",colnames(speciesAbun)),"Sample_ID"],rownames(group)),"Subject_ID"])]
colnames(speciesAbun) <- groupV2[gsub("^X","",colnames(speciesAbun)),"Sample_ID"]
speciesAbun <- speciesAbun[,rownames(group)]
testResult <- wilcoxon.FDR.mOTU4(speciesAbun,group)

##### analysis of CNV #####
rbind(
  read.table("00.rawdata/CNV/100129_Akkermansia muciniphila.CNV.txt",sep = "\t",header = T,row.names = 1) %>% as.data.frame() %>% tibble::rownames_to_column('ID') %>% reshape2::melt(by='ID') %>% dplyr::group_by(variable) %>% dplyr::summarise(Total=sum(value)) %>% dplyr::mutate(Group = group[variable,'Stage']) %>% dplyr::filter(Group != 'NA') %>% dplyr::mutate(Tax='A. muciniphila') %>% as.data.frame(),
  read.table("00.rawdata/CNV/100145_Bifidobacterium bifidum.CNV.txt",sep = "\t",header = T,row.names = 1) %>% as.data.frame() %>% tibble::rownames_to_column('ID') %>% reshape2::melt(by='ID') %>% dplyr::group_by(variable) %>% dplyr::summarise(Total=sum(value)) %>% dplyr::mutate(Group = group[variable,'Stage']) %>% dplyr::filter(Group != 'NA') %>% dplyr::mutate(Tax='B. bifidum') %>% as.data.frame(),
  read.table("00.rawdata/CNV/100396_Dorea_A longicatena.CNV.txt",sep = "\t",header = T,row.names = 1) %>% as.data.frame() %>% tibble::rownames_to_column('ID') %>% reshape2::melt(by='ID') %>% dplyr::group_by(variable) %>% dplyr::summarise(Total=sum(value)) %>% dplyr::mutate(Group = group[variable,'Stage']) %>% dplyr::filter(Group != 'NA') %>% dplyr::mutate(Tax='D. longicatena') %>% as.data.frame(),
  read.table("00.rawdata/CNV/100666_Bacteroides cellulosilyticus.CNV.txt",sep = "\t",header = T,row.names = 1) %>% as.data.frame() %>% tibble::rownames_to_column('ID') %>% reshape2::melt(by='ID') %>% dplyr::group_by(variable) %>% dplyr::summarise(Total=sum(value)) %>% dplyr::mutate(Group = group[variable,'Stage']) %>% dplyr::filter(Group != 'NA') %>% dplyr::mutate(Tax='B. cellulosilyticus') %>% as.data.frame(),
  read.table("00.rawdata/CNV/103102_Bacteroides nordii.CNV.txt",sep = "\t",header = T,row.names = 1) %>% as.data.frame() %>% tibble::rownames_to_column('ID') %>% reshape2::melt(by='ID') %>% dplyr::group_by(variable) %>% dplyr::summarise(Total=sum(value)) %>% dplyr::mutate(Group = group[variable,'Stage']) %>% dplyr::filter(Group != 'NA') %>% dplyr::mutate(Tax='B. nordii') %>% as.data.frame(),
  read.table("00.rawdata/CNV/139103_Lachnospira eligens.CNV.txt",sep = "\t",header = T,row.names = 1) %>% as.data.frame() %>% tibble::rownames_to_column('ID') %>% reshape2::melt(by='ID') %>% dplyr::group_by(variable) %>% dplyr::summarise(Total=sum(value)) %>% dplyr::mutate(Group = group[variable,'Stage']) %>% dplyr::filter(Group != 'NA') %>% dplyr::mutate(Tax='L. eligens') %>% as.data.frame() 
) %>% ggplot(aes(Group,Total,fill=Group)) + facet_wrap('Tax',nrow = 1,scales = 'free') + geom_boxplot(outlier.shape = NA) + geom_point(position = 'jitter',size=0.5) + geom_signif(comparisons = analysis.list.info$compaired$Stage, color="black",step_increase = 0.08, map_signif_level=function(p) if(p < 0.01){sprintf("%.3f",p)}else{if(p < 0.05){sprintf("%.3f",p)}else{sprintf("%.3f",p)}}, test = wilcox.test) + theme_classic() + xlab('') + ylab('CNV of gmGUS') + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none") + scale_fill_manual(values = analysis.list.info$colors$Stage)

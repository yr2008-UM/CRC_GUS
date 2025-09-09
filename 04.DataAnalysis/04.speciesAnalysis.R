# ######################################################################
#
## Species-Level Analysis
## Author: Junru Chen, junru.chen2019@hotmail.com
## Version 1.0
## Date: 09.09.2025
#
## Description:
## Analyzes species-level cumulative GUS abundance/number and copy number 
## variation (CNV), performs differential abundance testing.
## Generates Figure 2f-h. Key components include:
## 1. Calculation of cumulative GUS abundance/number
## 2. Differential analysis
## 3. Visualization
## 4. Load CNV data from outputs of MIDAS2 and perform visualization
#
## Inputs:
##   - 00.rawdata/GUSsStat.csv: gmGUS metadata
##   - 00.rawdata/GUSabun_TPM.csv: TPM-normalized abundance matrix
##   - 00.rawdata/group.csv: Sample metadata
##   - 00.rawdata/CNV/: CNV data
#
## Outputs:
##   - Heatmaps of species abundance (Fig. 2f,h)
##   - Differential abundance plots (Fig. 2g)
##   - Boxplot for copy number variation (Fig. S4)
#
## Dependencies: 
## dplyr, ggplot2, ggpubr, reshape2
#
## Usage:
## 1. Set the working directory using setwd("the directory of scripts")
## 2. Ensure all input files are in the specified paths
## 3. Ensure all required packages were installed
## 4. Run script sections sequentially in RStudio or R command line
## 5. Output visualizations will be generated in R/RStudio
#
## Tested on: R 4.3.2 (macOS Ventura 13.5)
#
##############################################################################

# Install required packages if missing
required_pkgs <- c("dplyr", "ggplot2", "ggpubr", "reshape2")
install_missing <- required_pkgs[!required_pkgs %in% installed.packages()]
if(length(install_missing)) install.packages(install_missing)

# Load libraries
library(dplyr)     # Data manipulation
library(ggplot2)   # Plotting
library(ggpubr)    # Plotting
library(reshape2)  # Data reshaping

### SECTION 1: definition of function ###

## Function: group_mean
## Purpose: Calculate group-wise statistics (median, mean, SD) for a data matrix
## Inputs:
##   - group.subgroup: Dataframe with group assignments
##   - data: Data matrix (features x samples)
##   - colInfo: Column name containing group information
## Output: Dataframe with group statistics
group_mean <- function(group.subgroup,data,colInfo="Group"){
  # Get unique groups
  tmpg=unique(group.subgroup[,colInfo])
  tmpn=length(tmpg)
  rows=nrow(data)
  
  # Initialize result matrix
  result <- matrix(nrow=rows,ncol=tmpn*3)
  name <- NULL
  
  # Create column names
  for(i in tmpg){ name <- c(name,paste(i,'Median',sep='_'),paste(i,'Mean',sep='_'),paste(i,'SD',sep='_'))    }
  colnames(result)=name
  
  # Calculate statistics for each group
  for(i in tmpg){
    # Subset samples for current group
    tmpdata=data[,rownames(group.subgroup[which(group.subgroup[,colInfo]==i),,drop=F])]

    # colnames
    name <- c(paste(i,'Median',sep='_'),paste(i,'Mean',sep='_'),paste(i,'SD',sep='_'))
    
    # Calculate statistics for each feature
    for (j in 1:rows){
      tmp.median=median(as.numeric(tmpdata[j,]))
      tmp.mean=mean(as.numeric(tmpdata[j,]))
      tmp.SD=sd(as.numeric(tmpdata[j,]))
      result[j,name] <- c(tmp.median,tmp.mean,tmp.SD)
    }
  }
  
  # Format and return results
  rownames(result) <- rownames(data)
  result <- data.frame(result)
  return(result)
}

## Function: wilcoxon.FDR.TAX
## Purpose: Perform differntial test based on wilcoxon rank sum test and FDR correction
## Inputs:
##   - group.subgroup: Dataframe with group assignments
##   - data: Data matrix (features x samples)
## Output: List with statistics
wilcoxon.FDR.TAX <- function(data,group){
  
  data=data[,rownames(group)]
  tmpg=unique(group$Stage)
  species.abun <- data
  
  #statistics for groups
  result <- group_mean(group,data,'Stage')
  
  #wilcoxon + FDR correction
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
  
  # cutoff for sign species
  cutoff <- 61.6 #10% of samples
  cutFC <- log2(2)

  # Format and exact sign species
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
  
  # Combine results and format output
  allSign <- unique(c(signTaxG0,signTaxG1,signTaxG2,signTaxG3))
  signFeatureStageFC <- data.frame(Tax=allSign,MP=signTaxMP[allSign,'log2FC'],S0=signTaxS0[allSign,'log2FC'],S12=signTaxS12[allSign,'log2FC'],S34=signTaxS34[allSign,'log2FC']) %>% tibble::column_to_rownames('Tax')
  signFeatureStageMean <- result[allSign,c("Healthy_Mean","MP_Mean","S0_Mean",'SI_II_Mean',"SIII_IV_Mean")]
  
  return(list(testResult=result,signFC=signFeatureStageFC,signMean=signFeatureStageMean,signSpecies=signSpecies))
}

## Function: wilcoxon.FDR.mOTU4
## Purpose: Perform differntial test from mOTU4/metaphlan
## Inputs:
##   - group.subgroup: Dataframe with group assignments
##   - data: Data matrix (features x samples)
## Output: Dataframe with statistics
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


### SECTION 2: data loading ###

# Load metadata and gmGUS abundance data
GUSsStat <- read.csv("00.rawdata/GUSsStat.csv",header = T,row.names = 1)
GUSabun_TPM <- read.csv("00.rawdata/GUSabun_TPM.csv",header = T,row.names = 1)
group <- read.csv("00.rawdata/group.csv",header = T,row.names = 1)
groupV2 <- group %>% dplyr::mutate(ID=Subject_ID) %>% tibble::remove_rownames() %>% tibble::column_to_rownames('ID')

# Analysis configuration
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


### SECTION 3: analysis for species-level cumulative abundance ###

# Select species with >1 gmGUS
sel.taxes <- (table(GUSsStat[,'species']) %>% as.data.frame() %>% dplyr::filter(Freq > 1) %>% dplyr::arrange(-Freq))$Var1 %>% as.vector()

# Calculate cumulative abundance per species
species.abun <- matrix(nrow=length(sel.taxes),ncol = ncol(GUSabun_TPM))
rownames(species.abun) <-sel.taxes
colnames(species.abun) <- colnames(GUSabun_TPM)
for (mid.sp in sel.taxes) {
  species.abun[mid.sp,] <- GUSabun_TPM[GUSsStat %>% dplyr::filter(species == mid.sp) %>% rownames(),] %>% colSums() %>% as.vector()
}

# differential abundance analysis ###
wilcoxonResult <- wilcoxon.FDR.TAX(species.abun,group)

## visualization for abundance, fig. 2f ###

# Prepare data for heatmap (Z-scores)
heat.in <- wilcoxonResult$signMean
colnames(heat.in) <- c("HC","MP","S0","S12","S34")
heat.in <- t(scale(t(heat.in),center=T,scale=T))

# Create heatmap
heat.in %>% as.data.frame() %>% tibble::rownames_to_column('Tax') %>% reshape2::melt(by = 'Tax') %>% dplyr::mutate(Zscore=value) %>% ggplot(aes(variable,Tax,fill=Zscore)) + geom_tile() + xlab("")  + theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "top") + ylab("") + scale_fill_gradient2(low='navy',mid='white',high='firebrick3') + ylim(rownames(wilcoxonResult$signMean))

# Create plot for log2FC
wilcoxonResult$signFC %>% tibble::rownames_to_column('Tax') %>% reshape2::melt(by='Tax') %>% ggplot(aes(variable,Tax,fill=value)) + geom_tile() + xlab("") + theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "top",axis.text.y = element_blank()) + ylab("") + scale_fill_gradient2(low='navy',mid='white',high='firebrick3')  + ylim(rownames(wilcoxonResult$signMean))

# Create plot for up or down
rbind(
  wilcoxonResult$signSpecies$Up %>% unlist() %>% as.data.frame() %>% dplyr::rename(Tax = 1) %>% tibble::rownames_to_column('Group') %>% dplyr::mutate(Group=gsub("S34.*","S34",Group),Type='Up'),
  wilcoxonResult$signSpecies$Down %>% unlist() %>% as.data.frame() %>% dplyr::rename(Tax = 1) %>% tibble::rownames_to_column('Group') %>% dplyr::mutate(Group=gsub("S12.*","S12",Group),Type='Down')
) %>% ggplot(aes(Group,Tax,shape=Type,fill=Type)) + geom_point(size=2) + scale_shape_manual(values = c("Down"=25,"Up"=24)) + scale_fill_manual(values = c("Down"="blue","Up"="red")) + xlab("") + ylab("") + ylim(rownames(wilcoxonResult$signMean)) + theme_classic() + theme(legend.position = "top",axis.text.x = element_text(angle = 45,hjust = 1))

### SECTION 4: analysis for species-level cumulative number, Fig. 2g ###

# Prepare data
sel.taxes <- (table(GUSsStat[,'species']) %>% as.data.frame() %>% dplyr::filter(Freq > 1) %>% dplyr::arrange(-Freq))$Var1 %>% as.vector()
data <- GUSabun_TPM
data[data > 0] <- 1
species.abun <- matrix(nrow=length(sel.taxes),ncol = ncol(GUSabun_TPM))
rownames(species.abun) <-sel.taxes
colnames(species.abun) <- colnames(GUSabun_TPM)
for (mid.sp in sel.taxes) {
  species.abun[mid.sp,] <- data[GUSsStat %>% dplyr::filter(species == mid.sp) %>% rownames(),] %>% colSums() %>% as.vector()
}

# visualization
species.abun[c("Bacteroides cellulosilyticus","Bacteroides faecium","Bacteroides nordii"),] %>% as.data.frame() %>% tibble::rownames_to_column('ID') %>% reshape2::melt(by='ID') %>% dplyr::mutate(Group=group[variable,'Stage']) %>% ggplot(aes(Group,value,fill=Group)) + geom_violin() + geom_boxplot(col='black',size=1,width =0.1,outlier.colour = NA) + xlab("")  + theme_pubclean() + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "top") + ylab("") + facet_wrap('ID',nrow = 1)  + geom_signif(comparisons = analysis.list.info$compaired$Stage, color="black",step_increase = 0.08, map_signif_level=function(p) if(p < 0.01){sprintf("%.3f",p)}else{if(p < 0.05){sprintf("%.3f",p)}else{sprintf("%.3f",p)}}, test = wilcox.test) + scale_fill_manual(values=analysis.list.info$colors$Stage)


### SECTION 4: analysis for species-level relative abundance, Fig. 2h ###

## Analysis based on silva LTP

# Prepare data
selSpecies <- rownames(wilcoxonResult$signMean)
selSpecies <- intersect(
  selSpecies,
  read.table("00.rawdata/supp/Supp.Silva.txt",sep = "\t",header = T,row.names = 1,quote = "") %>% rownames()
)
selSpecies <- c(selSpecies,"Eubacterium eligens")
speciesAbun <- read.table("00.rawdata/supp/Supp.Silva.txt",sep = "\t",header = T,row.names = 1,quote = "")[selSpecies,]
speciesAbun <- speciesAbun[,paste0("X",group[intersect(groupV2[gsub("^X","",colnames(speciesAbun)),"Sample_ID"],rownames(group)),"Subject_ID"])]
colnames(speciesAbun) <- groupV2[gsub("^X","",colnames(speciesAbun)),"Sample_ID"]
sample1 <- group %>% dplyr::filter(Stage == "Healthy") %>% rownames()

# Mean abundance and differential analysis for CRC stages vs. healthy
meanH <- apply(speciesAbun,1,function(x){ if(sum(x) == 0){test <- 0;}else{test <-mean(as.numeric(x[sample1]));}}) %>% as.numeric
midData <- NA
for (midGroup in c("MP","S0","SI_II","SIII_IV")) {
  sample2 <- group %>% dplyr::filter(Stage == midGroup) %>% rownames()
  test <- apply(speciesAbun,1,function(x){ if(sum(x) == 0){p <- NA;}else{test <- wilcox.test(as.numeric(x[sample1]),as.numeric(x[sample2]),conf.int = T);p <- test$p.value;}}) %>% as.numeric
  testG <- apply(speciesAbun,1,function(x){ if(sum(x) == 0){p <- NA;}else{test <- wilcox.test(as.numeric(x[sample1]),as.numeric(x[sample2]),conf.int = T,alternative="greater");p <- test$p.value;}}) %>% as.numeric
  testL <- apply(speciesAbun,1,function(x){ if(sum(x) == 0){p <- NA;}else{test <- wilcox.test(as.numeric(x[sample1]),as.numeric(x[sample2]),conf.int = T,alternative="less");p <- test$p.value;}}) %>% as.numeric
  meanC <- apply(speciesAbun,1,function(x){ if(sum(x) == 0){test <- 0;}else{test <-mean(as.numeric(x[sample2]));}}) %>% as.numeric
  midData <- cbind(midData,data.frame(row.names=rownames(speciesAbun),Pvalue=test,greaterP=testL,lessP=testG,meanC=meanC))
}
colnames(midData) <- c('Test',paste0(c("Pvalue","greaterP","lessP","mean"),c(rep("_MP",4),rep("_S0",4),rep("_S12",4),rep("_S34",4))))
midData <- midData %>% dplyr::select(-Test)
midData$mean_HC <- meanH

# Visualization
orderS <- c("Bacteroides helcogenes","Bifidobacterium bifidum","Eubacterium eligens","Dorea longicatena","Bacteroides nordii","Bacteroides cellulosilyticus")
heat.in <- midData %>% dplyr::filter(Pvalue_MP < 0.05 | Pvalue_S0 < 0.05 | Pvalue_S12 < 0.05 | Pvalue_S34 < 0.05) %>% dplyr::mutate(HC=mean_HC,MP=mean_MP,S0=mean_S0,S12=mean_S12,S34=mean_S34) %>% dplyr::select(HC,MP,S0,S12,S34)
heat.in <- t(scale(t(heat.in),center=T,scale=T))
p1 <- heat.in %>% as.data.frame() %>% tibble::rownames_to_column('Tax') %>% reshape2::melt(by = 'Tax') %>% dplyr::mutate(Zscore=value) %>% ggplot(aes(variable,Tax,fill=Zscore)) + geom_tile() + xlab("")  + theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "top") + ylab("") + scale_fill_gradient2(low='navy',mid='white',high='firebrick3') + ylim(orderS)
p2 <- midData %>% dplyr::filter(Pvalue_MP < 0.05 | Pvalue_S0 < 0.05 | Pvalue_S12 < 0.05 | Pvalue_S34 < 0.05) %>% dplyr::mutate(HC=mean_HC,MP=mean_MP,S0=mean_S0,S12=mean_S12,S34=mean_S34) %>% dplyr::select(HC,MP,S0,S12,S34) %>% dplyr::mutate(MP=ifelse(MP > HC,log2(MP/HC),-log2(HC/MP)),S0=ifelse(S0 > HC,log2(S0/HC),-log2(HC/S0)),S12=ifelse(S12 > HC,log2(S12/HC),-log2(HC/S12)),S34=ifelse(S34 > HC,log2(S34/HC),-log2(HC/S34))) %>% dplyr::select(-HC) %>% tibble::rownames_to_column('Tax') %>% reshape2::melt(by='Tax') %>% ggplot(aes(variable,Tax,fill=value)) + geom_tile() + xlab("") + theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "top",axis.text.y = element_blank()) + ylab("") + scale_fill_gradient2(low='navy',mid='white',high='firebrick3') + ylim(orderS)
p3 <- midData %>% dplyr::filter(Pvalue_MP < 0.05 | Pvalue_S0 < 0.05 | Pvalue_S12 < 0.05 | Pvalue_S34 < 0.05) %>% dplyr::mutate(MP=ifelse(Pvalue_MP < 0.05,ifelse(mean_MP > mean_HC,'Up','Down'),'NA'),S0=ifelse(Pvalue_S0 < 0.05,ifelse(mean_S0 > mean_HC,'Up','Down'),'NA'),S12=ifelse(Pvalue_S12 < 0.05,ifelse(mean_S12 > mean_HC,'Up','Down'),'NA'),S34=ifelse(Pvalue_S34 < 0.05,ifelse(mean_S34 > mean_HC,'Up','Down'),'NA')) %>% dplyr::select(MP,S0,S12,S34) %>% tibble::rownames_to_column('Tax') %>% as.data.frame() %>% reshape2::melt("Tax") %>% dplyr::filter(value != 'NA') %>% ggplot(aes(variable,Tax,shape=value,fill=value)) + geom_point(size=2) + scale_shape_manual(values = c("Down"=25,"Up"=24)) + scale_fill_manual(values = c("Down"="blue","Up"="red")) + xlab("") + ylab("") + ylim(orderS) + theme_classic() + theme(legend.position = "top",axis.text.y = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))
ggarrange(p1,p2,p3,nrow = 1)

## analysis based on mOTU4
dbInfo <- read.table("00.rawdata/mOTUs4/mOTUsv4.0.gtdb.taxonomy.80mv.tsv",sep = "\t",header = T,row.names = 2)
data <- read.table("00.rawdata/mOTUs4/merged.mOTUs.table.txt",sep = "\t",header = T,row.names = 1)
data$Species <- as.vector(dbInfo[rownames(data),'species'])
data <- data %>% as.data.frame() %>% reshape2::melt('Species') %>% dplyr::group_by(Species,variable) %>% dplyr::summarise(meanS=sum(value)) %>% as.data.frame() %>% reshape2::dcast(Species ~ variable)
data <- data %>% as.data.frame() %>% dplyr::filter(Species != '') %>% tibble::column_to_rownames('Species')
selSpecies <- rownames(data)[grep("Bacteroides.faecium|helcogenes|nordii|Subdoligranulum.variabile|bifidum|cellulosilyticus|muciniphila|longicatena|eligens",rownames(data))]
speciesAbun <- data[selSpecies,rownames(group)]
testResult <- wilcoxon.FDR.mOTU4(speciesAbun,group)

## analysis based on metaphlan
data <- read.table("00.rawdata/supp/Supp.metaphlan.txt",sep = "\t",header = T,row.names = 1)
selSpecies <- rownames(data)[grep("Bacteroides.faecium|helcogenes|nordii|Subdoligranulum.variabile|bifidum|cellulosilyticus|muciniphila|longicatena|eligens",rownames(data))]
groupV2 <- group %>% dplyr::mutate(ID=Subject_ID) %>% tibble::remove_rownames() %>% tibble::column_to_rownames('ID')
speciesAbun <- data[selSpecies,]
speciesAbun <- speciesAbun[,paste0("X",group[intersect(groupV2[gsub("^X","",colnames(speciesAbun)),"Sample_ID"],rownames(group)),"Subject_ID"])]
colnames(speciesAbun) <- groupV2[gsub("^X","",colnames(speciesAbun)),"Sample_ID"]
speciesAbun <- speciesAbun[,rownames(group)]
testResult <- wilcoxon.FDR.mOTU4(speciesAbun,group)


### SECTION 5: analysis for species-level CNV, Fig. S4 ###

# Prepare data
plotIN <- rbind(
  read.table("00.rawdata/CNV/100129_Akkermansia muciniphila.CNV.txt",sep = "\t",header = T,row.names = 1) %>% as.data.frame() %>% tibble::rownames_to_column('ID') %>% reshape2::melt(by='ID') %>% dplyr::group_by(variable) %>% dplyr::summarise(Total=sum(value)) %>% dplyr::mutate(Group = group[variable,'Stage']) %>% dplyr::filter(Group != 'NA') %>% dplyr::mutate(Tax='A. muciniphila') %>% as.data.frame(),
  read.table("00.rawdata/CNV/100145_Bifidobacterium bifidum.CNV.txt",sep = "\t",header = T,row.names = 1) %>% as.data.frame() %>% tibble::rownames_to_column('ID') %>% reshape2::melt(by='ID') %>% dplyr::group_by(variable) %>% dplyr::summarise(Total=sum(value)) %>% dplyr::mutate(Group = group[variable,'Stage']) %>% dplyr::filter(Group != 'NA') %>% dplyr::mutate(Tax='B. bifidum') %>% as.data.frame(),
  read.table("00.rawdata/CNV/100396_Dorea_A longicatena.CNV.txt",sep = "\t",header = T,row.names = 1) %>% as.data.frame() %>% tibble::rownames_to_column('ID') %>% reshape2::melt(by='ID') %>% dplyr::group_by(variable) %>% dplyr::summarise(Total=sum(value)) %>% dplyr::mutate(Group = group[variable,'Stage']) %>% dplyr::filter(Group != 'NA') %>% dplyr::mutate(Tax='D. longicatena') %>% as.data.frame(),
  read.table("00.rawdata/CNV/100666_Bacteroides cellulosilyticus.CNV.txt",sep = "\t",header = T,row.names = 1) %>% as.data.frame() %>% tibble::rownames_to_column('ID') %>% reshape2::melt(by='ID') %>% dplyr::group_by(variable) %>% dplyr::summarise(Total=sum(value)) %>% dplyr::mutate(Group = group[variable,'Stage']) %>% dplyr::filter(Group != 'NA') %>% dplyr::mutate(Tax='B. cellulosilyticus') %>% as.data.frame(),
  read.table("00.rawdata/CNV/103102_Bacteroides nordii.CNV.txt",sep = "\t",header = T,row.names = 1) %>% as.data.frame() %>% tibble::rownames_to_column('ID') %>% reshape2::melt(by='ID') %>% dplyr::group_by(variable) %>% dplyr::summarise(Total=sum(value)) %>% dplyr::mutate(Group = group[variable,'Stage']) %>% dplyr::filter(Group != 'NA') %>% dplyr::mutate(Tax='B. nordii') %>% as.data.frame(),
  read.table("00.rawdata/CNV/139103_Lachnospira eligens.CNV.txt",sep = "\t",header = T,row.names = 1) %>% as.data.frame() %>% tibble::rownames_to_column('ID') %>% reshape2::melt(by='ID') %>% dplyr::group_by(variable) %>% dplyr::summarise(Total=sum(value)) %>% dplyr::mutate(Group = group[variable,'Stage']) %>% dplyr::filter(Group != 'NA') %>% dplyr::mutate(Tax='L. eligens') %>% as.data.frame() 
)

# visualization
plotIN %>% ggplot(aes(Group,Total,fill=Group)) + facet_wrap('Tax',nrow = 1,scales = 'free') + geom_boxplot(outlier.shape = NA) + geom_point(position = 'jitter',size=0.5) + geom_signif(comparisons = analysis.list.info$compaired$Stage, color="black",step_increase = 0.08, map_signif_level=function(p) sprintf("%.3f",p), test = wilcox.test) + theme_classic() + xlab('') + ylab('CNV of gmGUS') + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none") + scale_fill_manual(values = analysis.list.info$colors$Stage)

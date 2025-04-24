# ######################################################################
#
## GUS analysis for Figure 3AB
## Author: JR C
## Date: 01.01.2025
#
# ######################################################################
library(dplyr)
library(ggplot2)
library(ggpubr)

##### Get Data #####
loop.colors <- c("No Loop"="#FABA39","Mini-Loop 1"="#1AE4B6","Loop 1"="#4686FB","Loop 2"="#7A0403", "Mini-Loop 2"=  "#E4460A", "Mini-Loop 1,2"="#A2FC3C","No coverage" ="#30123B")
GUSabun_TPM <- read.csv("00.rawdata/GUSabun_TPM.csv",header = T,row.names = 1)
group <- read.csv("00.rawdata/group.csv",header = T,row.names = 1)
GUSsStat <- read.csv("00.rawdata/GUSsStat.csv",header = T,row.names = 1)
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

##### diff of gmGUSs between groups #####
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
  cutoff <- nrow(group)*0.1
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
signGUSs <- wilcoxon.FDR.TAX(GUSabun_TPM,group)

##### venn analysis #####
library(VennDiagram)
vennPlot <- venn.diagram(list(
  Healthy.VS.MP=c(signGUSs$signSpecies$Up$MP,signGUSs$signSpecies$Down$MP),
  Healthy.VS.S0=c(signGUSs$signSpecies$Up$S0,signGUSs$signSpecies$Down$S0),
  Healthy.VS.SI_II=c(signGUSs$signSpecies$Up$S12,signGUSs$signSpecies$Down$S12),
  Healthy.VS.SIII_IV=c(signGUSs$signSpecies$Up$S34,signGUSs$signSpecies$Down$S34)
  ),
  filename = NULL,
  fill = c('blue','yellow', 'purple', 'green')
)
grid.draw(vennPlot)

##### figure 3B #####
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

# tree
library(ggtree)
library(treeio)
beast_tree <- read.newick("00.rawdata/38gmGUS.tre")
ggtree(beast_tree) + geom_tiplab(size=5)

# loop type
GUSsStat[ggtreeOrderGUS,] %>% dplyr::select(Loop) %>% tibble::rownames_to_column(var='GUSID') %>% dplyr::mutate(Lab=gsub("Loop ","L",gsub("Mini-Loop ","mL",gsub("No Loop","NL",Loop)))) %>% ggplot(aes(Loop,GUSID,fill=Loop,label=Lab)) + geom_tile()  + geom_text() + scale_fill_manual(values = loop.colors) + ylim(rev(ggtreeOrderGUS)) + theme_pubclean() + theme(legend.position = "null") + xlab("")

# heatmap of mean abundance
library(pheatmap)
pheatmap((signGUSs$testResult %>% dplyr::select(Healthy_Mean,MP_Mean,S0_Mean,SI_II_Mean,SIII_IV_Mean))[ggtreeOrderGUS,],scale="row",cluster_cols = F,cluster_rows = F,color = colorRampPalette(colors = c("navy","white","firebrick3"))(100),cellheight = 12,cellwidth = 15,filename = NA)

# up and down
signFeatureStage <- rbind(
  signGUSs$signSpecies$Up %>% unlist() %>% as.data.frame() %>% dplyr::rename(Tax = 1) %>% tibble::rownames_to_column('Group') %>% dplyr::mutate(Group=gsub("S34.*","S34",Group),Type='Up') %>% dplyr::mutate(Group=gsub("S12.*","S12",Group)),
  signGUSs$signSpecies$Down %>% unlist() %>% as.data.frame() %>% dplyr::rename(Tax = 1) %>% tibble::rownames_to_column('Group') %>% dplyr::mutate(Group=gsub("S12.*","S12",Group),Type='Down')  %>% dplyr::mutate(Group=gsub("S34.*","S34",Group))
) %>% dplyr::mutate(Group=gsub("S0.*","S0",gsub("MP.*","MP",Group)))
signFeatureStage %>% ggplot(aes(Group,Tax,shape=Type,fill=Type)) + geom_point(size=2) + scale_shape_manual(values = c("Down"=25,"Up"=24)) + scale_fill_manual(values = c("Down"="blue","Up"="red")) + xlab("") + ylab("") + ylim(rev(ggtreeOrderGUS)) + theme_classic() + theme(legend.position = "top",axis.text.x = element_text(angle = 45,hjust = 1))

# prevalence
calculatePrevalence <- function(inputData,groupCol){
  group <- group %>% dplyr::filter(! is.na(get(groupCol)))
  mid <- inputData
  mid[mid > 0] <- 1
  mark <- matrix(nrow = nrow(inputData),ncol = length(unique(group[,groupCol])))
  colnames(mark) <- unique(group[,groupCol])
  rownames(mark) <- rownames(inputData)
  for (i in 1:nrow(inputData)) {
    for (mid.group in unique(group[,groupCol])) {
      mid.group2 <- group %>% filter(get(groupCol) == mid.group) %>% rownames()
      mid.mark <- (sum( mid[i,mid.group2]) / length(mid.group2) ) * 100
      mark[i,mid.group] <- mid.mark
    }
  }
  return(mark)
}
preResult <- calculatePrevalence(GUSabun_TPM[ggtreeOrderGUS,],"Stage") %>% reshape2::melt(by="row.names") 
preResult$Var2 <- factor(preResult$Var2,levels = analysis.list.info$levels$Stage)
preResult %>% ggplot(aes(value,Var1,color=Var2)) + geom_vline(xintercept = c(10),linetype=2,color="grey") + geom_line(group="Var2") + geom_point() + facet_wrap("Var2",nrow = 1,ncol = 5) + theme_classic() + xlab("Prevalence (100%)") + ylab("") + scale_color_manual(values = analysis.list.info$colors$Stage)  + theme(legend.position = "none") + ylim(rev(ggtreeOrderGUS))

##### analysis for AUS, FRA, GER #####
cohortResult <- list()
for (midCountry in c('AUS','FRA','GER')) {
  groupMuti <- read.csv("00.rawdata/Cohorts/Cohorts.group.csv",row.names = 1,header = T)  %>% dplyr::mutate(Stage = ifelse(Group == 'CTR','CTR',ifelse(Stage %in% c('S3','S4'),'S34',Stage))) %>% dplyr::filter(Country == midCountry) %>% dplyr::filter(Stage != 'NA') 
  speciesAbun <- read.table("00.rawdata/Cohorts/Cohorts.BGUS.abun.TPM.txt",sep = "\t",header = T,row.names = 1,quote = "",check.names = F)[ggtreeOrderGUS,rownames(groupMuti)]
  sample1 <- groupMuti %>% dplyr::filter(Stage == "CTR") %>% rownames()
  meanH <- apply(speciesAbun,1,function(x){ if(sum(x) == 0){test <- 0;}else{test <-mean(as.numeric(x[sample1]));}}) %>% as.numeric
  midData <- NA
  for (midGroup in setdiff(unique(groupMuti$Stage),'CTR')) {
    sample2 <- groupMuti %>% dplyr::filter(Stage == midGroup) %>% rownames()
    test <- apply(speciesAbun,1,function(x){ if(sum(x) == 0){p <- NA;}else{test <- wilcox.test(as.numeric(x[sample1]),as.numeric(x[sample2]),conf.int = T);p <- test$p.value;}}) %>% as.numeric
    meanC <- apply(speciesAbun,1,function(x){ if(sum(x) == 0){test <- 0;}else{test <-mean(as.numeric(x[sample2]));}}) %>% as.numeric
    middf <- data.frame(row.names=rownames(speciesAbun),Pvalue=test,meanC=meanC)
    colnames(middf) <- paste0(c('Pvalue_','mean_'),midGroup)
    midData <- cbind(midData,middf)
  }
  midData <- midData %>% dplyr::select(-midData)
  midData$mean_HC <- meanH
  cohortResult[[midCountry]] <- midData
}

expectedGUSs <- unique(c(signGUSs$signSpecies$Up$S12,signGUSs$signSpecies$Up$S34,signGUSs$signSpecies$Down$S12,signGUSs$signSpecies$Down$S34))
p1 <- rbind(
  cohortResult$AUS[expectedGUSs,] %>% dplyr::mutate(Type=ifelse(Pvalue_S12 < 0.05 | Pvalue_S34  < 0.05,'sign','')) %>% dplyr::filter(Type == 'sign') %>% dplyr::mutate(S12=ifelse(Pvalue_S12 < 0.05,ifelse(mean_HC>mean_S12,'Down','Up'),''),S34=ifelse(Pvalue_S34 < 0.05,ifelse(mean_HC>mean_S34,'Down','Up'),'')) %>% dplyr::select(S12,S34) %>% tibble::rownames_to_column('gus') %>% reshape2::melt('gus') %>% dplyr::mutate(Country='AUS'),
  cohortResult$FRA[expectedGUSs,] %>% dplyr::mutate(Type=ifelse(Pvalue_S12 < 0.05 | Pvalue_S34  < 0.05,'sign','')) %>% dplyr::filter(Type == 'sign') %>% dplyr::mutate(S12=ifelse(Pvalue_S12 < 0.05,ifelse(mean_HC>mean_S12,'Down','Up'),''),S34=ifelse(Pvalue_S34 < 0.05,ifelse(mean_HC>mean_S34,'Down','Up'),'')) %>% dplyr::select(S12,S34) %>% tibble::rownames_to_column('gus') %>% reshape2::melt('gus') %>% dplyr::mutate(Country='FRA'),
  cohortResult$GER[expectedGUSs,] %>% dplyr::mutate(Type=ifelse(Pvalue_S12 < 0.05 | Pvalue_S34  < 0.05,'sign','')) %>% dplyr::filter(Type == 'sign') %>% dplyr::mutate(S12=ifelse(Pvalue_S12 < 0.05,ifelse(mean_HC>mean_S12,'Down','Up'),''),S34=ifelse(Pvalue_S34 < 0.05,ifelse(mean_HC>mean_S34,'Down','Up'),'')) %>% dplyr::select(S12,S34) %>% tibble::rownames_to_column('gus') %>% reshape2::melt('gus') %>% dplyr::mutate(Country='GER')
) %>% ggplot(aes(variable,gus,fill=value,shape=value)) + facet_wrap('Country',nrow = 1,ncol = 3) + geom_point() + scale_shape_manual(values = c("Down"=25,"Up"=24)) + scale_fill_manual(values = c("Down"="blue","Up"="red")) + ylim(rev(ggtreeOrderGUS)) + xlab("") + ylab("") + theme_bw()

p2 <- rbind(
  scale(cohortResult$AUS %>% dplyr::select(mean_HC,mean_S12,mean_S34) %>% t(),center = T,scale = T) %>% as.data.frame() %>% t() %>% as.data.frame() %>% tibble::rownames_to_column('gus') %>% reshape2::melt('gus') %>% dplyr::mutate(Country='AUS'),
  scale(cohortResult$FRA %>% dplyr::select(mean_HC,mean_S12,mean_S34) %>% t(),center = T,scale = T) %>% as.data.frame() %>% t() %>% as.data.frame() %>% tibble::rownames_to_column('gus') %>% reshape2::melt('gus') %>% dplyr::mutate(Country='FRA'),
  scale(cohortResult$GER %>% dplyr::select(mean_HC,mean_S12,mean_S34) %>% t(),center = T,scale = T) %>% as.data.frame() %>% t() %>% as.data.frame() %>% tibble::rownames_to_column('gus') %>% reshape2::melt('gus') %>% dplyr::mutate(Country='GER')
) %>% ggplot(aes(variable,gus,fill=value)) + facet_wrap('Country',nrow = 1,ncol = 3) + geom_tile() + scale_fill_gradient2(low = 'lightblue',mid = 'white',high='red') + ylim(rev(ggtreeOrderGUS)) + xlab("") + ylab("") + theme_bw()

ggarrange(p1,p2,nrow = 1,ncol = 2)

output <- cbind(
  cohortResult$AUS %>% dplyr::select(mean_HC,mean_S12,mean_S34,Pvalue_S12,Pvalue_S34),
  cohortResult$FRA %>% dplyr::select(mean_HC,mean_S12,mean_S34,Pvalue_S12,Pvalue_S34),
  cohortResult$GER %>% dplyr::select(mean_HC,mean_S12,mean_S34,Pvalue_S12,Pvalue_S34)
  )

##### analysis for gender differences #####
speciesAbun <- GUSabun_TPM[ggtreeOrderGUS,rownames(group)]
sample1 <- group %>% dplyr::filter(Stage %in% c('S0','SI_II','SIII_IV')) %>% dplyr::filter(Gender == "M") %>% rownames()
sample2 <- group %>% dplyr::filter(Stage %in% c('S0','SI_II','SIII_IV')) %>% dplyr::filter(Gender == "F") %>% rownames()
meanM <- apply(speciesAbun,1,function(x){ if(sum(x) == 0){test <- 0;}else{test <-mean(as.numeric(x[sample1]));}}) %>% as.numeric
meanF <- apply(speciesAbun,1,function(x){ if(sum(x) == 0){test <- 0;}else{test <-mean(as.numeric(x[sample2]));}}) %>% as.numeric
test <- apply(speciesAbun,1,function(x){ if(sum(x) == 0){p <- NA;}else{test <- wilcox.test(as.numeric(x[sample1]),as.numeric(x[sample2]),conf.int = T);p <- test$p.value;}}) %>% as.numeric
middf <- data.frame(gmGUS=rownames(speciesAbun),Pvalue=test,meanFemale=meanF,meanMale=meanM)
middf %>% dplyr::filter(Pvalue < 0.05)
# for healthy group, no significant differences were observed.

##### analysis for Age, BMI, Brinkman.Index, Alcohol differences #####
library(psych)
midGroup <- group[rownames(group %>% dplyr::filter(Stage %in% c('S0','SI_II','SIII_IV'))),c('Age','BMI',"Brinkman.Index", "Alcohol")]
speciesAbun <- GUSabun_TPM[ggtreeOrderGUS,rownames(midGroup)] %>% t() %>% as.data.frame()
corrRe <- corr.test(speciesAbun,midGroup,method="spearman")
colnames(corrRe$p) <- paste0(colnames(corrRe$p),'P')
outIN <- cbind(corrRe$r,corrRe$p) %>% as.data.frame() %>% dplyr::mutate(Age=paste0(sprintf("%.3f",Age)," (",sprintf("%.3f",AgeP),")"),BMI=paste0(sprintf("%.3f",BMI)," (",sprintf("%.3f",BMIP),")"),Brinkman.Index=paste0(sprintf("%.3f",Brinkman.Index)," (",sprintf("%.3f",Brinkman.IndexP),")"),Alcohol=paste0(sprintf("%.3f",Alcohol)," (",sprintf("%.3f",AlcoholP),")")) %>% dplyr::select(Age,BMI,Brinkman.Index,Alcohol)
corrRe$p %>% as.data.frame() %>% tibble::rownames_to_column('ID') %>% reshape2::melt(by='ID') %>% dplyr::filter(value < 0.05) %>% dplyr::mutate(R=corrRe$r[ID,])
# for healthy group, only signficant correlation between Lac_eligens.GUS2 and Age was observed, which was also detected in CRC patients.

##### blocked wilcoxon rank sum test after factor analysis #####
library(coin)
midGroup <- group %>% dplyr::filter(Stage %in% c('Healthy','MP'))
speciesAbun <- GUSabun_TPM[ggtreeOrderGUS,rownames(midGroup)] %>% t()
speciesAbun <- cbind(speciesAbun,midGroup)
speciesAbun$Stage <- factor(as.vector(speciesAbun$Stage),levels = c('Healthy','MP'))
speciesAbun$Gender <- factor(as.vector(speciesAbun$Gender),levels = c('F','M'))
summary(aov(Par_chongii.GUS1~Stage + Gender, data=speciesAbun))

midGroup <- group %>% dplyr::filter(Stage %in% c('Healthy','SIII_IV'))
speciesAbun <- GUSabun_TPM[ggtreeOrderGUS,rownames(midGroup)] %>% t()
speciesAbun <- cbind(speciesAbun,midGroup)
speciesAbun$Stage <- factor(as.vector(speciesAbun$Stage),levels = c('Healthy','SIII_IV'))
summary(aov(Lac_eligens.GUS2 ~ Stage+Gender+Age, speciesAbun))

midGroup <- group %>% dplyr::filter(Stage %in% c('Healthy','SI_II'))
speciesAbun <- GUSabun_TPM[ggtreeOrderGUS,rownames(midGroup)] %>% t()
speciesAbun <- cbind(speciesAbun,midGroup)
speciesAbun$Stage <- factor(as.vector(speciesAbun$Stage),levels = c('Healthy','SI_II'))
speciesAbun$Gender <- factor(as.vector(speciesAbun$Gender),levels = c('F','M'))
summary(aov(Stu_stutzeri.GUS~Stage + Alcohol, data=speciesAbun))

midGroup <- group %>% dplyr::filter(Stage %in% c('Healthy','S0'))
speciesAbun <- GUSabun_TPM[ggtreeOrderGUS,rownames(midGroup)] %>% t()
speciesAbun <- cbind(speciesAbun,midGroup)
speciesAbun$Stage <- factor(as.vector(speciesAbun$Stage),levels = c('Healthy','S0'))
speciesAbun$Gender <- factor(as.vector(speciesAbun$Gender),levels = c('F','M'))
summary(aov(Ech_strongylocentroti.GUS~Stage + Brinkman.Index, data=speciesAbun))

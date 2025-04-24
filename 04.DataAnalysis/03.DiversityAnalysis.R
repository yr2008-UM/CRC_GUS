# ######################################################################
#
## Diversity analysis (Figure 2DE)
## Author: JR C
## Date: 01.01.2025
#
# ######################################################################
library(vegan)
library(ape)
library(dplyr)
library(ggplot2)
library(ggpubr)

##### Get Data #####
group <- read.csv("00.rawdata/group.csv",header = T,row.names = 1)
GUSsStat <- read.csv("00.rawdata/GUSsStat.csv",header = T,row.names = 1)
GUSabun_TPM <- read.csv("00.rawdata/GUSabun_TPM.csv",header = T,row.names = 1)
loop.colors <- c("No Loop"="#FABA39","Mini-Loop 1"="#1AE4B6","Loop 1"="#4686FB","Loop 2"="#7A0403", "Mini-Loop 2"=  "#E4460A", "Mini-Loop 1,2"="#A2FC3C","No coverage" ="#30123B")
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

##### PCoA Analysis #####
dataTB <- vegdist(t(GUSabun_TPM), method="bray")
group.info <- group[,'Stage']
names(group.info) <- rownames(group)
#we dist  
adonis(dataTB~group.info,permutations = 999)
#PCOA
PCOA <- pcoa(dataTB, correction="none", rn=NULL)
result <-PCOA$values[,"Relative_eig"]
pco1 = as.numeric(sprintf("%.3f",result[1]))*100
pco2 = as.numeric(sprintf("%.3f",result[2]))*100
pc = as.data.frame(PCOA$vectors)
pc$names = rownames(pc)
xlab=paste("PCoA1 (",pco1,"%)",sep="")
ylab=paste("PCoA2 (",pco2,"%)",sep="")
pc$Group <- factor(as.vector(group.info),levels = analysis.list.info$levels$Stage)
ggplot(pc, aes(Axis.1,Axis.2, fill=Group,color=Group,shape=Group)) +
  labs(x=xlab,y=ylab) +
  geom_hline(yintercept=0,linetype=4,color="grey") +
  geom_vline(xintercept=0,linetype=4,color="grey") +
  geom_point(size=4,alpha=0.7) +
  scale_fill_manual(values = analysis.list.info$colors$Stage)+
  scale_color_manual(values = analysis.list.info$colors$Stage)+
  stat_ellipse(show.legend = F,level = 0.95)+
  theme(axis.text.x=element_text(colour = "black",angle=45,vjust=1,hjust=1,size = 7), axis.text.y=element_text(colour = "black",size = 7),panel.background = element_rect(fill="white",color="black",linetype=1,size=1),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text = element_text(size=7))

##### adonis adjusted by confounding factors #####
midGroup <- group %>% dplyr::filter(Stage %in% c("Healthy","MP"))
midGroup$Stage <- factor(as.vector(midGroup$Stage),levels = c('Healthy','MP'))
adonis2(t(GUSabun_TPM[,rownames(midGroup)])~Stage + Age + Gender + BMI + Brinkman.Index + Alcohol , data = midGroup,permutations = 999, method = "bray", by = "margin",na.action=na.omit)

midGroup <- group %>% dplyr::filter(Stage %in% c("Healthy","S0"))
midGroup$Stage <- factor(as.vector(midGroup$Stage),levels = c('Healthy','S0'))
adonis2(t(GUSabun_TPM[,rownames(midGroup)])~Stage + Age + Gender + BMI + Brinkman.Index + Alcohol , data = midGroup,permutations = 999, method = "bray", by = "margin",na.action=na.omit)

midGroup <- group %>% dplyr::filter(Stage %in% c("Healthy","SI_II"))
midGroup$Stage <- factor(as.vector(midGroup$Stage),levels = c('Healthy','SI_II'))
adonis2(t(GUSabun_TPM[,rownames(midGroup)])~Stage + Age + Gender + BMI + Brinkman.Index + Alcohol , data = midGroup,permutations = 999, method = "bray", by = "margin",na.action=na.omit)

midGroup <- group %>% dplyr::filter(Stage %in% c("Healthy","SIII_IV"))
midGroup$Stage <- factor(as.vector(midGroup$Stage),levels = c('Healthy','SIII_IV'))
adonis2(t(GUSabun_TPM[,rownames(midGroup)])~Stage + Age + Gender + BMI + Brinkman.Index + Alcohol , data = midGroup,permutations = 999, method = "bray", by = "margin",na.action=na.omit)

##### difference of total number of BGUS #####
data <- GUSabun_TPM
data[data > 0] <- 1
alpha.result <- colSums(data) %>% as.data.frame() %>% dplyr::rename(No.Of.BGUS=1)
alpha.result %>% tibble::rownames_to_column('ID') %>% dplyr::mutate(Group=group[ID,'Stage']) %>% dplyr::filter(Group != 'NA') %>% ggplot(aes(Group,No.Of.BGUS,fill=Group)) + geom_violin() + geom_boxplot(col='black',size=1,width =0.1,outlier.colour = NA) + theme_classic() + xlab('') + ylab('Number of gmGUS') + scale_fill_manual(values = analysis.list.info$colors$Stage) + theme(legend.position = "null") + geom_signif(comparisons = analysis.list.info$compaired$Stage, color="black",step_increase = 0.08, map_signif_level=function(p) if(p < 0.01){'+'}else{if(p < 0.05){'*'}else{sprintf("%.2f",p)}}, test = wilcox.test) + xlab("") + ylab("Number of gmGUS") + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none") + xlim(analysis.list.info$levels$Stage)

##### difference of total abundance of BGUS #####
data <- GUSabun_TPM
alpha.result <- colSums(data) %>% as.data.frame() %>% dplyr::rename(Abundance.Of.BGUS=1)
alpha.result <- cbind(alpha.result,group)
alpha.result %>% dplyr::select(Abundance.Of.BGUS,Stage) %>% ggplot(aes(Stage,Abundance.Of.BGUS,fill=Stage)) + geom_violin() + geom_boxplot(col='black',size=1,width =0.1,outlier.colour = NA) + theme_classic() + xlab('') + ylab('Abundance of gmGUS') + scale_fill_manual(values = analysis.list.info$colors$Stage) + theme(legend.position = "null")  + geom_signif(comparisons = analysis.list.info$compaired$Stage, color="black",step_increase = 0.08, map_signif_level=function(p) if(p < 0.01){'+'}else{if(p < 0.05){'*'}else{sprintf("%.2f",p)}}, test = wilcox.test) + xlab("") + ylab("Abundance of gmGUS") + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none") + xlim(analysis.list.info$levels$Stage)

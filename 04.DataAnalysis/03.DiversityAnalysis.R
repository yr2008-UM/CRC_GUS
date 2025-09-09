# ######################################################################
#
## Diversity Analysis
## Author: Junru Chen, junru.chen2019@hotmail.com
## Version 1.0
## Date: 09.09.2025
#
## Description:
##   Performs beta-diversity analysis (PCoA) and compares alpha diversity 
##   metrics across CRC stages. Generates Figure 2d-e panels.
#
## Inputs:
##   - 00.rawdata/group.csv: Sample metadata
##   - 00.rawdata/GUSsStat.csv: gmGUS metadata
##   - 00.rawdata/GUSabun_TPM.csv: TPM-normalized abundance matrix
#
## Outputs:
##   - PCoA plot (Fig 2d)
##   - Alpha diversity boxplots (Fig 2e)
#
## Dependencies: vegan, ape, dplyr, ggplot2, ggpubr
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
required_pkgs <- c("vegan", "ape", "dplyr", "ggplot2", "ggpubr")
install_missing <- required_pkgs[!required_pkgs %in% installed.packages()]
if(length(install_missing)) install.packages(install_missing)

# Load libraries
library(vegan)     # Ecological diversity analysis
library(ape)       # Phylogenetics and evolution
library(dplyr)     # Data manipulation
library(ggplot2)   # Plotting
library(ggpubr)    # Plotting

### SECTION 1: data loading ###

# Load metadata and gmGUS abundance data
group <- read.csv("00.rawdata/group.csv",header = T,row.names = 1)
GUSsStat <- read.csv("00.rawdata/GUSsStat.csv",header = T,row.names = 1)
GUSabun_TPM <- read.csv("00.rawdata/GUSabun_TPM.csv",header = T,row.names = 1)

# Analysis configuration
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

### SECTION 2: beta diversity analysis (Figure 2d) ###

# Calculate Bray-Curtis dissimilarity
dataTB <- vegdist(t(GUSabun_TPM), method="bray")

# PERMANOVA test
group.info <- group[,'Stage']
names(group.info) <- rownames(group)
adtest <- adonis(dataTB~group.info,permutations = 999)
print(adtest$aov.tab)

# PCoA analysis & Prepare plotting data
PCOA <- pcoa(dataTB, correction="none", rn=NULL)
result <-PCOA$values[,"Relative_eig"]
pco1 = as.numeric(sprintf("%.3f",result[1]))*100
pco2 = as.numeric(sprintf("%.3f",result[2]))*100
pc = as.data.frame(PCOA$vectors)
pc$names = rownames(pc)
xlab=paste("PCoA1 (",pco1,"%)",sep="")
ylab=paste("PCoA2 (",pco2,"%)",sep="")
pc$Group <- factor(as.vector(group.info),levels = analysis.list.info$levels$Stage)

# Create PCoA plot
ggplot(pc, aes(Axis.1,Axis.2, fill=Group,color=Group,shape=Group)) +
  labs(x=xlab,y=ylab) +
  geom_hline(yintercept=0,linetype=4,color="grey") +
  geom_vline(xintercept=0,linetype=4,color="grey") +
  geom_point(size=4,alpha=0.7) +
  scale_fill_manual(values = analysis.list.info$colors$Stage)+
  scale_color_manual(values = analysis.list.info$colors$Stage)+
  stat_ellipse(show.legend = F,level = 0.95)+
  theme(axis.text.x=element_text(colour = "black",angle=45,vjust=1,hjust=1,size = 7), axis.text.y=element_text(colour = "black",size = 7),panel.background = element_rect(fill="white",color="black",linetype=1,size=1),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text = element_text(size=7))

# PERMANOVA adjusted by confounding factors
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

### SECTION 3: alpha diversity analysis (Figure 2e) ###

# Calculate total abundance per sample
data <- GUSabun_TPM
alpha.result <- colSums(data) %>% as.data.frame() %>% dplyr::rename(Abundance.Of.BGUS=1)
alpha.result <- cbind(alpha.result,group)
alpha.result <- alpha.result %>% dplyr::select(Abundance.Of.BGUS,Stage) 

# visualization
alpha.result %>% ggplot(aes(Stage,Abundance.Of.BGUS,fill=Stage)) + geom_violin() + geom_boxplot(col='black',size=1,width =0.1,outlier.colour = NA) + theme_classic() + xlab('') + ylab('Abundance of gmGUS') + scale_fill_manual(values = analysis.list.info$colors$Stage) + theme(legend.position = "null")  + geom_signif(comparisons = analysis.list.info$compaired$Stage, color="black",step_increase = 0.08, map_signif_level=function(p) if(p < 0.01){'+'}else{if(p < 0.05){'*'}else{sprintf("%.2f",p)}}, test = wilcox.test) + xlab("") + ylab("Abundance of gmGUS") + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none") + xlim(analysis.list.info$levels$Stage)

# Calculate richness (number of gmGUSs per sample)
data <- GUSabun_TPM
data[data > 0] <- 1
alpha.result <- colSums(data) %>% as.data.frame() %>% dplyr::rename(No.Of.BGUS=1)
alpha.result <- alpha.result %>% tibble::rownames_to_column('ID') %>% dplyr::mutate(Group=group[ID,'Stage']) %>% dplyr::filter(Group != 'NA')

# visualization
alpha.result %>% ggplot(aes(Group,No.Of.BGUS,fill=Group)) + geom_violin() + geom_boxplot(col='black',size=1,width =0.1,outlier.colour = NA) + theme_classic() + xlab('') + ylab('Number of gmGUS') + scale_fill_manual(values = analysis.list.info$colors$Stage) + theme(legend.position = "null") + geom_signif(comparisons = analysis.list.info$compaired$Stage, color="black",step_increase = 0.08, map_signif_level=function(p) if(p < 0.01){'+'}else{if(p < 0.05){'*'}else{sprintf("%.2f",p)}}, test = wilcox.test) + xlab("") + ylab("Number of gmGUS") + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none") + xlim(analysis.list.info$levels$Stage)
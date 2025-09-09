##############################################################################
#
## Title: Loop Category Analysis
## Author: Junru Chen, junru.chen2019@hotmail.com
## Version 1.0
## Date: 09.09.2025
#
## Description:
##   Analyzes loop category distribution, abundance differences, and 
##   taxonomic composition. Generates Fig. 2a,b,c and Fig. S3.
#
## Inputs:
##   - 00.rawdata/group.csv: Sample metadata
##   - 00.rawdata/GUSsStat.csv: gmGUS metadata with loop annotations
##   - 00.rawdata/GUSabun_TPM.csv: TPM-normalized gmGUS abundance matrix
#
## Outputs:
##   - Length distribution plot (Fig 2a)
##   - Taxonomic pie charts (Fig 2b)
##   - Loop abundance violin plots (Fig 2c)
##   - Loop abundance and GUS number analysis (Supplementary Fig. 3)
#
## Dependencies: dplyr, ggplot2, ggpubr, viridis, reshape2
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
required_pkgs <- c("dplyr", "ggplot2", "ggpubr", "viridis", "reshape2")
install_missing <- required_pkgs[!required_pkgs %in% installed.packages()]
if(length(install_missing)) install.packages(install_missing)

# Load required libraries
library(dplyr)      # Data manipulation
library(ggplot2)    # Plotting
library(ggpubr)     # Plotting
library(viridis)    # Color scales
library(reshape2)   # Data reshaping

### SECTION 1: data loading and initialization ###

# Load metadata and abundance matrix
group <- read.csv("00.rawdata/group.csv",header = T,row.names = 1)
GUSsStat <- read.csv("00.rawdata/GUSsStat.csv",header = T,row.names = 1)
GUSabun_TPM <- read.csv("00.rawdata/GUSabun_TPM.csv",header = T,row.names = 1)

# Define color scheme for loop categories
loop.colors <- c("No Loop"="#FABA39","Mini-Loop 1"="#1AE4B6","Loop 1"="#4686FB","Loop 2"="#7A0403", "Mini-Loop 2"=  "#E4460A", "Mini-Loop 1,2"="#A2FC3C","No coverage" ="#30123B")

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

### SECTION 2: length ditribution (Figure 2a) ###
GUSsStat %>% tibble::rownames_to_column(var="ID") %>% ggplot(aes(reorder(ID,Length),weight=Length,fill=Loop)) + geom_bar(width = 1) + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "top") + xlab(paste0(nrow(data)," gmGUS Sequences Ordered by Length")) + ylab("Length (AA)")  + scale_fill_manual(values = loop.colors)
prop.table(table(GUSsStat$Loop)) * 100 

### SECTION 3: taxonomic pie charts (Figure 2b) ###
gg.result <- list()
i<- 1
for (mid.loop in unique(GUSsStat$Loop)) {
  # Filter data for current loop type
  mid.data <- GUSsStat %>% dplyr::filter(Loop == mid.loop) %>% dplyr::select(genus) %>% table %>% as.data.frame() %>% dplyr::arrange(-Freq)
  
  # Handle unclassified genera
  if(length(which(as.vector(mid.data$genus) == 'Unclassified')) > 0){
    mid.data <- rbind(
      mid.data[-which(as.vector(mid.data$genus) == 'Unclassified'),],
      mid.data %>% dplyr::filter(genus == 'Unclassified')
    )
  }
  
  # Prepare data for plotting
  mid.others <- sum(mid.data$Freq) - ((head(mid.data,n=5))$Freq  %>% sum())
  mid.data <- head(mid.data,n=5)
  mid.data$genus <- factor(mid.data$genus %>% as.vector(),levels = c(mid.data$genus %>% as.vector(),"Others") %>% rev())
  mid.data[6,] <- c("Others",mid.others)
  mid.col <- c(colorRampPalette(c(loop.colors[mid.loop],"white"))(18)[c(1,5,9,13,17)],"grey62") %>% rev
  mid.data$Freq <-mid.data$Freq %>% as.numeric()
  mid.label <- paste(mid.data$genus %>% as.vector(), paste('(',sprintf("%.1f",(mid.data$Freq / sum(mid.data$Freq)) * 100), '%)', sep = ''), sep = ' ') %>% rev
  
  # Create pie chart
  gg.result[[mid.loop]] <- ggplot(mid.data, mapping = aes(x = 'Content', y = Freq, fill = genus)) + geom_bar(stat = 'identity', position = 'stack', width = 1) + coord_polar(theta = 'y') + theme_void() + scale_fill_manual(labels = mid.label,values =mid.col) + labs(title = mid.loop)  #c( palette()[2:6],"grey62") %>% rev
  i <- i+1
}

# Arrange all pie charts in a grid
ggarrange(plotlist = gg.result, nrow = 3, ncol = 3)

### SECTION 4: analysis for loop abundance (Figure 2c, Supplementary Fig. 3a) ###

# Calculate abundance per loop category
loop.abun <- matrix(nrow=length(unique(GUSsStat$Loop)),ncol = ncol(GUSabun_TPM))
rownames(loop.abun) <- unique(GUSsStat$Loop)
colnames(loop.abun) <- colnames(GUSabun_TPM)
for (mid.loop in unique(GUSsStat$Loop)) {
  loop.abun[mid.loop,] <- GUSabun_TPM[GUSsStat %>% dplyr::filter(Loop == mid.loop) %>% rownames(),] %>% colSums() %>% as.vector()
}
loop.abun <- loop.abun %>% t() %>% as.data.frame()

# Kruskal-Wallis tests for group differences
loop.abunKW <- loop.abun
loop.abunKW$Group <- group[rownames(loop.abun),'Stage']
loop.abunKW <- loop.abunKW %>% dplyr::filter(Group != 'NA')
kruskal.test(`No Loop` ~ Group, data = loop.abunKW)
kruskal.test(`Loop 1` ~ Group, data = loop.abunKW) 
kruskal.test(`Loop 2` ~ Group, data = loop.abunKW) 
kruskal.test(`Mini-Loop 1` ~ Group, data = loop.abunKW) 
kruskal.test(`Mini-Loop 2` ~ Group, data = loop.abunKW) 
kruskal.test(`Mini-Loop 1,2` ~ Group, data = loop.abunKW) 

# Visualization
mid.data <- loop.abun
mid.data$Group <- group[rownames(mid.data),'Stage']
mid.data <- mid.data %>% dplyr::filter(Group != "NA" & ! is.na(Group))
mid.data <- reshape2::melt(mid.data,by="Group")
mid.data$Group <- factor(mid.data$Group,levels = analysis.list.info$levels$Stage)
mid.data$variable <- factor(as.vector(mid.data$variable),levels = names(loop.colors))
mid.data %>% ggplot(aes(Group,value,fill=Group)) + geom_violin() + geom_boxplot(col='black',size=1,width =0.1,outlier.colour = NA) + scale_fill_manual(values = analysis.list.info$colors$Stage)+ theme_classic() + facet_wrap("variable",nrow = 1,ncol = 7,scales = "free") + geom_signif(comparisons = analysis.list.info$compaired$Stage, color="black",step_increase = 0.08, map_signif_level=function(p) sprintf("%.3f",p), test = wilcox.test) + xlab("") + ylab("Abundance") + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")

### SECTION 5: analysis for GUS number of loops  (Supplementary Fig. 3b) ###

# prepare data for visualization
data <- GUSabun_TPM
data[data > 0] <- 1
result <- matrix(nrow = ncol(data),ncol = 7)
rownames(result) <- colnames(data)
colnames(result) <- names(loop.colors)[1:7]
result[is.na(result)] <- 0
for (i in colnames(data)) {
  for (j in names(loop.colors)[1:7]) {
    result[i,j] <- sum(data[GUSsStat %>% dplyr::filter(Loop == j) %>% rownames(),i])
  }
}
result <- as.data.frame(result)
result$Group <- as.vector(group[rownames(result),'Stage'])

#visualization
result %>% reshape2::melt("Group") %>% ggplot(aes(Group,value,fill=Group)) + geom_violin() + geom_boxplot(col='black',size=1,width =0.1,outlier.colour = NA) + scale_fill_manual(values = analysis.list.info$colors[['Stage']])+ theme_classic() + facet_wrap("variable",nrow = 1,ncol = 7,scales = "free") + geom_signif(comparisons = analysis.list.info$compaired[['Stage']], color="black",step_increase = 0.08, map_signif_level=function(p) sprintf("%.3f",p), test = wilcox.test) + xlab("") + ylab("Number of gmGUS") + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")

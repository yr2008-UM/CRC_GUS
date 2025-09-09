##############################################################################
#
## Rarefaction Curves and Plots for References and 550 gmGUSs
## Author: Junru Chen, junru.chen2019@hotmail.com
## Version 1.0
## Date: 09.09.2025
#
## Description:
##   Generates rarefaction curves for gmGUS diversity and bar plots for
##   reference sequences and identified gmGUSs. Corresponds to Figure 1
##   and Supplementary Figures 1,2 in the manuscript.
#
## Inputs:
##   - 00.rawdata/group.csv: Sample metadata with rownames as sample IDs
##   - 00.rawdata/GUSabun_ABS.csv: Absolute abundance matrix of gmGUSs (rows=gmGUSs, columns=samples)
##   - 00.rawdata/referenceStat.csv: Reference GUS sequences metadata
##   - 00.rawdata/GUSsStat.csv: Identified gmGUS sequences metadata
#
## Outputs:
##   - Rarefaction curves (Figure 1b)
##   - Bar plots for references and 550 gmGUSs (Supplementary Figures 1,2)
#
## Dependencies: vegan, dplyr, ggplot2
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
required_pkgs <- c("vegan", "dplyr", "ggplot2")
install_missing <- required_pkgs[!required_pkgs %in% installed.packages()]
if(length(install_missing)) install.packages(install_missing)

# Load libraries
library(vegan)    # For diversity analysis
library(dplyr)    # For data manipulation
library(ggplot2)  # For plotting

### SECTION 1: data loading ###
# Read sample metadata (row names = sample IDs)
group <- read.csv("00.rawdata/group.csv", header = TRUE, row.names = 1)

# Read absolute abundance matrix (rows = gmGUSs, columns = samples)
GUSabun_ABS <- read.csv("00.rawdata/GUSabun_ABS.csv", header = TRUE, row.names = 1)

# Read reference and gmGUS sequence metadata
referenceStat <- read.csv("00.rawdata/referenceStat.csv", header = TRUE, row.names = 1)
GUSsStat <- read.csv("00.rawdata/GUSsStat.csv", header = TRUE, row.names = 1)

### SECTION 2: rarefraction curves, Fig. 1b ###

# Generate rarefaction curve for all samples
data.spec <- specaccum(t(GUSabun_ABS), method = "random")  
plot(data.spec, ci.type = "poly", col = "blue", lwd = 2, ci.lty = 0,
     ci.col = "lightblue", xlab = "Number of samples",
     ylab = "Number of gmGUS")

# Add healthy group curve
data.spec <- specaccum(t(GUSabun_ABS[,group %>% dplyr::filter(Stage == "Healthy") %>% rownames()]), method = "random")  
plot(data.spec, ci.type = "poly", col = '#cecccb', lwd = 2, ci.lty = 0,
     ci.col = "lightblue", xlab = "Number of samples",
     ylab = "Number of gmGUS",add=T)

# Add CRC group curve
data.spec <- specaccum(t(GUSabun_ABS[,group %>% dplyr::filter(Stage %in% c('S0','SI_II','SIII_IV ')) %>% rownames()]), method = "random")  
plot(data.spec, ci.type = "poly", col = '#f3764a', lwd = 2, ci.lty = 0,
     ci.col = "lightblue", xlab = "Number of samples",
     ylab = "Number of gmGUS",add=T)

### SECTION 3: barplot for loops of references (Supplementary Figures 1a) ###
referenceStat %>% tibble::rownames_to_column(var="ID") %>% ggplot(aes(reorder(ID,Length),weight=Length,fill=Loop)) + geom_bar(width = 1,color="grey") + theme_classic() + theme(legend.position = "top",axis.text.x = element_text(angle = 45,hjust = 1)) + xlab("114 references Ordered by Length") + ylab("Length (AA)")
prop.table(table(referenceStat$Loop)) * 100

### SECTION 4: species distribution of references (Supplementary Figures 1b) ###
table((referenceStat %>% dplyr::filter(LastRank == 's'))$species) %>% as.data.frame() %>% ggplot(aes(reorder(Var1,-Freq),weight=Freq)) + geom_bar(fill="springgreen4") + xlab("") + theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + geom_text(aes(y=Freq,label=Freq)) + ylab('Number')

### SECTION 5: species distribution of 550 GUSs (Supplementary Figures 2) ###
table(GUSsStat$species) %>% as.data.frame %>% dplyr::arrange(-Freq) %>% dplyr::filter(Freq >= 5) %>% ggplot(aes(reorder(Var1,-Freq),weight=Freq)) + geom_bar(fill="springgreen4") + xlab("") + theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + geom_text(aes(y=Freq,label=Freq))  + ylab("Number of gmGUS")

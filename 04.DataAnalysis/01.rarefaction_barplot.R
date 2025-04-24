# ######################################################################
#
## Rarefaction curves
## Barplot for 114 references and 550 gmGUSs
## Author: JR C
## Date: 01.01.2025
#
# ######################################################################
library(vegan)
library(dplyr)
library(ggplot2)

##### Get Data #####
group <- read.csv("00.rawdata/group.csv",header = T,row.names = 1)
GUSabun_ABS <- read.csv("00.rawdata/GUSabun_ABS.csv",header = T,row.names = 1)
referenceStat <- read.csv("00.rawdata/referenceStat.csv",header = T,row.names = 1)
GUSsStat <- read.csv("00.rawdata/GUSsStat.csv",header = T,row.names = 1)

##### Rarefaction Curves #####
data.spec <- specaccum(t(GUSabun_ABS), method = "random")  
plot(data.spec, ci.type = "poly", col = "blue", lwd = 2, ci.lty = 0,
     ci.col = "lightblue", xlab = "Number of samples",
     ylab = "Number of gmGUS")
data.spec <- specaccum(t(GUSabun_ABS[,group %>% dplyr::filter(Stage == "Healthy") %>% rownames()]), method = "random")  
plot(data.spec, ci.type = "poly", col = '#cecccb', lwd = 2, ci.lty = 0,
     ci.col = "lightblue", xlab = "Number of samples",
     ylab = "Number of gmGUS",add=T)
data.spec <- specaccum(t(GUSabun_ABS[,group %>% dplyr::filter(Stage %in% c('MP','S0','SI_II','SIII_IV ')) %>% rownames()]), method = "random")  
plot(data.spec, ci.type = "poly", col = '#f3764a', lwd = 2, ci.lty = 0,
     ci.col = "lightblue", xlab = "Number of samples",
     ylab = "Number of gmGUS",add=T)

##### barPlot for 114 references ######
referenceStat %>% tibble::rownames_to_column(var="ID") %>% ggplot(aes(reorder(ID,Length),weight=Length,fill=Loop)) + geom_bar(width = 1,color="grey") + theme_classic() + theme(legend.position = "top",axis.text.x = element_text(angle = 45,hjust = 1)) + xlab("114 references Ordered by Length") + ylab("Length (AA)")
table((referenceStat %>% dplyr::filter(LastRank == 's'))$species) %>% as.data.frame() %>% ggplot(aes(reorder(Var1,-Freq),weight=Freq)) + geom_bar(fill="springgreen4") + xlab("") + theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + geom_text(aes(y=Freq,label=Freq)) + ylab('Number')
prop.table(table(referenceStat$Loop)) * 100 

##### barPlot for 550 gmGUss ######
table(GUSsStat$species) %>% as.data.frame %>% dplyr::arrange(-Freq) %>% dplyr::filter(Freq >= 5) %>% ggplot(aes(reorder(Var1,-Freq),weight=Freq)) + geom_bar(fill="springgreen4") + xlab("") + theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + geom_text(aes(y=Freq,label=Freq))  + ylab("Number of gmGUS")


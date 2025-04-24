# ######################################################################
#
## Analysis for Loop Category
## Author: JR C
## Date: 01.01.2025
#
# ######################################################################
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

##### Length Distribution #####
GUSsStat %>% tibble::rownames_to_column(var="ID") %>% ggplot(aes(reorder(ID,Length),weight=Length,fill=Loop)) + geom_bar(width = 1) + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "top") + xlab(paste0(nrow(data)," gmGUS Sequences Ordered by Length")) + ylab("Length (AA)")  + scale_fill_manual(values = loop.colors)
prop.table(table(GUSsStat$Loop)) * 100 

##### Pie Chart #####
library(viridis)
gg.result <- list()
i<- 1
for (mid.loop in unique(GUSsStat$Loop)) {
  mid.data <- GUSsStat %>% dplyr::filter(Loop == mid.loop) %>% dplyr::select(genus) %>% table %>% as.data.frame() %>% dplyr::arrange(-Freq)
  if(length(which(as.vector(mid.data$genus) == 'Unclassified')) > 0){
    mid.data <- rbind(
      mid.data[-which(as.vector(mid.data$genus) == 'Unclassified'),],
      mid.data %>% dplyr::filter(genus == 'Unclassified')
    )
  }
  mid.others <- sum(mid.data$Freq) - ((head(mid.data,n=5))$Freq  %>% sum())
  mid.data <- head(mid.data,n=5)
  if(mid.loop == "Mini-Loop"){
    #mid.data$genus <- factor(mid.data$genus %>% as.vector(),levels = mid.data$genus %>% as.vector() %>% rev())
    #mid.col <- colorRampPalette(c(loop.colors[mid.loop],"white"))(10)[c(1,5,9)] %>% rev
  }else{
    mid.data$genus <- factor(mid.data$genus %>% as.vector(),levels = c(mid.data$genus %>% as.vector(),"Others") %>% rev())
    mid.data[6,] <- c("Others",mid.others)
    #mid.col <- c(viridis(n=8,alpha = 1,begin=0,end=1,direction = 1,option = "A")[3:7],"grey62") %>% rev
    mid.col <- c(colorRampPalette(c(loop.colors[mid.loop],"white"))(18)[c(1,5,9,13,17)],"grey62") %>% rev
  }
  mid.data$Freq <-mid.data$Freq %>% as.numeric()
  mid.label <- paste(mid.data$genus %>% as.vector(), paste('(',sprintf("%.1f",(mid.data$Freq / sum(mid.data$Freq)) * 100), '%)', sep = ''), sep = ' ') %>% rev
  gg.result[[mid.loop]] <- ggplot(mid.data, mapping = aes(x = 'Content', y = Freq, fill = genus)) + geom_bar(stat = 'identity', position = 'stack', width = 1) + coord_polar(theta = 'y') + theme_void() + scale_fill_manual(labels = mid.label,values =mid.col) + labs(title = mid.loop)  #c( palette()[2:6],"grey62") %>% rev
  #ggsave(paste0(mid.loop,".pdf"),gg.result[[mid.loop]],width = 4,height = 3)
  i <- i+1
}
ggarrange(gg.result$`No Loop`,gg.result$`Mini-Loop 1`,gg.result$`Loop 1`,gg.result$`Loop 2`,gg.result$`Mini-Loop 2`,gg.result$`Mini-Loop 1,2`,nrow = 3,ncol = 2)

##### Difference of Abundance #####
loop.abun <- matrix(nrow=length(unique(GUSsStat$Loop)),ncol = ncol(GUSabun_TPM))
rownames(loop.abun) <- unique(GUSsStat$Loop)
colnames(loop.abun) <- colnames(GUSabun_TPM)
for (mid.loop in unique(GUSsStat$Loop)) {
  loop.abun[mid.loop,] <- GUSabun_TPM[GUSsStat %>% dplyr::filter(Loop == mid.loop) %>% rownames(),] %>% colSums() %>% as.vector()
}
loop.abun <- loop.abun %>% t() %>% as.data.frame()

#KW test
loop.abunKW <- loop.abun
loop.abunKW$Group <- group[rownames(loop.abun),'Stage']
loop.abunKW <- loop.abunKW %>% dplyr::filter(Group != 'NA')
kruskal.test(`No Loop` ~ Group, data = loop.abunKW)
kruskal.test(`Loop 1` ~ Group, data = loop.abunKW) 
kruskal.test(`Loop 2` ~ Group, data = loop.abunKW) 
kruskal.test(`Mini-Loop 1` ~ Group, data = loop.abunKW) 
kruskal.test(`Mini-Loop 2` ~ Group, data = loop.abunKW) 
kruskal.test(`Mini-Loop 1,2` ~ Group, data = loop.abunKW) 

mid.data <- loop.abun
mid.data$Group <- group[rownames(mid.data),'Stage']
mid.data <- mid.data %>% dplyr::filter(Group != "NA" & ! is.na(Group))
mid.data <- reshape2::melt(mid.data,by="Group")
mid.data$Group <- factor(mid.data$Group,levels = analysis.list.info$levels$Stage)
mid.data$variable <- factor(as.vector(mid.data$variable),levels = names(loop.colors))
mid.data %>% ggplot(aes(Group,value,fill=Group)) + geom_violin() + geom_boxplot(col='black',size=1,width =0.1,outlier.colour = NA) + scale_fill_manual(values = analysis.list.info$colors$Stage)+ theme_classic() + facet_wrap("variable",nrow = 1,ncol = 7,scales = "free") + geom_signif(comparisons = analysis.list.info$compaired$Stage, color="black",step_increase = 0.08, map_signif_level=function(p) if(p < 0.01){sprintf("%.3f",p)}else{if(p < 0.05){sprintf("%.3f",p)}else{sprintf("%.3f",p)}}, test = wilcox.test) + xlab("") + ylab("Abundance") + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")

##### Difference of GUS number #####
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
result %>% reshape2::melt("Group") %>% ggplot(aes(Group,value,fill=Group)) + geom_violin() + geom_boxplot(col='black',size=1,width =0.1,outlier.colour = NA) + scale_fill_manual(values = analysis.list.info$colors[['Stage']])+ theme_classic() + facet_wrap("variable",nrow = 1,ncol = 7,scales = "free") + geom_signif(comparisons = analysis.list.info$compaired[['Stage']], color="black",step_increase = 0.08, map_signif_level=function(p) if(p < 0.01){'+'}else{if(p < 0.05){'*'}else{sprintf("%.2f",p)}}, test = wilcox.test) + xlab("") + ylab("Number of gmGUS") + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")

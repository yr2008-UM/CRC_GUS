# ######################################################################
#
## enzyme assay + cell experiment + RNAseq (Figure 6)
## Author: JR C
## Date: 01.01.2025
#
# ######################################################################
library(ggplot2)
library(dplyr)
library(ggsci)
library(ggsignif)
library(ggpubr)
library(tibble)
library(readxl)
library(ggbreak)
library(Rmisc)
library(ggprism)

##### enzyme assay #####
data <- read_xlsx("./00.rawdata/in vitro/17GUSformat.xlsx",sheet = 1,col_names = T) %>% as.data.frame() 
data <- data %>% dplyr::mutate(Dot1=(((Dot1/100)*300)/CC)/(Time/60),Dot2=(((Dot2/100)*300)/CC)/(Time/60),Dot3=(((Dot3/100)*300)/CC)/(Time/60)) %>% dplyr::select(-Time,-CC) %>% reshape2::melt(by='GUS') %>% summarySE( measurevar="value", groupvars=c("GUS","Sub"))
orderGUS <- (data %>% dplyr::group_by(GUS) %>% dplyr::summarise(value=mean(value)) %>% dplyr::arrange(desc(value)))$GUS
data$GUS <- factor(as.vector(data$GUS),levels = orderGUS)
ggplot(data, aes(x = GUS, y = (value), fill = Sub)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin = (value), ymax = (value) + (se)),
                width=.2, position=position_dodge(.9)) +
  xlab("") + ylab("uM/h per uM GUS") +
  theme_pubclean() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggplot(data, aes(x = GUS, y = (value), fill = Sub)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin = (value), ymax = (value) + (se)),
                width=.2, position=position_dodge(.9)) +
  xlab("") + ylab("uM/h per uM GUS") +
  theme_pubclean() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0,1600)

data <- read_xlsx("./00.rawdata/in vitro/17GUSformat.xlsx",sheet = 2,col_names = T) %>% as.data.frame() %>% reshape2::melt(id=c('GUS','Sub','Time'))
data$Time <- as.numeric(data$Time)
data <- summarySE(data, measurevar="value", groupvars=c("GUS","Sub","Time"))
data$GUS <- factor(as.vector(data$GUS),levels = orderGUS)
ggplot(data,aes(Time,value,fill=Sub,shape=Sub,color=Sub)) + geom_point() + geom_line() + geom_errorbar(aes(ymin = value-se, ymax = value + se), width=.2, position=position_dodge(.9)) + xlab("") + ylab("uM/h per uM GUS") + facet_wrap('GUS',ncol = 5) + theme_bw(base_size = 12) + xlim(0,300) + scale_color_manual(values = pal_lancet()(2))

##### cell experiment #####
#HCT116
data <- read_xlsx("./00.rawdata/in vitro/BC-S3.cell experiment.xlsx",sheet = 1,col_names = T) %>% as.data.frame()
d1 <- data[,c('Time','Rep','Day','Control')] %>% dplyr::mutate(Mark = 'Control')
colnames(d1) <- c('Time','Rep','Day','Value','Mark')
d2 <- data[,c('Time','Rep','Day','BC.S3')] %>% dplyr::mutate(Mark = 'BC.G3')
colnames(d2) <- c('Time','Rep','Day','Value','Mark')
data <- rbind(d1,d2)
compare_means(Value~Mark, data=data, group.by = "Day",method = 'anova')
data$Mark <- factor(as.vector(data$Mark),levels = c("Control",'BC.G3'))
ggline(data, x = "Day", y = "Value", color = "Mark",shape='Mark', palette = 'lancet',add = "mean_sd",ylab = "OD450 value",ggtheme = theme_prism(base_size = 12)) +stat_compare_means(label = "p.signif",aes(group=Mark),method = "anova",  hide.ns = TRUE)
data %>% dplyr::group_by(Day,Mark) %>% dplyr::summarise(Value=mean(Value)) %>% as.data.frame()

#HT29
data <- read_xlsx("./00.rawdata/in vitro/BC-S3.cell experiment.xlsx",sheet = 2,col_names = T) %>% as.data.frame()
d1 <- data[,c('Time','Rep','Day','Control')] %>% dplyr::mutate(Mark = 'Control')
colnames(d1) <- c('Time','Rep','Day','Value','Mark')
d2 <- data[,c('Time','Rep','Day','BC.S3')] %>% dplyr::mutate(Mark = 'BC.G3')
colnames(d2) <- c('Time','Rep','Day','Value','Mark')
data <- rbind(d1,d2)
compare_means(Value~Mark, data=data, group.by = "Day",method = 'anova')
data$Mark <- factor(as.vector(data$Mark),levels = c("Control",'BC.G3'))
ggline(data, x = "Day", y = "Value", color = "Mark",shape='Mark', palette = 'lancet',add = "mean_sd",ylab = "OD450 value",ggtheme = theme_prism(base_size = 12)) +stat_compare_means(label = "p.signif",aes(group=Mark),method = "anova",  hide.ns = TRUE)
data %>% dplyr::group_by(Day,Mark) %>% dplyr::summarise(Value=mean(Value)) %>% as.data.frame()

##### volcano plots #####
library(EnhancedVolcano)
library(ggbreak)
df <- read.table("./00.rawdata/RNAseq/BC_G3vsControl_deg_all.xls",sep = "\t",header = T,row.names = 1,quote = "") 
keyvals <- ifelse(
  df$log2FoldChange < -log2(1.2) & df$padj < 0.05, 'royalblue3',
  ifelse(df$log2FoldChange > log2(1.2) & df$padj <0.05, 'red3','grey87'))
keyvals[is.na(keyvals)] <- 'grey87'
names(keyvals)[keyvals == 'red3'] <- 'Up'
names(keyvals)[keyvals == 'grey87'] <- 'nosig'
names(keyvals)[keyvals == 'royalblue3'] <- 'down'
EnhancedVolcano(df,
                lab = df$gene_name,
                legendPosition = 'none',
                x = 'log2FoldChange',
                y = 'padj',
                subtitle = NULL,
                title = "BC_G3 vs Control",
                # ylim = c(0, 50),
                xlim = c(-2.5, 2.5),
                axisLabSize = 12,
                FCcutoff = log2(1.2), 
                pCutoff = 0.05,
                pointSize = c(ifelse(abs(df$log2FoldChange) > log2(1.2) & df$padj < 0.05, 2, 1)),
                labSize = 2,
                caption = paste0("total = ", nrow(df), "; Up = ",length(names(keyvals)[names(keyvals) == 'Up']),"; Down = ",length(names(keyvals)[names(keyvals) == 'down'])),
                colCustom = keyvals,
                colAlpha = 0.8,
                gridlines.major = F, 
                gridlines.minor = F,
                border = 'full', 
                borderWidth = 0.5, 
                borderColour = 'black')+
  scale_y_break(c(13,38))

##### metascape #####
data <- read_xlsx("./00.rawdata/RNAseq/BC.G3.Up.metascape_result.xlsx",sheet = 2,col_names = T) %>% as.data.frame() %>% dplyr::filter(GroupID %in% paste0(c(1:20),"_Summary")) %>% dplyr::filter(`Log(q-value)` < log10(0.05))
colnames(data)[6] <- 'log10Q'
data$InTerm_InList <- as.numeric(gsub("\\/-","",data$InTerm_InList))
data$Type <- 'Up'
data2 <- read_xlsx("./00.rawdata/RNAseq/BC.G3.Down.metascape_result.xlsx",sheet = 2,col_names = T) %>% as.data.frame() %>% dplyr::filter(GroupID %in% paste0(c(1:20),"_Summary")) %>% dplyr::filter(`Log(q-value)` < log10(0.05))
colnames(data2)[6] <- 'log10Q'
data2$InTerm_InList <- as.numeric(gsub("\\/-","",data2$InTerm_InList))
data2$Type <- 'Down'
data <- rbind(data,data2)
data$Type <- factor(as.vector(data$Type),levels = c('Up','Down'))

ggplot(data,aes(reorder(Description,-log10Q), -log10Q,fill=log10Q))+
  geom_bar(stat = "identity")+
  facet_grid('Type',scales = "free",space = "fixed") + 
  geom_text(aes(label=InTerm_InList, y=-log10Q+1),size=3)+
  coord_flip()+
  labs(x='',y='-log10(q value)')+
  #scale_fill_manual(values = c('#852f88','#eb990c','#0f8096'))+
  scale_fill_gradient(low = "red3",high = "royalblue3")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        plot.margin=unit(x=c(top.mar=0.2,right.mar=0.2,bottom.mar=0.2,left.mar=0.2),units="inches"))

##### single cell #####
expressionMat <- read.table("./00.rawdata/RNAseq/SC_seq/heatmap.gct",sep = "\t",skip = 2,header = T,row.names = 1)
expressionMat <- t(expressionMat) %>% as.data.frame()
test <- expressionMat %>% dplyr::select(-SpecimenType)
Pvalue <- apply(test,2,function(x){ P <- wilcox.test(as.numeric(x[1:50000]),as.numeric(x[50001:100000]));return(P$p.value)}) %>% as.numeric()

meanExp <- read.table("./00.rawdata/RNAseq/SC_seq/expression.gct",sep = "\t",skip = 2,header = T,row.names = 1) %>% as.data.frame() %>% dplyr::filter(N > 0 | T > 0)
log2FC <- meanExp %>% dplyr::mutate(log2FC=ifelse(T > N,log2(T/N),-log2(N/T)))
meanExp <- as.data.frame(t(scale(t(meanExp))))%>% tibble::rownames_to_column('ID') %>% reshape2::melt(by='ID') %>% dplyr::rename(Expression=3)
perc <- read.table("./00.rawdata/RNAseq/SC_seq/percent.gct",sep = "\t",skip = 2,header = T,row.names = 1) %>% as.data.frame() %>% tibble::rownames_to_column('ID') %>% reshape2::melt(by='ID') %>% dplyr::rename(Perc=3)
selGenes <- rownames(log2FC)

inTable <- read.table("./00.rawdata/RNAseq/BC_G3vsControl_deg_all.xls",sep = "\t",header = T,quote = "") %>% dplyr::filter(gene_name %in% selGenes)  %>% dplyr::arrange(desc(log2FoldChange))
orderGenes <- as.vector(inTable$gene_name)
p1 <- inTable %>% ggplot(aes(gene_name,weight=log2FoldChange,fill=log2FoldChange)) + geom_bar(color='black') + theme_pubclean() + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + scale_fill_gradient2(low = "royalblue3",mid = 'white',high = "red3") + ylab("log2Fold Change") + xlim(orderGenes)
inTable <- (read.table("./00.rawdata/RNAseq/BC_G3vsControl_deg_all.xls",sep = "\t",header = T,quote = "") %>% dplyr::filter(gene_name %in% selGenes) %>% tibble::column_to_rownames('gene_name'))[,c("Control_1","Control_2","Control_3" ,"BC_G3_1","BC_G3_2","BC_G3_3")]
p2 <- as.data.frame(t(scale(t(inTable))))%>% tibble::rownames_to_column('ID') %>% reshape2::melt(by='ID') %>% ggplot(aes(ID,variable,fill=value)) + geom_tile() + theme_pubclean() + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + scale_fill_gradient2(low = "royalblue3",mid = 'white',high = "red3")+ xlim(orderGenes)
p3<- merge(meanExp,perc,by=c('ID','variable'),all=T) %>% dplyr::filter(Expression != 'NA') %>% ggplot(aes(ID,variable,color=Expression,size=Perc)) + geom_point() + theme_pubclean() + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + scale_color_gradient(low = "royalblue3",high = "red3") + xlim(orderGenes)
p4 <- log2FC %>% tibble::rownames_to_column('ID') %>% ggplot(aes(ID,weight=log2FC,fill=log2FC)) + geom_bar(color='black') + theme_pubclean() + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + scale_fill_gradient2(low = "royalblue3",mid = 'white',high = "red3") + ylab("log2Fold Change") + xlim(orderGenes)
cowplot::plot_grid(p1,p2,p3,p4,align = 'v',ncol = 1)

##### GSEA analysis #####
library(org.Hs.eg.db)
library(clusterProfiler)

inTable <- read.table("./00.rawdata/RNAseq/BC_G3vsControl_deg.xls",sep = "\t",header = T,row.names = 1,quote = "")
BCvsCTL <- inTable %>% dplyr::arrange(desc(log2FoldChange)) %>% dplyr::select(log2FoldChange)  %>% tibble::rownames_to_column('ENSEMBL')
geneFC <- as.vector(BCvsCTL$log2FoldChange)
names(geneFC) <- as.vector(BCvsCTL$ENSEMBL)
GOresult <- gseGO(geneFC,ont='ALL',OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL', pvalueCutoff = 1, pAdjustMethod = 'fdr')
inputTable <- GOresult %>% as.data.frame() %>% dplyr::select(-core_enrichment,-leading_edge) %>% dplyr::filter(p.adjust < 0.05 & abs(NES) > 1)  %>% dplyr::arrange(ONTOLOGY,desc(NES))
inputTable$Description <- factor(as.vector(inputTable$Description),levels = rev(as.vector(inputTable$Description)))
col <- rev(c(rep('#852f88',inputTable %>% dplyr::filter(ONTOLOGY == 'BP') %>% nrow()),rep('#eb990c',inputTable %>% dplyr::filter(ONTOLOGY == 'CC') %>% nrow()),rep('#0f8096',inputTable %>% dplyr::filter(ONTOLOGY == 'MF') %>% nrow())))
ggplot(inputTable,aes(Description, NES))+
  geom_bar(aes(fill=ONTOLOGY),stat = "identity")+
  #geom_text(aes(label=enrichmentScore, y=enrichmentScore+5),size=3)+
  coord_flip()+
  labs(x='',y='NES', title = 'BC_G3 vs Control')+
  scale_fill_manual(values = c('#852f88','#eb990c','#0f8096'))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.y = element_text(size=rel(0.85),colour =col),
        plot.margin=unit(x=c(top.mar=0.2,right.mar=0.2,bottom.mar=0.2,left.mar=0.2),units="inches"))


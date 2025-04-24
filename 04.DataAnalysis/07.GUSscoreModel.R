# ######################################################################
#
## GUSscore analysis
## Author: JR C
## Date: 01.01.2025
#
# ######################################################################
library(dplyr)
library(caret)
library(Boruta)
library(survival)
library(survminer)
library(forestploter)
library(glmnet)
library(timeROC)
library(readxl)
library(vegan)

##### Get Data #####
group <- read.csv("00.rawdata/group.csv",header = T,row.names = 1)
GUSabun_TPM <- read.csv("00.rawdata/GUSabun_TPM.csv",header = T,row.names = 1)
signGUSs <- readRDS("00.rawdata/signGUSs.RDS")
OutComeData <- read_excel("00.rawdata/clinical_supp.xlsx")

##### prepare data #####
id1 <- group %>% dplyr::filter( Stage %in% c('SI_II','SIII_IV') & Subject_ID %in% (OutComeData)$MetagenomeID) %>% rownames()
id2 <- (group %>% dplyr::filter(Stage %in% c('SI_II','SIII_IV') & Subject_ID %in% (OutComeData)$MetagenomeID))$Subject_ID  %>% as.character()
inGUSs <- unique(c(signGUSs$signSpecies$Up$S12,signGUSs$signSpecies$Up$S34,signGUSs$signSpecies$Down$S12,signGUSs$signSpecies$Down$S34))
dataBGUS <- GUSabun_TPM[inGUSs,id1] %>% t()
inputTable <- cbind(
  dataBGUS,
  (OutComeData %>% tibble::column_to_rownames('MetagenomeID'))[id2,] %>% dplyr::rename(Time=13,Event=14) %>% dplyr::select(Time,Event) %>% dplyr::mutate(Status=ifelse(Event == '-',1,0))
) 

##### coxph #####
beta_co <- c()
P <- c()
HR <- c()
HR_lower <- c()
HR_higher <- c()
z_p <- c()
Wald_p <- c()
Likelihood_p <- c()
for (midGene in inGUSs) {
  coxMid <- coxph(as.formula(paste0('Surv(Time,Status)~',midGene)),inputTable)
  surv_mid <- summary(coxMid)
  P <- append(P, cox.zph(coxMid)$table[1,3])
  beta_co <- append(beta_co,surv_mid$coefficients[,1])
  HR <- append(HR,exp(surv_mid$coefficients[,1]))
  HR_lower <- append(HR_lower,surv_mid$conf.int[,3])
  HR_higher <- append(HR_higher,surv_mid$conf.int[,4])
  z_p <- append(z_p,surv_mid$coefficients[,5])
  Wald_p <- append(Wald_p,as.numeric(surv_mid$waldtest[3]))
  Likelihood_p <- append(Likelihood_p,as.numeric(surv_mid$logtest[3]))
}
uni_cox_re <- data.frame(SYMBOL=inGUSs,P=P,beta=beta_co,'Hazard_Ratio'=HR,'HR_lower'=HR_lower,HR_higher=HR_higher,'z_pvalue'=z_p,'Wald_pvalue'=Wald_p,'Likelihood_pvalue'=Likelihood_p)
uni_cox_re %>% dplyr::filter(z_pvalue < 0.05 | Wald_pvalue  < 0.05 | Likelihood_pvalue < 0.05 )

##### forest plot #####
uni_cox_re_sig <- uni_cox_re %>% dplyr::filter(z_pvalue < 0.05 & Hazard_Ratio != 1 & P > 0.05) %>% dplyr::arrange(-Hazard_Ratio) %>% dplyr::mutate(HR=round(Hazard_Ratio,2),pValue=round(z_pvalue,3)) %>% dplyr::select(SYMBOL,HR,pValue,Hazard_Ratio,HR_lower,HR_higher) 
uni_cox_re_sig$` ` <- paste(rep(" ", nrow(uni_cox_re_sig)), collapse = " ")
tm <- forest_theme(core = list(bg_params=list(fill = c("white"))), 
                   base_size=5,
                   summary_col = "black", 
                   arrow_label_just = "end", 
                   ci_pch=19,
                   ci_lty=3,
                   arrow_type = "closed") 
p1 <- forest(uni_cox_re_sig[,c(1,2,3,7)], 
             est = uni_cox_re_sig$Hazard_Ratio, 
             lower = uni_cox_re_sig$HR_lower, 
             upper = uni_cox_re_sig$HR_higher, 
             ci_column = 4,ref_line = 1,xlab='Hazard Ratio',
             theme = tm) 
p1$heights <- rep(unit(3, "mm"), nrow(p1)) 

##### model construction #####
featureSelection <- function(inputTable,seeds){
  result <- NULL
  for(seedIN in seeds){
    set.seed(seedIN)
    print(seedIN)
    train_index <- createDataPartition(inputTable$Status, p = 0.8, list = FALSE)
    train_data <- inputTable[train_index, ]
    test_data <- inputTable[-train_index, ]
    
    boruta <- Boruta(
      x = as.matrix(train_data[,setdiff(colnames(train_data),c('Time','Status'))]),
      y = Surv(as.numeric(inputTable[rownames(train_data),'Time']),as.numeric(inputTable[rownames(train_data),'Status'])),
      pValue=0.05, mcAdj=T,maxRuns=500)
    
    feature <- boruta$finalDecision %>% as.data.frame() %>% dplyr::rename(Type = 1) %>% dplyr::filter(Type %in% c('Tentative','Confirmed')) %>% tibble::rownames_to_column('ID')
    if(nrow(feature) > 0){
      if(length(as.vector(feature$ID)) == 1){
        feature$Imp <- mean(boruta$ImpHistory[,as.vector(feature$ID)])
      }else{
        feature$Imp <- as.numeric(apply(boruta$ImpHistory[,as.vector(feature$ID)],2,mean))
      }
      feature$Seed <- seedIN
      if(is.null(result)){
        result <- feature
      }else{
        result <- rbind(result,feature)
      }
    }
  }
  return(result)
}
checkModel <- function(inputTable,seeds){
  resultMatrix <- matrix(nrow = length(seeds),ncol = 9)
  i <- 1
  for(seedIN in seeds){
    set.seed(seedIN)
    print(seedIN)
    
    train_index <- createDataPartition(inputTable$Status, p = 0.8, list = FALSE)
    train_data <- inputTable[train_index, ]
    test_data <- inputTable[-train_index, ]
    
    train_x <- as.matrix(train_data[,setdiff(colnames(train_data),c('Time','Status'))]) 
    train_y <- survival::Surv(as.numeric(inputTable[rownames(train_data),'Time']),as.numeric(inputTable[rownames(train_data),'Status']))
    
    test_x <- as.matrix(test_data[,setdiff(colnames(test_data),c('Time','Status'))]) 
    test_y <- survival::Surv(as.numeric(inputTable[rownames(test_data),'Time']),as.numeric(inputTable[rownames(test_data),'Status']))
    
    cv_fit = cv.glmnet(train_x,train_y,nfolds=10,family="cox")
    cf <- coef(cv_fit,s=cv_fit$lambda.min)
    cf <- as.matrix(cf) %>% as.data.frame() %>% dplyr::rename(coef=1) %>% dplyr::filter(coef!=0)
    selectGenes <- rownames(cf)
    
    if(length(selectGenes) > 1){
      riskScore <- cbind(predict(cv_fit,train_x,s=cv_fit$lambda.min,family='cox'),as.numeric(inputTable[rownames(train_x),'Status']),as.numeric(inputTable[rownames(train_x),'Time']))
      colnames(riskScore)=c('riskScore','Status','Time')
      riskScore <- riskScore %>% as.data.frame() %>% dplyr::mutate(Event=ifelse(Status == 0,'Alive','Dead'))
      
      res.cut <- surv_cutpoint(riskScore, time = "Time", event = "Status",variables = c("riskScore"))
      cutPoint <- res.cut$cutpoint$cutpoint
      riskScore <- riskScore %>% dplyr::mutate(score=riskScore,riskScore=ifelse(riskScore > cutPoint,'high','low'))
      riskScore$Time <-  as.numeric(riskScore$Time)
      
      fit <- survdiff(Surv(Time, Status) ~ riskScore, data = riskScore)
      riskScore$score <- as.numeric(riskScore$score)
      riskScore <- riskScore %>% dplyr::arrange(score)  %>% tibble::rownames_to_column("sampleID") %>% tibble::rowid_to_column('id')
      ROC <- timeROC(T=riskScore$Time, #生存时间
                     delta=riskScore$Status,   #生存状态
                     marker=riskScore$score, #计算timeROC的变量
                     cause=1,                #阳性结局指标数值(1表示死亡)
                     weighting="marginal",   #计算方法，默认为marginal
                     times=c(12*5,12*6),       #时间点，选取1年，3年和5年的生存率
                     iid=TRUE)
      resultTrain <- c(as.numeric(fit$pvalue),as.numeric(ROC$AUC))
      
      riskScore <- cbind(predict(cv_fit,test_x,s=cv_fit$lambda.min,family='cox'),as.numeric(inputTable[rownames(test_x),'Status']),as.numeric(inputTable[rownames(test_x),'Time']))
      colnames(riskScore)=c('riskScore','Status','Time')
      riskScore <- riskScore %>% as.data.frame %>% dplyr::mutate(score=riskScore,riskScore=ifelse(riskScore > cutPoint,'high','low'))
      
      if(length(unique(riskScore$riskScore)) ==2){
        
        riskScore$Time <-  as.numeric(riskScore$Time)
        fit <- survdiff(Surv(Time, Status) ~ riskScore, data = riskScore)
        riskScore$score <- as.numeric(riskScore$score)
        riskScore <- riskScore %>% dplyr::arrange(score)  %>% tibble::rownames_to_column("sampleID") %>% tibble::rowid_to_column('id')
        
        e1 = tryCatch({
          timeROC(T=riskScore$Time, #生存时间
                  delta=riskScore$Status,   #生存状态
                  marker=riskScore$score, #计算timeROC的变量
                  cause=1,                #阳性结局指标数值(1表示死亡)
                  weighting="marginal",   #计算方法，默认为marginal
                  times=c(12*5,12*6),       #时间点，选取1年，3年和5年的生存率
                  iid=TRUE)
        },error = function(e){
          3
        })
        if(is.list(e1)){
          ROC <- e1
          resultTest <- c(as.numeric(fit$pvalue),as.numeric(ROC$AUC))
          resultMatrix[i,] <- c(resultTrain,resultTest,seedIN,length(selectGenes),paste0(selectGenes,collapse = '; '))
          i <- i+1
        }
      }
    }
  }
  colnames(resultMatrix) <- c("p.train",'AUC.Y5','AUC.Y6','p.test','testAUC.Y5','testAUC.Y6','seed','Num','Feature')
  resultMatrix <- as.data.frame(resultMatrix)
  return(resultMatrix)
}
data <- t(GUSabun_TPM[,rownames(inputTable)])
test <- data
test[test > 0] <- 1
data <- data[,colSums(test) > nrow(inputTable) * 0.1]
inputSelection <- cbind(data,inputTable[,c('Time','Status')])
featureSelectionResult <- featureSelection(inputSelection,c(1:500))
selTax <- as.vector((featureSelectionResult$ID %>% table %>% as.data.frame() %>% dplyr::rename(ID=1) %>% dplyr::filter(Freq > 500*0.1))$ID)

inputAfterSelection <- inputSelection[,c(selTax,'Time','Status')]
ResultModel <- checkModel(inputAfterSelection,c(1:10000))
ResultModel %>%  dplyr::filter(p.train < 0.05 & p.test < 0.05 & AUC.Y6 > testAUC.Y6 & testAUC.Y6 > 0.8) %>% dplyr::select(-Feature) %>% dplyr::filter(Num >= 5)

##### visualization #####
lastModelVis <- function(inputTable,seedIN){
  set.seed(seedIN)
  
  train_index <- createDataPartition(inputTable$Status, p = 0.8, list = FALSE)
  train_data <- inputTable[train_index, ]
  test_data <- inputTable[-train_index, ]
  
  train_x <- as.matrix(train_data[,setdiff(colnames(train_data),c('Time','Status'))]) 
  train_y <- survival::Surv(as.numeric(inputTable[rownames(train_data),'Time']),as.numeric(inputTable[rownames(train_data),'Status']))
  test_x <- as.matrix(test_data[,setdiff(colnames(test_data),c('Time','Status'))]) 
  test_y <- survival::Surv(as.numeric(inputTable[rownames(test_data),'Time']),as.numeric(inputTable[rownames(test_data),'Status']))
  
  cv_fit = cv.glmnet(train_x,train_y,nfolds=10,family="cox")
  cf <- coef(cv_fit,s=cv_fit$lambda.min)
  cf <- as.matrix(cf) %>% as.data.frame() %>% dplyr::rename(coef=1) %>% dplyr::filter(coef!=0)
  selectGenes <- rownames(cf)
  
  #visualization for model
  fitmodel <- glmnet(train_x, train_y,family = "cox")
  
  #train set
  riskScore <- cbind(predict(cv_fit,train_x,s=cv_fit$lambda.min,family='cox'),as.numeric(inputTable[rownames(train_x),'Status']),as.numeric(inputTable[rownames(train_x),'Time']))
  colnames(riskScore)=c('riskScore','Status','Time')
  riskScore <- riskScore %>% as.data.frame() %>% dplyr::mutate(Event=ifelse(Status == 0,'Alive','Dead'))
  res.cut <- surv_cutpoint(riskScore, time = "Time", event = "Status",variables = c("riskScore"))
  cutPoint <- res.cut$cutpoint$cutpoint
  riskScore <- riskScore %>% dplyr::mutate(score=riskScore,riskScore=ifelse(riskScore > cutPoint,'high','low'))
  riskScore$Time <-  as.numeric(riskScore$Time)
  riskScore$score <- as.numeric(riskScore$score)
  riskScore <- riskScore %>% dplyr::arrange(score)  %>% tibble::rownames_to_column("sampleID") %>% tibble::rowid_to_column('id')
  riskScoreTrain <- riskScore
  
  #test set
  riskScore <- cbind(predict(cv_fit,test_x,s=cv_fit$lambda.min,family='cox'),as.numeric(inputTable[rownames(test_x),'Status']),as.numeric(inputTable[rownames(test_x),'Time']))
  colnames(riskScore)=c('riskScore','Status','Time')
  riskScore <- riskScore %>% as.data.frame() %>% dplyr::mutate(Event=ifelse(Status == 0,'Alive','Dead'))
  riskScore <- riskScore %>% as.data.frame %>% dplyr::mutate(score=riskScore,riskScore=ifelse(riskScore > cutPoint,'high','low'))
  riskScore$Time <-  as.numeric(riskScore$Time)
  riskScore$score <- as.numeric(riskScore$score)
  riskScore <- riskScore %>% dplyr::arrange(score)  %>% tibble::rownames_to_column("sampleID") %>% tibble::rowid_to_column('id')
  riskScoreTest <- riskScore
  
  #vis for cutpoint
  p1 <- plot(res.cut, "riskScore", palette = "npg",xlab="riskScore")
  ggsave("14.GUSscore/lasso.cutpoint.pdf",cowplot::plot_grid(p1$riskScore$distribution,p1$riskScore$maxstat, ncol = 1, align = "v",rel_heights=c(1,1)),width = 6,height = 5)
  
  resultF1 <- matrix(nrow = 2,ncol = 8)
  for (midType in c('train','test')) {
    if(midType == 'train'){
      riskScore = riskScoreTrain
    }else{
      riskScore = riskScoreTest
    }
    
    fit <- survdiff(Surv(Time, Status) ~ riskScore, data = riskScore)
    
    p1 <- ggsurvplot(survfit(Surv(Time, Status) ~ riskScore, data = riskScore),data=riskScore,
                     pval = T, conf.int = TRUE,
                     linetype = "strata", # Change line type by groups
                     risk.table = TRUE,
                     surv.median.line = "hv", # Specify median survival
                     ggtheme = theme_bw()# Change ggplot2 theme
    )
    ggsave(paste0("14.GUSscore/",midType,"survplot.pdf"),cowplot::plot_grid(p1$plot,p1$table, ncol = 1, align = "v",rel_heights=c(3,1)),width = 6,height = 6)
    
    p1 <- riskScore %>% ggplot(aes(x=id,y=score,color=riskScore)) + 
      geom_vline(xintercept = min(which(riskScore$riskScore == 'high')),lty=4,lwd=0.25) + 
      geom_point() + theme_bw() + 
      theme(axis.text.x =element_text(size=5), axis.text.y=element_text(size=5),axis.title.x =element_text(size=6), axis.title.y=element_text(size=6) ) +
      theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.ticks.length.x  = unit(0, "pt"),plot.margin = margin()) + labs(x="",y="riskScore")  + scale_color_manual(values = c('#ff4757','#546de5'))
    p2 <- riskScore %>% ggplot(aes(id,Time,color=Event)) + geom_point(alpha=0.8) + theme_bw() + 
      geom_vline(xintercept = min(which(riskScore$riskScore == 'high')),lty=4,lwd=0.25) + 
      theme(axis.text.x =element_text(size=5), axis.text.y=element_text(size=5),axis.title.x =element_text(size=6), axis.title.y=element_text(size=6) ) + scale_color_manual(values = c('#ff4757','#546de5')) +
      theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.ticks.length.x  = unit(0, "pt"),plot.margin = margin()) + labs(x="",y="Time") 
    inHeatmap <- t(inputTable)[selectGenes,as.vector(riskScore$sampleID)]
    inHeatmap <-t(scale(t(inHeatmap)))
    #midclust <- pheatmap::pheatmap(inHeatmap,scale='none',show_colnames = F,filename = NA,clustering_method='average')
    midclust <- hclust(vegdist(t(inputTable)[selectGenes,as.vector(riskScore$sampleID)],method = 'bray'))
    print(midclust)
    p3 <- inHeatmap %>% reshape2::melt(by='row.names') %>% dplyr::rename(SYMBOL=1,sampleID=2) %>% ggplot(aes(sampleID,SYMBOL,fill=value)) + geom_tile() + theme_bw() + theme(axis.text.x =element_text(size=5), axis.text.y=element_text(size=5),axis.title.x =element_text(size=6), axis.title.y=element_text(size=6) ) + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.ticks.length.x  = unit(0, "pt"),plot.margin = margin()) + xlim(as.vector(riskScore$sampleID)) + scale_fill_gradient2(low = "#546de5",mid = "white",high = "#ff4757") + labs(x='',y='',fill='Expression') + ylim(midclust$labels[midclust$order])  #ylim(midclust$tree_row$labels[midclust$tree_row$order]) 
    ggsave(paste0("14.GUSscore/",midType,"multiplot.pdf"),cowplot::plot_grid(p1,p2,p3, ncol = 1, align = "v",rel_heights = c(1,1,0.7)),width = 4,height = 3)
    
    ROC <- timeROC(T=riskScore$Time, #生存时间
                   delta=riskScore$Status,   #生存状态
                   marker=riskScore$score, #计算timeROC的变量
                   cause=1,                #阳性结局指标数值(1表示死亡)
                   weighting="marginal",   #计算方法，默认为marginal
                   times=c(12*5,12*6),       #时间点，选取1年，3年和5年的生存率
                   iid=TRUE)
    pdf(paste0("14.GUSscore/",midType,"ROC.pdf"),width = 5,height = 5)
    plot(ROC,time=12*5, col="red", lty=1,lwd=2, title = "") 
    plot(ROC,time=12*6, col="blue", add=TRUE, lty=1,lwd=2)
    legend("bottomright",#图例画在右下角
           c(paste0("5-Year (AUC = ",round(ROC[["AUC"]][1],3)," )"),
             paste0("6-Year (AUC = ",round(ROC[["AUC"]][2],3)," )")),
           col=c("red","blue"), 
           lty=1,  lwd=2,  bty = "n" #o表示用框框把图例部分框起来，为默认。n表示不画框框
    )
    dev.off()
    
    ### recall, F1, accuracy #####
    riskScore <-  riskScore %>% dplyr::mutate(riskScore=ifelse(riskScore == 'low',0,1))
    stat <- SeSpPPVNPV(res.cut$cutpoint$cutpoint,T=riskScore$Time,delta=riskScore$Status,marker=riskScore$score,cause=1,  weighting="marginal",times=12*6,iid=T)
    
    Se <- stat$TP[2] #sensitity
    Sp <- 1- stat$FP[2] #specificity
    Pos <- as.numeric(stat$Stats[,'Cases']['t=72'])     #阳性病例数
    TotalNum <- as.numeric(stat$Stats[,'Cases']['t=72']) + as.numeric(stat$Stats[,'survivor at t']['t=72']) #总病例数
    TP = as.numeric(sprintf("%.0f",Se*Pos))
    FN = Pos - TP
    TN = as.numeric( Sp * (TotalNum - Pos))
    FP = TotalNum - Pos - TN
    precision <- TP / (TP + FP)
    recall <- TP / (TP + FN)
    f1_score <- 2 * (precision * recall) / (precision + recall)
    specificity <- TN / (TN + FP)
    accuracy <- (TP + TN) / (TP + TN + FP + FN)
    
    if(midType == 'train'){
      resultF1[1,] <- c(precision,recall,f1_score,accuracy,TP,FN,TN,FP)
    }else{
      resultF1[2,] <- c(precision,recall,f1_score,accuracy,TP,FN,TN,FP)
    }
  }
  rownames(resultF1) <- c('Train','Test')
  colnames(resultF1) <- c('precision','recall','f1','accuracy','TP','FN','TN','FP')
  resultF1 <- as.data.frame(resultF1)
  
  return(list(cv_fit=cv_fit,fit=fitmodel,res.cut=res.cut,riskScoreTrain=riskScoreTrain,riskScoreTest=riskScoreTest,features=selectGenes,resultF1=resultF1))
}
lastModel <- lastModelVis(inputAfterSelection,1108)
pdf("14.GUSscore/lasso.partialLD.pdf",width = 6,height = 5);plot(lastModel$cv_fit);dev.off()
pdf("14.GUSscore/lasso.coef.pdf",width = 6,height = 5);plot(lastModel$fit,xvar = "lambda");dev.off()
lastModel$resultF1

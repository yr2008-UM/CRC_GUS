# ######################################################################
#
## GUSscore Model Construction
## Author: Junru Chen, junru.chen2019@hotmail.com
## Version 1.0
## Date: 09.09.2025
#
## Description:
## This script implements the GUSscore model for predicting CRC survival 
## outcomes. Key components include:
## 1. Cox proportional hazards regression
## 2. GUSscore model construction using LASSO
## 3. Survival analysis and riskScore stratification
## 4. Model visualization: survival curves, time-dependent ROC curves, etc.
#
## Input Files:
## 1. 00.rawdata/group.csv - Sample metadata
## 2. 00.rawdata/GUSabun_TPM.csv - gmGUS abundance matrix
## 3. 00.rawdata/signGUSs.RDS - Significant gmGUSs
## 4. 00.rawdata/clinical_supp.xlsx - Clinical outcome data
#
## Outputs:
## 1. GUSscore model
## 2. Survival curves and riskScore stratification plots
## 3. Time-dependent ROC curves
## 4. fig. 3d-g, fig. S9, S10
#
## Dependencies:
## dplyr, survival, survminer, glmnet, timeROC, readxl,
## vegan, caret, Boruta, forestploter
#
## Usage:
## 1. Set the working directory using setwd("the directory of scripts")
## 2. Ensure all input files are in the specified paths
## 3. Ensure all required packages were installed
## 4. Run script sections sequentially in RStudio or R command line
## 5. Output visualizations will be generated in the working directory
#
# Tested on: R 4.3.2 (macOS Ventura 13.5)
#
# ######################################################################

# Install required packages if missing
required_pkgs <- c("dplyr", "survival", "survminer", "glmnet", "timeROC","readxl","vegan","caret","Boruta","forestploter")
install_missing <- required_pkgs[!required_pkgs %in% installed.packages()]
if(length(install_missing)) install.packages(install_missing)
# if failed, try BiocManager::install('package')

# Load required libraries
library(dplyr)        # Data manipulation
library(survival)     # Survival analysis
library(survminer)    # Calculation of cutpoint
library(glmnet)       # LASSO regression
library(timeROC)      # Time-dependent ROC curves
library(readxl)       # Excel file reading
library(vegan)        # Ecological diversity analysis
library(caret)        # Machine learning
library(Boruta)       # Feature selection
library(forestploter) # Forest plot

### SECTION 1: definition of functions ###

## Function: featureSelection
## Purpose: Select features using Boruta
## Inputs:
##   - inputTable: Dataframe with features and survival data
##   - seeds: Vector of random seeds
## Output: Dataframe of selected features
featureSelection <- function(inputTable,seeds){
  result <- NULL
  for(seedIN in seeds){
    set.seed(seedIN)
    print(seedIN)
    
    # Split data into training and testing sets (80:20)
    train_index <- createDataPartition(inputTable$Status, p = 0.8, list = FALSE)
    train_data <- inputTable[train_index, ]
    test_data <- inputTable[-train_index, ]
    
    # Perform Boruta feature selection
    boruta <- Boruta(
      x = as.matrix(train_data[,setdiff(colnames(train_data),c('Time','Status'))]),
      y = Surv(as.numeric(inputTable[rownames(train_data),'Time']),as.numeric(inputTable[rownames(train_data),'Status'])),
      pValue=0.05, mcAdj=T,maxRuns=500)
    
    # Extract confirmed and tentative features
    feature <- boruta$finalDecision %>% as.data.frame() %>% dplyr::rename(Type = 1) %>% dplyr::filter(Type %in% c('Tentative','Confirmed')) %>% tibble::rownames_to_column('ID')
    
    if(nrow(feature) > 0){
      # Calculate mean importance
      if(length(as.vector(feature$ID)) == 1){
        feature$Imp <- mean(boruta$ImpHistory[,as.vector(feature$ID)])
      }else{
        feature$Imp <- as.numeric(apply(boruta$ImpHistory[,as.vector(feature$ID)],2,mean))
      }
      feature$Seed <- seedIN
      # Combine results
      if(is.null(result)){
        result <- feature
      }else{
        result <- rbind(result,feature)
      }
    }
  }
  return(result)
}

## Function: checkModel
## Purpose: Train and evaluate model
## Inputs:
##   - inputTable: Dataframe with features and survival data
##   - seeds: Vector of random seeds
## Output: Dataframe of evaluation results
checkModel <- function(inputTable,seeds){
  resultMatrix <- matrix(nrow = length(seeds),ncol = 9)
  i <- 1
  for(seedIN in seeds){
    set.seed(seedIN)
    print(seedIN)
    
    # Split data into training and testing sets (80:20)
    train_index <- createDataPartition(inputTable$Status, p = 0.8, list = FALSE)
    train_data <- inputTable[train_index, ]
    test_data <- inputTable[-train_index, ]
    
    # Prepare data matrices
    train_x <- as.matrix(train_data[,setdiff(colnames(train_data),c('Time','Status'))]) 
    train_y <- survival::Surv(as.numeric(inputTable[rownames(train_data),'Time']),as.numeric(inputTable[rownames(train_data),'Status']))
    
    test_x <- as.matrix(test_data[,setdiff(colnames(test_data),c('Time','Status'))]) 
    test_y <- survival::Surv(as.numeric(inputTable[rownames(test_data),'Time']),as.numeric(inputTable[rownames(test_data),'Status']))
    
    # Perform cross-validated LASSO regression
    cv_fit = cv.glmnet(train_x,train_y,nfolds=10,family="cox")
    
    # Obtain coefficents and features
    cf <- coef(cv_fit,s=cv_fit$lambda.min)
    cf <- as.matrix(cf) %>% as.data.frame() %>% dplyr::rename(coef=1) %>% dplyr::filter(coef!=0)
    selectGenes <- rownames(cf)
    
    if (length(selectGenes) > 1) {
      # Calculate risk scores
      riskScore <- cbind(predict(cv_fit,train_x,s=cv_fit$lambda.min,family='cox'),as.numeric(inputTable[rownames(train_x),'Status']),as.numeric(inputTable[rownames(train_x),'Time']))
      colnames(riskScore)=c('riskScore','Status','Time')
      riskScore <- riskScore %>% as.data.frame() %>% dplyr::mutate(Event=ifelse(Status == 0,'Alive','Dead'))
      
      # Determine optimal cutpoint
      res.cut <- surv_cutpoint(riskScore, time = "Time", event = "Status",variables = c("riskScore"))
      cutPoint <- res.cut$cutpoint$cutpoint
      
      # Cut the riskScore according to cutpoint
      riskScore <- riskScore %>% dplyr::mutate(score=riskScore,riskScore=ifelse(riskScore > cutPoint,'high','low'))
      riskScore$Time <-  as.numeric(riskScore$Time)
      riskScore$score <- as.numeric(riskScore$score)
      
      # Survival analysis
      fit <- survdiff(Surv(Time, Status) ~ riskScore, data = riskScore)
      
      # Time-dependent ROC
      riskScore <- riskScore %>% dplyr::arrange(score)  %>% tibble::rownames_to_column("sampleID") %>% tibble::rowid_to_column('id')
      ROC <- timeROC(T=riskScore$Time, 
                     delta=riskScore$Status,   
                     marker=riskScore$score, 
                     cause=1,                
                     weighting="marginal",  
                     times=c(12*5,12*6),      
                     iid=TRUE)
      
      resultTrain <- c(as.numeric(fit$pvalue),as.numeric(ROC$AUC))
      
      # Test set evaluation
      riskScore <- cbind(predict(cv_fit,test_x,s=cv_fit$lambda.min,family='cox'),as.numeric(inputTable[rownames(test_x),'Status']),as.numeric(inputTable[rownames(test_x),'Time']))
      colnames(riskScore)=c('riskScore','Status','Time')
      riskScore <- riskScore %>% as.data.frame %>% dplyr::mutate(score=riskScore,riskScore=ifelse(riskScore > cutPoint,'high','low'))
      
      if(length(unique(riskScore$riskScore)) ==2){
        riskScore$Time <-  as.numeric(riskScore$Time)
        
        # Survival analysis for test set
        fit <- survdiff(Surv(Time, Status) ~ riskScore, data = riskScore)
        
        # Time-dependent ROC for test set
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
        
        # Save result
        if(is.list(e1)){
          ROC <- e1
          resultTest <- c(as.numeric(fit$pvalue),as.numeric(ROC$AUC))
          resultMatrix[i,] <- c(resultTrain,resultTest,seedIN,length(selectGenes),paste0(selectGenes,collapse = '; '))
          i <- i+1
        }
      }
    }
  }
  
  # Format results
  colnames(resultMatrix) <- c("p.train",'AUC.Y5','AUC.Y6','p.test','testAUC.Y5','testAUC.Y6','seed','Num','Feature')
  resultMatrix <- as.data.frame(resultMatrix)
  return(resultMatrix)
}

## Function: lastModelVis
## Purpose: Visualize final GUSscore model
## Inputs:
##   - inputTable: Dataframe with features and survival data
##   - seedIN: Random seed for reproducibility
## Output: Model visualization plots and List of results
lastModelVis <- function(inputTable,seedIN){
  set.seed(seedIN)
  
  # Split data into training and testing sets (80:20)
  train_index <- createDataPartition(inputTable$Status, p = 0.8, list = FALSE)
  train_data <- inputTable[train_index, ]
  test_data <- inputTable[-train_index, ]
  
  # Prepare data matrices
  train_x <- as.matrix(train_data[,setdiff(colnames(train_data),c('Time','Status'))]) 
  train_y <- survival::Surv(as.numeric(inputTable[rownames(train_data),'Time']),as.numeric(inputTable[rownames(train_data),'Status']))
  test_x <- as.matrix(test_data[,setdiff(colnames(test_data),c('Time','Status'))]) 
  test_y <- survival::Surv(as.numeric(inputTable[rownames(test_data),'Time']),as.numeric(inputTable[rownames(test_data),'Status']))
  
  # Perform cross-validated LASSO regression & Obtain coefficients
  cv_fit = cv.glmnet(train_x,train_y,nfolds=10,family="cox")
  cf <- coef(cv_fit,s=cv_fit$lambda.min)
  cf <- as.matrix(cf) %>% as.data.frame() %>% dplyr::rename(coef=1) %>% dplyr::filter(coef!=0)
  selectGenes <- rownames(cf)
  
  # Visualize for model
  fitmodel <- glmnet(train_x, train_y,family = "cox")
  
  ## Training Set Visualization ##
  # Calculate risk scores
  riskScore <- cbind(predict(cv_fit,train_x,s=cv_fit$lambda.min,family='cox'),as.numeric(inputTable[rownames(train_x),'Status']),as.numeric(inputTable[rownames(train_x),'Time']))
  colnames(riskScore)=c('riskScore','Status','Time')
  riskScore <- riskScore %>% as.data.frame() %>% dplyr::mutate(Event=ifelse(Status == 0,'Alive','Dead'))
  
  # Determine optimal cutpoint
  res.cut <- surv_cutpoint(riskScore, time = "Time", event = "Status",variables = c("riskScore"))
  cutPoint <- res.cut$cutpoint$cutpoint
  riskScore <- riskScore %>% dplyr::mutate(score=riskScore,riskScore=ifelse(riskScore > cutPoint,'high','low'))
  riskScore$Time <-  as.numeric(riskScore$Time)
  riskScore$score <- as.numeric(riskScore$score)
  
  # Arrange by score
  riskScore <- riskScore %>% dplyr::arrange(score)  %>% tibble::rownames_to_column("sampleID") %>% tibble::rowid_to_column('id')
  riskScoreTrain <- riskScore
  
  ## Test Set Visualization ##
  # Calculate risk scores
  riskScore <- cbind(predict(cv_fit,test_x,s=cv_fit$lambda.min,family='cox'),as.numeric(inputTable[rownames(test_x),'Status']),as.numeric(inputTable[rownames(test_x),'Time']))
  colnames(riskScore)=c('riskScore','Status','Time')
  riskScore <- riskScore %>% as.data.frame() %>% dplyr::mutate(Event=ifelse(Status == 0,'Alive','Dead'))
  riskScore <- riskScore %>% as.data.frame %>% dplyr::mutate(score=riskScore,riskScore=ifelse(riskScore > cutPoint,'high','low'))
  riskScore$Time <-  as.numeric(riskScore$Time)
  riskScore$score <- as.numeric(riskScore$score)
  
  # Arrange by score
  riskScore <- riskScore %>% dplyr::arrange(score)  %>% tibble::rownames_to_column("sampleID") %>% tibble::rowid_to_column('id')
  riskScoreTest <- riskScore
  
  ## Visualization: Cutpoint, fig. S9c ##
  p1 <- plot(res.cut, "riskScore", palette = "npg",xlab="riskScore")
  ggsave("lasso.cutpoint.pdf",cowplot::plot_grid(p1$riskScore$distribution,p1$riskScore$maxstat, ncol = 1, align = "v",rel_heights=c(1,1)),width = 6,height = 5)
  
  ## Visualization: Survival Analysis ##
  resultF1 <- matrix(nrow = 2,ncol = 8)
  for (midType in c('train','test')) {
    if(midType == 'train'){
      riskScore = riskScoreTrain
    }else{
      riskScore = riskScoreTest
    }
    
    # Survival difference
    fit <- survdiff(Surv(Time, Status) ~ riskScore, data = riskScore)
    
    # Kaplan-Meier plot, fig. 3f, fig. S10b
    p1 <- ggsurvplot(survfit(Surv(Time, Status) ~ riskScore, data = riskScore),data=riskScore,
                     pval = T, conf.int = TRUE,
                     linetype = "strata", # Change line type by groups
                     risk.table = TRUE,
                     surv.median.line = "hv", # Specify median survival
                     ggtheme = theme_bw()# Change ggplot2 theme
    )
    ggsave(paste0(midType,"survplot.pdf"),cowplot::plot_grid(p1$plot,p1$table, ncol = 1, align = "v",rel_heights=c(3,1)),width = 6,height = 6)
    
    # Risk score distribution, fig. 3e, fig. S10a
    p1 <- riskScore %>% ggplot(aes(x=id,y=score,color=riskScore)) + 
      geom_vline(xintercept = min(which(riskScore$riskScore == 'high')),lty=4,lwd=0.25) + 
      geom_point() + theme_bw() + 
      theme(axis.text.x =element_text(size=5), axis.text.y=element_text(size=5),axis.title.x =element_text(size=6), axis.title.y=element_text(size=6) ) +
      theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.ticks.length.x  = unit(0, "pt"),plot.margin = margin()) + labs(x="",y="riskScore")  + scale_color_manual(values = c('#ff4757','#546de5'))
    
    # Survival time distribution, fig. 3e, fig. S10a
    p2 <- riskScore %>% ggplot(aes(id,Time,color=Event)) + geom_point(alpha=0.8) + theme_bw() + 
      geom_vline(xintercept = min(which(riskScore$riskScore == 'high')),lty=4,lwd=0.25) + 
      theme(axis.text.x =element_text(size=5), axis.text.y=element_text(size=5),axis.title.x =element_text(size=6), axis.title.y=element_text(size=6) ) + scale_color_manual(values = c('#ff4757','#546de5')) +
      theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.ticks.length.x  = unit(0, "pt"),plot.margin = margin()) + labs(x="",y="Time") 
    
    # Feature heatmap, fig. 3e, fig. S10a
    inHeatmap <- t(inputTable)[selectGenes,as.vector(riskScore$sampleID)]
    inHeatmap <-t(scale(t(inHeatmap)))
    midclust <- hclust(vegdist(t(inputTable)[selectGenes,as.vector(riskScore$sampleID)],method = 'bray'))
    print(midclust)
    p3 <- inHeatmap %>% reshape2::melt(by='row.names') %>% dplyr::rename(SYMBOL=1,sampleID=2) %>% ggplot(aes(sampleID,SYMBOL,fill=value)) + geom_tile() + theme_bw() + theme(axis.text.x =element_text(size=5), axis.text.y=element_text(size=5),axis.title.x =element_text(size=6), axis.title.y=element_text(size=6) ) + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.ticks.length.x  = unit(0, "pt"),plot.margin = margin()) + xlim(as.vector(riskScore$sampleID)) + scale_fill_gradient2(low = "#546de5",mid = "white",high = "#ff4757") + labs(x='',y='',fill='Expression') + ylim(midclust$labels[midclust$order])  #ylim(midclust$tree_row$labels[midclust$tree_row$order]) 
    
    # Save combined plot
    ggsave(paste0(midType,"multiplot.pdf"),cowplot::plot_grid(p1,p2,p3, ncol = 1, align = "v",rel_heights = c(1,1,0.7)),width = 4,height = 3)
    
    ## Time-dependent ROC, fig. 3g ##
    ROC <- timeROC(T=riskScore$Time,
                   delta=riskScore$Status,
                   marker=riskScore$score, 
                   cause=1,  
                   weighting="marginal",
                   times=c(12*5,12*6), 
                   iid=TRUE)
    pdf(paste0(midType,"ROC.pdf"),width = 5,height = 5)
    plot(ROC,time=12*5, col="red", lty=1,lwd=2, title = "") 
    plot(ROC,time=12*6, col="blue", add=TRUE, lty=1,lwd=2)
    legend("bottomright",
           c(paste0("5-Year (AUC = ",round(ROC[["AUC"]][1],3)," )"),
             paste0("6-Year (AUC = ",round(ROC[["AUC"]][2],3)," )")),
           col=c("red","blue"), 
           lty=1,  lwd=2,  bty = "n" 
    )
    dev.off()
    
    ## Performance Metrics, fig. S10c ##
    riskScore <-  riskScore %>% dplyr::mutate(riskScore=ifelse(riskScore == 'low',0,1))
    stat <- SeSpPPVNPV(res.cut$cutpoint$cutpoint,T=riskScore$Time,delta=riskScore$Status,marker=riskScore$score,cause=1,  weighting="marginal",times=12*6,iid=T)
    Se <- stat$TP[2] #sensitity
    Sp <- 1- stat$FP[2] #specificity
    Pos <- as.numeric(stat$Stats[,'Cases']['t=72']) 
    TotalNum <- as.numeric(stat$Stats[,'Cases']['t=72']) + as.numeric(stat$Stats[,'survivor at t']['t=72'])
    TP = as.numeric(sprintf("%.0f",Se*Pos))
    FN = Pos - TP
    TN = as.numeric( Sp * (TotalNum - Pos))
    FP = TotalNum - Pos - TN
    precision <- TP / (TP + FP)
    recall <- TP / (TP + FN)
    f1_score <- 2 * (precision * recall) / (precision + recall)
    specificity <- TN / (TN + FP)
    accuracy <- (TP + TN) / (TP + TN + FP + FN)
    
    # Store results
    if(midType == 'train'){
      resultF1[1,] <- c(precision,recall,f1_score,accuracy,TP,FN,TN,FP)
    }else{
      resultF1[2,] <- c(precision,recall,f1_score,accuracy,TP,FN,TN,FP)
    }
  }
  
  # Format performance metrics
  rownames(resultF1) <- c('Train','Test')
  colnames(resultF1) <- c('precision','recall','f1','accuracy','TP','FN','TN','FP')
  resultF1 <- as.data.frame(resultF1)
  
  return(list(cv_fit=cv_fit,fit=fitmodel,res.cut=res.cut,riskScoreTrain=riskScoreTrain,riskScoreTest=riskScoreTest,features=selectGenes,resultF1=resultF1))
}

### SECTION 2: data loading ###

# Load input data and metadata
group <- read.csv("00.rawdata/group.csv",header = T,row.names = 1)
GUSabun_TPM <- read.csv("00.rawdata/GUSabun_TPM.csv",header = T,row.names = 1)
signGUSs <- readRDS("00.rawdata/signGUSs.RDS")
OutComeData <- read_excel("00.rawdata/clinical_supp.xlsx")

# Prepare analysis dataset for cox analysis
id1 <- group %>% dplyr::filter( Stage %in% c('SI_II','SIII_IV') & Subject_ID %in% (OutComeData)$MetagenomeID) %>% rownames()
id2 <- (group %>% dplyr::filter(Stage %in% c('SI_II','SIII_IV') & Subject_ID %in% (OutComeData)$MetagenomeID))$Subject_ID  %>% as.character()

# Select significant gmGUSs for analysis
inGUSs <- unique(c(signGUSs$signSpecies$Up$S12,signGUSs$signSpecies$Up$S34,signGUSs$signSpecies$Down$S12,signGUSs$signSpecies$Down$S34))

# Prepare abundance matrix
dataBGUS <- GUSabun_TPM[inGUSs,id1] %>% t()

# Prepare outcome data
inputTable <- cbind(
  dataBGUS,
  (OutComeData %>% tibble::column_to_rownames('MetagenomeID'))[id2,] %>% dplyr::rename(Time=13,Event=14) %>% dplyr::select(Time,Event) %>% dplyr::mutate(Status=ifelse(Event == '-',1,0))
)

### SECTION 3: cox regression analysis ###

# Initialize result vectors
beta_co <- c()
P <- c()
HR <- c()
HR_lower <- c()
HR_higher <- c()
z_p <- c()
Wald_p <- c()
Likelihood_p <- c()

# Perform Cox regression for each gmGUS
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

# Compile results
uni_cox_re <- data.frame(SYMBOL=inGUSs,P=P,beta=beta_co,'Hazard_Ratio'=HR,'HR_lower'=HR_lower,HR_higher=HR_higher,'z_pvalue'=z_p,'Wald_pvalue'=Wald_p,'Likelihood_pvalue'=Likelihood_p)

# Filter significant results (p < 0.05)
uni_cox_re %>% dplyr::filter(z_pvalue < 0.05 | Wald_pvalue  < 0.05 | Likelihood_pvalue < 0.05 )

# Prepare data for forest plot
uni_cox_re_sig <- uni_cox_re %>% dplyr::filter(z_pvalue < 0.05 & Hazard_Ratio != 1 & P > 0.05) %>% dplyr::arrange(-Hazard_Ratio) %>% dplyr::mutate(HR=round(Hazard_Ratio,2),pValue=round(z_pvalue,3)) %>% dplyr::select(SYMBOL,HR,pValue,Hazard_Ratio,HR_lower,HR_higher) 
uni_cox_re_sig$` ` <- paste(rep(" ", nrow(uni_cox_re_sig)), collapse = " ")

# Define plot theme
tm <- forest_theme(core = list(bg_params=list(fill = c("white"))), 
                   base_size=5,
                   summary_col = "black", 
                   arrow_label_just = "end", 
                   ci_pch=19,
                   ci_lty=3,
                   arrow_type = "closed") 

# Generate forest plot, fig. 3d
p1 <- forest(uni_cox_re_sig[,c(1,2,3,7)], 
             est = uni_cox_re_sig$Hazard_Ratio, 
             lower = uni_cox_re_sig$HR_lower, 
             upper = uni_cox_re_sig$HR_higher, 
             ci_column = 4,ref_line = 1,xlab='Hazard Ratio',
             theme = tm) 
p1$heights <- rep(unit(3, "mm"), nrow(p1)) # Adjust plot heights
ggsave("forestPlot.pdf",p1,width=5,height=10)

### SECTION 4: construction of GUSscore model ###

## Prepare Data for Model ##
data <- t(GUSabun_TPM[,rownames(inputTable)])
test <- data
test[test > 0] <- 1
data <- data[,colSums(test) > nrow(inputTable) * 0.1] # Exsited > 10% of sampeles
inputSelection <- cbind(data,inputTable[,c('Time','Status')])

## Feature Selection ##
featureSelectionResult <- featureSelection(inputSelection,c(1:500))
selTax <- as.vector((featureSelectionResult$ID %>% table %>% as.data.frame() %>% dplyr::rename(ID=1) %>% dplyr::filter(Freq > 500*0.1))$ID)
inputAfterSelection <- inputSelection[,c(selTax,'Time','Status')]

## Model Training and Evaluation ##
ResultModel <- checkModel(inputAfterSelection,c(1:10000))

## Filter Successful Models ##
ResultModel %>%  dplyr::filter(p.train < 0.05 & p.test < 0.05 & AUC.Y6 > testAUC.Y6 & testAUC.Y6 > 0.8) %>% dplyr::select(-Feature) %>% dplyr::filter(Num >= 5)
best_seed <- 1108

### SECTION 5: Model visualization ###

## Generate Final Visualizations to current working directory ##
lastModel <- lastModelVis(inputAfterSelection,1108)

# Save plots for model, fig. S9a
pdf("lasso.partialLD.pdf",width = 6,height = 5);plot(lastModel$cv_fit);dev.off()
pdf("lasso.coef.pdf",width = 6,height = 5);plot(lastModel$fit,xvar = "lambda");dev.off()

# Display performance metrics
lastModel$resultF1

# ######################################################################
#
## Model construction for binary classifiers
## Author: JR C
## Date: 01.01.2025
#
# ######################################################################
library(dplyr)
library(caret)
library(pROC)
library(Boruta)
library(randomForest)

##### Get Data #####
group <- read.csv("00.rawdata/group.csv",header = T,row.names = 1)
GUSabun_TPM <- read.csv("00.rawdata/GUSabun_TPM.csv",header = T,row.names = 1)

##### Functions #####
onlyBoruta <- function(data,seeds){
  result <- NULL
  for (seedIN in seeds) {
    set.seed(seedIN)
    print(seedIN)
    # 分层划分训练集和测试集（8:2）
    train_index <- createDataPartition(data$Group, p = 0.8, list = FALSE)
    train_data <- data[train_index, ]
    test_data <- data[-train_index, ]
    
    # 转换数据为矩阵格式
    boruta <- Boruta(Group~.,data=train_data,pValue=0.05, mcAdj=T,maxRuns=500)
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
RFforBinary <- function(data,seeds,ntree=500){
  result <- matrix(nrow = length(seeds),ncol = 5)
  i <- 1
  for (seedIN in seeds) {
    set.seed(seedIN)
    # 分层划分训练集和测试集（8:2）
    train_index <- createDataPartition(data$Group, p = 0.8, list = FALSE)
    train_data <- data[train_index, ]
    test_data <- data[-train_index, ]
    
    selected_vars <- setdiff(colnames(train_data),c('Group'))
    train_selected <- train_data
    test_selected <- test_data
    
    # 定义参数网格（优化mtry）
    tune_grid <- expand.grid(.mtry = seq(1,length(selected_vars),by=1))  # 根据变量数量调整
    ctrl <- trainControl(method = "cv",number = 5,classProbs = TRUE,summaryFunction = twoClassSummary)
    rf_model <- train( Group ~ .,data = train_selected,  method = "rf", trControl = ctrl,tuneGrid = tune_grid,ntree = ntree,metric = "ROC")
    
    # 使用最优参数构建最终模型
    final_rf <- randomForest(Group ~ .,data = train_selected,ntree = ntree,mtry = rf_model$bestTune$mtry,importance = TRUE)
    test_pred <- predict(final_rf, newdata = test_selected)
    test_prob <- predict(final_rf, newdata = test_selected, type = "prob")
    # conf_matrix <- confusionMatrix(test_pred, test_data$Group)
    # ROC曲线和AUC
    roc_obj <- roc(test_data$Group, test_prob[, "CRC"])
    roc1 <- as.numeric(auc(roc_obj))
    roc_obj <- roc(train_selected$Group, predict(final_rf, newdata = train_selected, type = "prob")[, "CRC"])
    roc2 <- as.numeric(auc(roc_obj))
    result[i,] <- c(seedIN,roc1,roc2,length(selected_vars),paste0(selected_vars,collapse = "; "))
    i <- i+1
  }
  colnames(result) <- c('SeedIN','Test','Train','Num','Feature')
  result <- as.data.frame(result)
  return(result)
}
AfterRFforBinary <- function(data,seeds,GUSin,ret=F,ntree=500,cohort=T){
  result <- matrix(nrow = length(seeds),ncol = 8)
  i <- 1
  for (seedIN in seeds) {
    set.seed(seedIN)
    # 分层划分训练集和测试集（8:2）
    train_index <- createDataPartition(data$Group, p = 0.8, list = FALSE)
    train_data <- data[train_index, ]
    test_data <- data[-train_index, ]
    
    selected_vars <- setdiff(colnames(train_data),c('Group'))
    train_selected <- train_data
    test_selected <- test_data
    
    # 定义参数网格（优化mtry）
    tune_grid <- expand.grid(.mtry = seq(1,length(selected_vars),by=2))  # 根据变量数量调整
    ctrl <- trainControl(method = "cv",number = 5,classProbs = TRUE,summaryFunction = twoClassSummary)
    rf_model <- train( Group ~ .,data = train_selected,  method = "rf", trControl = ctrl,tuneGrid = tune_grid,ntree = ntree,metric = "ROC")
    
    # 使用最优参数构建最终模型
    final_rf <- randomForest(Group ~ .,data = train_selected,ntree = ntree,mtry = rf_model$bestTune$mtry,importance = TRUE)
    test_pred <- predict(final_rf, newdata = test_selected)
    test_prob <- predict(final_rf, newdata = test_selected, type = "prob")
    conf_matrixTest <- confusionMatrix(test_pred, test_data$Group,mode = "everything",positive = 'CRC')
    conf_matrixTrain <- confusionMatrix(predict(final_rf, newdata = train_selected), train_selected$Group,mode = "everything",positive = 'CRC')
    # ROC曲线和AUC
    roc_obj <- roc(test_data$Group, test_prob[, "CRC"])
    roc1 <- as.numeric(auc(roc_obj))
    roc_obj <- roc(train_selected$Group, predict(final_rf, newdata = train_selected, type = "prob")[, "CRC"])
    roc2 <- as.numeric(auc(roc_obj))
    
    rocValidation <- c()
    rocs <- list()
    
    if(cohort){
      for (midCountry in c("FRA","AUS","GER")) {
        groupMuti <- read.csv("00.rawdata/Cohorts/Cohorts.group.csv",check.names = F,row.names = 1) %>% dplyr::filter(Country == midCountry) %>% dplyr::filter(Group != 'NA') %>% dplyr::mutate(Group = ifelse(Group == 'CTR','Healthy','CRC'))
        groupMuti$Group <- factor(as.vector(groupMuti$Group),levels = c('Healthy','CRC'))
        speciesAbun <- read.table("00.rawdata/Cohorts/Cohorts.BGUS.abun.TPM.txt",sep = "\t",header = T,row.names = 1,quote = "",check.names = F)[GUSin,rownames(groupMuti)] %>% t() %>% as.data.frame()
        roc_obj <- roc(groupMuti$Group, predict(final_rf, newdata = speciesAbun, type = "prob")[, "CRC"])
        midroc <- as.numeric(auc(roc_obj))
        rocValidation <- c(rocValidation,midroc)
        rocs[[midCountry]] = list(
          Group = groupMuti$Group,
          data = speciesAbun
        )
      }
    }
    
    result[i,] <- c(seedIN,roc1,roc2,length(selected_vars),rocValidation,paste0(selected_vars,collapse = "; "))
    i <- i+1
  }
  colnames(result) <- c('SeedIN','Test','Train','Num',"FRA","AUS","GER",'Feature')
  result <- as.data.frame(result)
  
  if(ret){
    return(list(model=final_rf,conf_matrixTrain=conf_matrixTrain,conf_matrixTest=conf_matrixTest,plot=list(
      testGroup = test_data$Group, 
      testProb = test_prob[, "CRC"],
      trainGroup = train_selected$Group, 
      trainData = train_selected,
      validation=rocs
    )))
  }else{
    return(result)
  }
}
RFforBinaryAdenoma <- function(data,seeds,ntree=500,ret=F){
  result <- matrix(nrow = length(seeds),ncol = 5)
  i <- 1
  for (seedIN in seeds) {
    set.seed(seedIN)
    # 分层划分训练集和测试集（8:2）
    train_index <- createDataPartition(data$Group, p = 0.8, list = FALSE)
    train_data <- data[train_index, ]
    test_data <- data[-train_index, ]
    
    selected_vars <- setdiff(colnames(train_data),c('Group'))
    train_selected <- train_data
    test_selected <- test_data
    
    # 定义参数网格（优化mtry）
    tune_grid <- expand.grid(.mtry = seq(1,length(selected_vars),by=1))  # 根据变量数量调整
    ctrl <- trainControl(method = "cv",number = 5,classProbs = TRUE,summaryFunction = twoClassSummary)
    rf_model <- train( Group ~ .,data = train_selected,  method = "rf", trControl = ctrl,tuneGrid = tune_grid,ntree = ntree,metric = "ROC")
    
    # 使用最优参数构建最终模型
    final_rf <- randomForest(Group ~ .,data = train_selected,ntree = ntree,mtry = rf_model$bestTune$mtry,importance = TRUE)
    test_pred <- predict(final_rf, newdata = test_selected)
    test_prob <- predict(final_rf, newdata = test_selected, type = "prob")
    # conf_matrix <- confusionMatrix(test_pred, test_data$Group)
    # ROC曲线和AUC
    roc_obj <- roc(test_data$Group, test_prob[, "Adenoma"])
    roc1 <- as.numeric(auc(roc_obj))
    roc_obj <- roc(train_selected$Group, predict(final_rf, newdata = train_selected, type = "prob")[, "Adenoma"])
    roc2 <- as.numeric(auc(roc_obj))
    result[i,] <- c(seedIN,roc1,roc2,length(selected_vars),paste0(selected_vars,collapse = "; "))
    i <- i+1
  }
  colnames(result) <- c('SeedIN','Test','Train','Num','Feature')
  result <- as.data.frame(result)
  if(ret){
    return(list(model=final_rf,plot=list(
      testGroup = test_data$Group, 
      testProb = test_prob[, "Adenoma"],
      trainGroup = train_selected$Group, 
      trainData = train_selected
    )))
  }else{
    return(result)
  }
}
AfterRFforBinarySpecies <- function(data,seeds,GUSin,ret=F,ntree=500){
  result <- matrix(nrow = length(seeds),ncol = 5)
  i <- 1
  for (seedIN in seeds) {
    set.seed(seedIN)
    # 分层划分训练集和测试集（8:2）
    train_index <- createDataPartition(data$Group, p = 0.8, list = FALSE)
    train_data <- data[train_index, ]
    test_data <- data[-train_index, ]
    
    selected_vars <- setdiff(colnames(train_data),c('Group'))
    train_selected <- train_data
    test_selected <- test_data
    
    # 定义参数网格（优化mtry）
    tune_grid <- expand.grid(.mtry = seq(1,length(GUSin),by=1))  # 根据变量数量调整
    ctrl <- trainControl(method = "cv",number = 5,classProbs = TRUE,summaryFunction = twoClassSummary)
    rf_model <- train( Group ~ .,data = train_selected,  method = "rf", trControl = ctrl,tuneGrid = tune_grid,ntree = ntree,metric = "ROC")
    
    # 使用最优参数构建最终模型
    final_rf <- randomForest(Group ~ .,data = train_selected,ntree = ntree,mtry = rf_model$bestTune$mtry,importance = TRUE)
    test_pred <- predict(final_rf, newdata = test_selected)
    test_prob <- predict(final_rf, newdata = test_selected, type = "prob")
    conf_matrixTest <- confusionMatrix(test_pred, test_data$Group,mode = "everything",positive = 'CRC')
    conf_matrixTrain <- confusionMatrix(predict(final_rf, newdata = train_selected), train_selected$Group,mode = "everything",positive = 'CRC')
    # ROC曲线和AUC
    roc_obj <- roc(test_data$Group, test_prob[, "CRC"])
    roc1 <- as.numeric(auc(roc_obj))
    roc_obj <- roc(train_selected$Group, predict(final_rf, newdata = train_selected, type = "prob")[, "CRC"])
    roc2 <- as.numeric(auc(roc_obj))
    
    result[i,] <- c(seedIN,roc1,roc2,length(selected_vars),paste0(selected_vars,collapse = "; "))
    i <- i+1
  }
  colnames(result) <- c('SeedIN','Test','Train','Num','Feature')
  result <- as.data.frame(result)
  
  if(ret){
    return(list(model=final_rf,conf_matrixTrain=conf_matrixTrain,conf_matrixTest=conf_matrixTest,plot=list(
      testGroup = test_data$Group, 
      testProb = test_prob[, "CRC"],
      trainGroup = train_selected$Group, 
      trainData = train_selected
    )))
  }else{
    return(result)
  }
}

# ######################################################################
# Healthy (n = 247) vs. CRC (n = 258)
# ######################################################################

##### Prepare Data #####
data <- GUSabun_TPM[,group %>% dplyr::filter(Stage %in% c('Healthy','S0','SI_II','SIII_IV')) %>% rownames()] %>% t() %>% as.data.frame()
data$Group <- factor(as.vector((group %>% dplyr::mutate(Stage = ifelse(Stage %in% c('S0','SI_II',"SIII_IV"),'CRC','Healthy')))[rownames(data),'Stage']), levels = c("Healthy",'CRC'))

##### Feature Selection ######
resultForOnlyBoruta <- onlyBoruta(data,c(1:2000))
selTax <- as.vector(((resultForOnlyBoruta %>% dplyr::filter(Seed %in% c(1:500)))$ID %>% table %>% as.data.frame() %>% dplyr::rename(ID=1) %>% dplyr::filter(Freq > 500*0.1) %>% dplyr::arrange(Freq))$ID)
data <- data[,c(selTax,'Group')]
selTaxCRC <- selTax

##### model evaluation #####
BinaryResultselTax <- RFforBinary(data,c(1500:2000))

##### test for AUS,FRA,GER #####
AfterRFforBinaryResult <- AfterRFforBinary(data,as.numeric((BinaryResultselTax %>% dplyr::filter(Train > Test & Test > 0.8))$Seed),selTax)
# seed = 1783

##### ntree:500 ~ 1000#####
modelLast <- AfterRFforBinary(data,1783,selTax,ret = T,ntree = 500)
modelLast1 <- AfterRFforBinary(data,1783,selTax,ret = T,ntree = 600)
modelLast2 <- AfterRFforBinary(data,1783,selTax,ret = T,ntree = 700)
modelLast3 <- AfterRFforBinary(data,1783,selTax,ret = T,ntree = 800)
modelLast4 <- AfterRFforBinary(data,1783,selTax,ret = T,ntree = 900)
modelLast5 <- AfterRFforBinary(data,1783,selTax,ret = T,ntree = 1000) #choose ntree = 1000

##### ROC #####
roc1 <- roc(modelLast5$plot$testGroup, modelLast5$plot$testProb)
roc2 <- roc(modelLast5$plot$validation$AUS$Group, predict(modelLast5$model, newdata = modelLast5$plot$validation$AUS$data, type = "prob")[, "CRC"])
roc3 <- roc(modelLast5$plot$validation$GER$Group, predict(modelLast5$model, newdata = modelLast5$plot$validation$GER$data, type = "prob")[, "CRC"])
roc4 <- roc(modelLast5$plot$validation$FRA$Group, predict(modelLast5$model, newdata = modelLast5$plot$validation$FRA$data, type = "prob")[, "CRC"])
plot(roc1,col='red')
plot(roc2,col='orange',add=T)
plot(roc3,col='yellow',add=T)
plot(roc4,col='green',add=T)

##### var importance #####
varImpPlot(modelLast5$model,n=length(selTax))

# ######################################################################
# Healthy (n = 247) vs. Adenoma (n = 66)
# ######################################################################

##### Prepare Data #####
data <- GUSabun_TPM[,group %>% dplyr::filter(Stage %in% c('Healthy','MP')) %>% rownames()] %>% t() %>% as.data.frame()
data$Group <- factor(as.vector((group %>% dplyr::mutate(Stage = ifelse(Stage == 'MP','Adenoma','Healthy')))[rownames(data),'Stage']), levels = c("Healthy",'Adenoma'))

##### feature selection ######
AdenomaResultForOnlyBoruta <- onlyBoruta(data,c(1:500))
selTax <- as.vector(((AdenomaResultForOnlyBoruta )$ID %>% table %>% as.data.frame() %>% dplyr::rename(ID=1) %>% dplyr::filter(Freq > 500*0.1) %>% dplyr::arrange(Freq))$ID)
data <- data[,c(selTax,'Group')]

##### model evaluation #####
AdenomaRFforBinaryResult <- RFforBinaryAdenoma(data,c(1:500))
AdenomaRFforBinaryResult %>% dplyr::filter(Train > Test & Test > 0.8 ) %>% head(n=1) #seed = 7

##### 优化模型纳入的变量及ntree，ntree = 500 #####
modelLast <- RFforBinaryAdenoma(data,7,ret = T,ntree = 500) #last
modelLast1 <- RFforBinaryAdenoma(data,7,ret = T,ntree = 600)
modelLast2 <- RFforBinaryAdenoma(data,7,ret = T,ntree = 700)
modelLast3 <- RFforBinaryAdenoma(data,7,ret = T,ntree = 800)
modelLast4 <- RFforBinaryAdenoma(data,7,ret = T,ntree = 900)
modelLast5 <- RFforBinaryAdenoma(data,7,ret = T,ntree = 1000)

##### 绘制ROC曲线 #####
roc1 <- roc(modelLast$plot$testGroup, modelLast$plot$testProb)
plot(roc1,col='blue')

##### 绘制var importance #####
varImpPlot(modelLast$model,n=length(selTax))


# ######################################################################
# Healthy (n = 247) vs. CRC (n = 258)
# selTaxCRC + Species
# ######################################################################

##### Prepare Data #####
dbInfo <- read.table("00.rawdata/mOTUs4/mOTUsv4.0.gtdb.taxonomy.80mv.tsv",sep = "\t",header = T,row.names = 2)
data <- read.table("00.rawdata/mOTUs4/merged.mOTUs.table.txt",sep = "\t",header = T,row.names = 1)[,group %>% dplyr::filter(Stage %in% c('Healthy','S0','SI_II','SIII_IV')) %>% rownames()]
data$Species <- as.vector(dbInfo[rownames(data),'species'])
data <- data %>% as.data.frame() %>% reshape2::melt('Species') %>% dplyr::group_by(Species,variable) %>% dplyr::summarise(meanS=sum(value)) %>% as.data.frame() %>% reshape2::dcast(Species ~ variable)
data <- data %>% as.data.frame() %>% dplyr::filter(Species != '') %>% tibble::column_to_rownames('Species')

test <- data
test[test > 0] <- 1
data <- data[rowSums(test) >= ncol(test) * 0.1,]

dataGUS <- GUSabun_TPM[selTaxCRC,group %>% dplyr::filter(Stage %in% c('Healthy','S0','SI_II','SIII_IV')) %>% rownames()] 
data <- rbind(data,dataGUS[,colnames(data)]) %>% t() %>% as.data.frame()
data$Group <- factor(as.vector((group %>% dplyr::mutate(Stage = ifelse(Stage %in% c('S0','SI_II',"SIII_IV"),'CRC','Healthy')))[rownames(data),'Stage']), levels = c("Healthy",'CRC'))
colnames(data) <- gsub(" ","_",colnames(data))
colnames(data) <- gsub("-","_",colnames(data))
colnames(data) <- gsub("\\(","_",colnames(data))
colnames(data) <- gsub("\\)","_",colnames(data))

##### feature selection #####
resultForOnlyBorutaSpecies <- onlyBoruta(data,c(1:500))
selTax <- as.vector(((resultForOnlyBorutaSpecies)$ID %>% table %>% as.data.frame() %>% dplyr::rename(ID=1) %>% dplyr::filter(Freq > 500*0.1) %>% dplyr::arrange(Freq))$ID)
data <- data[,c(selTax,'Group')]

##### model evaluation #####
SpeciesRFforBinaryResult <- RFforBinary(data,c(1:500))
SpeciesRFforBinaryResult %>% dplyr::filter(Train > Test & Test > 0.8 )  #seed = 155
modelLast <- AfterRFforBinarySpecies(data,155,selTax,ret = T,ntree = 500)

##### ROC #####
roc1 <- roc(modelLast$plot$testGroup, modelLast$plot$testProb)
plot(roc1,col='red')

##### var importance #####
varImpPlot(modelLast$model,n=length(selTax))

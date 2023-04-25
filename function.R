draw_confusion_matrix <- function(cm,title) {
  
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title(title, cex.main=2)
  
  # create the matrix 
  classes = colnames(cm$table)
  rect(150, 430, 240, 370, col='#3F97D0')
  text(195, 435, classes[1], cex=1.2)
  rect(250, 430, 340, 370, col='#F7AD50')
  text(295, 435, classes[2], cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col='#F7AD50')
  rect(250, 305, 340, 365, col='#3F97D0')
  text(140, 400, classes[1], cex=1.2, srt=90)
  text(140, 335, classes[2], cex=1.2, srt=90)
  
  # add in the cm results 
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}  




model_compare = function(train,test, seed){
  
  ## obtain train and test y 
  
  train.y = as.factor(grepl("B",rownames(train)))
  
  test.y = as.factor(grepl("B",rownames(test)))
  
  train = train[,which(colMeans(train == 0) < 1)]
  
  ## using LRT1D from ZINB Package
  
  ## get p values and q values after fdr correction
  
  ## fitting ZINB model
  # y: Species count
  # X: Status of ASD
  
  ## LRT for multiple responses to find important predictors

  
  q.value = LRTnD(as.numeric(as.logical(train.y)),train)$q.value
  
  ## get filter train data 
  
  train.x = train[,which(q.value < 0.05)]
  
  ## fit cross validated lasso model and predict with test data
  
  set.seed(seed)
  cv.lasso.fit <- cv.glmnet(x = train.x, y = train.y, alpha = 1, nfolds = 10,family = "binomial",type.measure = "class")
  fit.lasso = glmnet(x = train.x, y = train.y, alpha = 1,family = "binomial",lambda = cv.lasso.fit$lambda.min)
  pred.lasso = predict(fit.lasso,newx = test[,colnames(train.x)],type = "response") 
  
  ## fit cross validated random forest model and predict with test data
  # wrap in 10 fold CV
  # random forest with entire dataset using "caret" package
  # grid search for tuning parameters: # trees, # predictors considered at each split
  set.seed(seed)
  trControl = trainControl(method = "cv", number = 10, search ="grid")
  
  fit.rf = train(x = train, y = train.y,
                 method = "rf",
                 metric = "Accuracy",
                 trControl = trControl)
  pred.rf = predict(fit.rf,data.frame(test),type = "prob") 
  
  ## summarize the prediction result from lasso and rf model
  
  result = data.frame(test.y,pred.lasso,pred.rf[,2])
  colnames(result) = c("status","pred.lasso","pred.rf")
  
  return(result)
  
}
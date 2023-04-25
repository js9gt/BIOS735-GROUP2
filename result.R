## Loading Package
library(Rcpp)
library(dplyr)
library(ZINB)
library(caret)
library(glmnet)
library(pROC)

## load function
source("function.R")

## read in data
microbiome = t(read.csv("Data/microbiome.csv",row.names = 1)[,-1])

# split data into train and test
set.seed(1)
train.num = sort(sample(1:nrow(microbiome),0.8*nrow(microbiome)))
train = microbiome[train.num,]
test = microbiome[-train.num,]

## analyze the data
result = model_compare(train,test, seed = 1)


## generate ROC curve and AUC score

tiff("Figures/ROC.tiff", units="in", width=5, height=5, res=200)
roc(result$status,as.numeric(result$pred.lasso),plot=T,percent = T,legacy.axes = T,print.auc = T,col = "#377eb8",
    ylab = "True positive (%)",xlab = "False Positive %",main = "ROC")
plot.roc(result$status,as.numeric(result$pred.rf),percent = T,add=T,print.auc = T,print.auc.y = 40,col="#4daf4a")
legend("bottomright",legend = c("Lasso","Random forest"),col=c("#377eb8","#4daf4a"),lwd=4)
dev.off()

## generate information of confusion matrix
## calculate accuracy, specificity, sensitivity metrics for both model

result$status = relevel(result$status,"TRUE")
result$pred.lasso = as.factor(result$pred.lasso > 0.5)
result$pred.lasso  = relevel(result$pred.lasso,"TRUE")
result$pred.rf = as.factor(result$pred.rf > 0.5)
result$pred.rf  = relevel(result$pred.rf,"TRUE")

cm.lasso <- confusionMatrix(data = result$pred.lasso, reference = result$status)
cm.rf <- confusionMatrix(data = result$pred.rf, reference = result$status)

tiff("Figures/CM_Lasso.tiff", units="in", width=5, height=5, res=200)
draw_confusion_matrix(cm.lasso,title = "Lasso")
dev.off()

tiff("Figures/CM_Rf.tiff", units="in", width=5, height=5, res=200)
draw_confusion_matrix(cm.rf,title = "Random Forest")
dev.off()

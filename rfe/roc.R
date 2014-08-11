library("ROCR")
library("caret")
library("doMC")

loocv <- function(x, y, ...){
  require("randomForest")
  require("plyr")
  
  loocv.iter <- function(x, y, x.test, ...){
    rf = randomForest(x, y, ...)
    p = predict(rf, newdata=x.test, type="prob")
    return(p)
  }
  
  p = laply(1:nrow(x), function(i){
    x.train = x[-i, ]
    x.test = x[i, , drop=F]
    y.train = y[-i]
    
    loocv.iter(x.train, y.train, x.test)
  }, .parallel=T)
}

roc.plot <- function(p1, p2, y, ...){
    require("ROCR")
    require("RColorBrewer")
    
    roc1 = p1[,2] %>% prediction(y) %>% performance("tpr", "fpr")
    roc2 = p2[,2] %>% prediction(y) %>% performance("tpr", "fpr")
    
    auc1 = y %>% roc(p1[,2]) %>% auc
    auc2 = y %>% roc(p2[,2]) %>% auc
    
    col = brewer.pal(5, "Set2")[4:5]
    
    plot(roc1, 
         col=col[1], 
         lwd=2, ...)
    
    plot(roc2, 
         col=col[2],
         lwd=2,
         add=T)
    
    legend.text = paste(c("All Probe Model", "10 Probe Model"),
                         " (",
                         round(c(auc1, auc2), 2),
                         ")", sep="")
    
    legend(0.5, 
           0.5, 
          legend.text, 
          lwd=2, 
          col=col)
}

registerDoMC(10)
setwd("/gpfs/home/ekramer/Projects/CEC/rfe")
load("../data/rfe_data.Rdata")
load("../data/cec.Rdata")

ten.gene.training = genes[genes$Cohort != "VALIDATION", ]
ten.gene.validation = genes[genes$Cohort == "VALIDATION" & genes$SN > 100, ]

## ROC on training set


# 10 gene set training
y = training$Status
x = training[gene.rfe$denovo$optVariables]
p.10.genes = loocv(x, y, ntree=1000)

# full gene set training
y = training$Status

x = select(training, -Array, -Cohort, -Status, -SN)
p.all.genes = loocv(x, y, ntree=1000)

# plotting

roc.plot(p.all.genes,
         p.10.genes,
         y,
         main="Training Set")

## ROCR on validation set

rf.10.gene = gene.rfe$denovo
rf.all.genes = randomForest(x, y, ntree=1000)

p.10.genes = predict(rf.10.gene, newdata=validation)
p.all.genes = predict(rf.all.genes, newdata=validation, type="prob")

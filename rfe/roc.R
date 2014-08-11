library("ROCR")
library("caret")
library("doMC")

loocv <- function(x, y, ...){
  
  ## Performs leave-one-out cross-validation
  ## ... are parameters passed to randomForest
  
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
  
    ## Plots ROC curves for two sets of predictions
    ## p1 -- "All probes"
    ## p2 -- "10 probe model"
    ## y -- common set of labels for predictions
  
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
                         ")", 
                        sep="")
    
    legend(0.5, 
           0.5, 
          legend.text, 
          lwd=2, 
          col=col)
}

get.predictions <- function(x, y.train, ...){
  require("randomForest")
  
  x.train = x[x$Cohort != "VALIDATION", ]
  x.validation = x[x$Cohort == "VALIDATION", ]
  
  p.train = loocv(x.train, y.train, ...)  
  
  rf = randomForest(x.train, y.train, ...)
  p.validation = predict(rf, newdata=x.validation, type="prob")
  
  return(list(p.train=p.train, p.validation=p.validation))
}



registerDoMC(10)
setwd("/gpfs/home/ekramer/Projects/CEC/rfe")
load("../data/rfe_data.Rdata")
load("../data/cec.Rdata")

y.train = genes$Status

ten.gene.x = genes[gene.rfe$denovo$optVariables]
all.probes.x = select(probes, -Status, -Cohort, -Array, -SN)

  

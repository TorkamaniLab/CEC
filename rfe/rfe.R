library("caret")
library("randomForest")
library("doMC")
library("plyr")
library("dplyr")

run.rfe <- function(d){
  
  ### CUSTOM RFE FUNCTIONS
  bigSummary = function(...) c(twoClassSummary(...), defaultSummary(...))
  select10 = function(...) 10
  
  ### CUSTROM LIST FOR CARET'S RFE FUNCTION
  denovoFuncs = rfFuncs
  denovoFuncs$summary = bigSummary
  denovoFuncs$selectSize = select10
  
  ## SUBSET INTO TRAINING AND VALIDATION
  training = d[d$Cohort %in% c("TRAINING", "TESTING"), ]
  validation = d[d$Cohort == "VALIDATION", ]
  
  ## PREPARE TRAINING DATA 
  x = training %>% select(-Array, -Cohort, -Status, -SN)
  y = training$Status
  
  ind = sapply(x, function(z) median(z[y == "AMI"]) > median(z[y != "AMI"]))
  x = x[ , ind]
  
  ### RFE CONTROL OBJECTS
  index = createMultiFolds(y, times=5)
  
  denovo.ctrl = rfeControl(method="repeatedcv",
                           repeats=5,
                           verbose=T,
                           functions=denovoFuncs,
                           returnResamp="all",
                           saveDetails=F,
                           index=index)
  
  ### RUN RFE
  denovo = rfe(x=x,
               y=y,
               sizes=1:20,
               metric="ROC",
               rfeControl=denovo.ctrl,
               ntree=1000)
  
  p = predict(denovo, newdata=validation)
  auc.all = auc(roc(validation$Status, p[,2]))
  
  ind = validation$SN > 100
  auc.highSN = auc(roc(validation$Status[ind], p[ind,2]))
  
  return(list(denovo=denovo, auc.all=auc.all, auc.highSN))
}


### SCRIPT 
load("~/Projects/CEC/data/cec.Rdata")
registerDoMC(8)

probe.rfe = run.rfe(probes)
print(probe.rfe$denovo)
cat(paste("AUC on entire validation set:", probe.rfe$auc.all))
cat(paste("AUC on high SN validation set:", probe.rfe$auc.highSN))

gene.rfe = run.rfe(genes)
print(gene.rfe$denovo)
cat(paste("AUC on entire validation set:", gene.rfe$auc.all))
cat(paste("AUC on high SN validation set:", probe.rfe$auc.highSN))


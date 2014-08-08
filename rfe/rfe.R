library("caret")
library("randomForest")
library("doMC")
library("plyr")
library("dplyr")

### CUSTOM RFE FUNCTIONS

bigSummary = function(...) c(twoClassSummary(...), defaultSummary(...))

upregulatedFit = function(x, y, first, last, ...){
  require("data.table")
  require("reshape2")
  
  # x holds features
  # y holds labels
  if(class(x)[1] != "data.table") x = as.data.table(x)
  
  ## using data.table to find 
  ## down-regulated genes
  tmp = copy(x)
  tmp[, Status:=y]
  
  # melts data.table by Status
  tmp2 = melt(tmp, id.vars="Status", variable.name="Gene") 
  
  # finds mean expression for 
  # each gene by status
  tmp3 = tmp2[ ,
              list(mean.exp=mean(value)), 
              by=c("Status", "Gene")]
  tmp3[, max.exp:=max(mean.exp), by="Gene"]
  
  # finds down-regulated genes
  down.genes = tmp3[Status == "Ctrl", ][mean.exp == max.exp , Gene]
  down.genes = as.character(down.genes)
  
  # set down-regulated genes to zeros
  if(length(down.genes) > 0) x[ , (down.genes) := 0]
  
  # pass modified data to normal fit function
  rfFuncs$fit(x, y, first, last, ...)
}

select10 = function(x, metric, maximize) 10


### CUSTROM LIST FOR CARET'S RFE FUNCTION

denovoFuncs = rfFuncs
denovoFuncs$summary = bigSummary
∂denovoFuncs$fit = upregulatedFit
denovoFuncs$selectSize = select10

### RUN RFE

# process datas
load("~/Projects/CEC/data/exprs.Rdata")

d = merge(probe.exprs, pheno[c("PTID", "status", "cohort")])

training = d[d$cohort %in% c("TRAINING", "TESTING"), ]
validation = d[d$cohort == "VALIDATION", ]

x = training %>% select(-PTID, -cohort, -status)
y = factor(training$status)

index = createMultiFolds(y, times=5)

denovo.ctrl = rfeControl(method="repeatedcv",
                         repeats=5,
                         verbose=T,
                         functions=denovoFuncs,
                         returnResamp="all",
                         saveDetails=F,
                         index=index)
registerDoMC(15)
denovo = rfe(x=x,
             y=y,
             sizes=1:20,
             metric="ROC",
             rfeControl=denovo.ctrl,
             ntree=1000)

d = merge(probe.exprs, pheno[c("PTID", "cohort", "status", "sn")])
validation = d[d$cohort == "VALIDATION", ]
validation = validation[validation$sn > 100, ]

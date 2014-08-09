library("caret")
library("randomForest")
library("doMC")
library("plyr")
library("dplyr")

### CUSTOM RFE FUNCTIONS

bigSummary = function(...) c(twoClassSummary(...), defaultSummary(...))

select10 = function(x, metric, maximize) 10


### CUSTROM LIST FOR CARET'S RFE FUNCTION

denovoFuncs = rfFuncs
denovoFuncs$summary = bigSummary
denovoFuncs$selectSize = select10

### RUN RFE

# process datas
load("~/Projects/CEC/data/exprs.Rdata")

d = merge(probe.exprs, pheno[c("PTID", "status", "cohort")])

training = d[d$cohort %in% c("TRAINING", "TESTING"), ]
validation = d[d$cohort == "VALIDATION", ]

x = training %>% select(-PTID, -cohort, -status)
y = factor(training$status)

ind = sapply(x, function(z) median(z[y == "AMI"]) > median(z[y != "AMI"]))
x = x[ , ind]

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
validation2 = validation[validation$sn > 100, ]

p = predict(denovo, newdata=validation)
p2 = predict(denovo, newdata=validation2)

auc(roc(validation$status, p[,2]))
auc(roc(validation2$status, p2[,2]))


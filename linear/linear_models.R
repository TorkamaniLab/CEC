library("kernlab")
library("glmnet")
library("plyr")
library("dplyr")
library("pROC")

setwd("/gpfs/home/ekramer/Projects/CEC/linear/")

load("~/Projects/CEC/data/cec_filtered2.Rdata")

set.seed(1)

select = dplyr::select

## delete anti-sense RNA genes
genes = genes[ , !grepl("-AS[0-9]+$", colnames(genes))]

## get appropriate subsets
## x.t - training features
## x.v - validation features
## y.t - training labels
## y.v - validation labels

x.t = genes[genes$Cohort != "VALIDATION", ] %>%
  select(-Cohort, -SN, -Status, -Array)
y.t = genes[genes$Cohort != "VALIDATION", ]$Status

x.v = genes[genes$Cohort == "VALIDATION" & genes$SN > 100, ] %>%
  select(-Cohort, -SN, -Status, -Array)
y.v = genes[genes$Cohort == "VALIDATION" & genes$SN > 100, ]$Status

## select probes up-regulated in AMI
ind = sapply(x.t, function(z) median(z[y.t == "AMI"]) > 1 + median(z[y.t !="AMI"]))
x.t = x.t[,ind]
x.v = x.v[,ind]

## train lasso 
m.lasso = cv.glmnet(as.matrix(x.t), y.t, family="binomial", pmax=5)
p.t.lasso = predict(m.lasso, newx=as.matrix(x.t)) # predict on training set
p.v.lasso = predict(m.lasso, newx=as.matrix(x.v)) # predict on validation set
cat(paste("Performance of LASSO:", 
          pROC::auc(roc(y.v, p.v.lasso)),
          "\n"))


## train elastic net
m.elastic = cv.glmnet(as.matrix(x.t), y.t, family="binomial", pmax=15, alpha=0.5)
p.t.net = predict(m.elastic, newx=as.matrix(x.t)) # predict on training set
p.v.net = predict(m.elastic, newx=as.matrix(x.v)) # predict on validation set
cat(paste("Performance of Elastic Net:", 
          pROC::auc(roc(y.v, p.v.lasso)),
          "\n"))

save(m.lasso, m.elastic, file="../data/linear_results_filtered.Rdata")

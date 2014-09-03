library("kernlab")
library("glmnet")
library("plyr")
library("dplyr")
library("pROC")

setwd("/gpfs/home/ekramer/Projects/CEC/linear/")

load("~/Projects/CEC/data/cec_filtered2.Rdata")

set.seed(1)

select = dplyr::select


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

## turn down alpha
net = cv.glmnet(as.matrix(x.t), y.t, family="binomial",  pmax=20, alpha=0.5)
p.t.net = predict(net, newx=as.matrix(x.t)) # predict on training set
p.v.net = predict(net, newx=as.matrix(x.v)) # predict on validation set
a = pROC::auc(roc(y.v, p.v.net))

co = coefficients(net)
co = as.matrix(co)
co = co[co!=0, ]

## try on qpcr
net = cv.glmnet(as.matrix(x.t[c("EDN1", "HBEGF", "NR4A2")]), y.t,
                family="binomial",
                alpha=0)
p.t.net = predict(net, newx=as.matrix(x.t[c("EDN1", "HBEGF", "NR4A2")])) # predict on training set
p.v.net = predict(net, newx=as.matrix(x.v[c("EDN1", "HBEGF", "NR4A2")])) # predict on validation set
a = pROC::auc(roc(y.v, p.v.net))

## testing small set
qpcr = read.delim("./data/qpcr2.txt")
genes$Microarray.ID = as.integer(sapply(strsplit(as.character(genes$Array), "_"), function(x) x[2]))

q.genes = colnames(qpcr)[7:length(colnames(qpcr))]
q.genes = q.genes[q.genes %in% colnames(genes)]

d = merge(qpcr, genes[c("Microarray.ID", q.genes)], by="Microarray.ID")

x = d[ , grepl("x$", colnames(d), perl=T)]
y = d[ , grepl("y$", colnames(d), perl=T)]

colnames(x) = gsub(".x$", "", colnames(x))
colnames(y) = gsub(".y$", "", colnames(y))

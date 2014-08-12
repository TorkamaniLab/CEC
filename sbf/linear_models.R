library("kernlab")
library("glmnet")
library("plyr")
library("dplyr")

setwd("/gpfs/home/ekramer/Projects/CEC/sbf/")

load("~/Projects/CEC/data/cec.Rdata")

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
ind = sapply(x.t, function(z) median(z[y.t == "AMI"]) > median(z[y.t !="AMI"]))
x.t = x.t[,ind]
x.v = x.v[,ind]

## train elastic net
net = cv.glmnet(as.matrix(x.t), y.t, pmax=10, family="binomial", alpha=1)
p.t.net = predict(net, newx=as.matrix(x.t)) # predict on training set
p.v.net = predict(net, newx=as.matrix(x.v)) # predict on validation set

## find selected genes
co = coefficients(net)
gene.names = rownames(co)[co[,1] != 0] # find non-zero coefficients
gene.names = gene.names[grep("Intercept", gene.names, invert=T)] # delete intercept term

## train gaussian process
gp = gausspr(x.t[gene.names], y.t)
p.t.gp = predict(gp, newdata=x.t[gene.names], type="probabilities")
p.v.gp = predict(gp, newdata=x.v[gene.names], type="probabilities")

save(net, p.t.net, p.v.net, gp, p.t.gp, p.v.gp, gene.names,
     file="../data/linear_results.Rdata")

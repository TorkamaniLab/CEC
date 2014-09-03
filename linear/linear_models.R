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

## train lasso 
m.lasso = cv.glmnet(as.matrix(x.t), y.t, family="binomial", pmax=5)
p.t.lasso = predict(m.lasso, newx=as.matrix(x.t)) # predict on training set
p.v.lasso = predict(m.lasso, newx=as.matrix(x.v)) # predict on validation set
pROC::auc(roc(y.v, p.v.lasso))

## find selected genes
co = coefficients(m.lasso)
genes.lasso = rownames(co)[co[,1] != 0] # find non-zero coefficients
genes.lasso = genes.lasso[grep("Intercept", genes.lasso, invert=T)] # delete intercept term

## train elastic net
m.elastic = cv.glmnet(as.matrix(x.t), y.t, family="binomial", pmax=20, alpha=0.5)
p.t.net = predict(m.elastic, newx=as.matrix(x.t)) # predict on training set
p.v.net = predict(m.elastic, newx=as.matrix(x.v)) # predict on validation set
pROC::auc(roc(y.v, p.v.net))

co = coefficients(m.elastic)
genes.elastic = rownames(co)[co[,1] != 0] # find non-zero coefficients
genes.elastic = genes.elastic[grep("Intercept", genes.elastic, invert=T)] # delete intercept term

## post-hoc fold changes
qpcr = read.delim("../data/qpcr.txt")

fc.lasso = data.frame(Gene=genes.lasso) %>%
  mutate(Coefficient=-co[genes.lasso,1]) %>%
  mutate(log.FC.discovery=sapply(x.t[genes.lasso], function(z) median(z[y.t=="AMI"]) - median(z[y.t!="AMI"]))) %>%
  mutate(log.FC.validation=sapply(x.v[genes.lasso], function(z) median(z[y.v=="AMI"]) - median(z[y.v!="AMI"]))) %>%
  mutate(FC.discovery=2^log.FC.discovery) %>%
  mutate(FC.validation=2^log.FC.validation) %>%
  arrange(-Coefficient)


fc.elastic = data.frame(Gene=genes.elastic) %>%
  mutate(Coefficient=-co[genes.elastic,1]) %>%
  mutate(log.FC.discovery=sapply(x.t[genes.elastic], function(z) median(z[y.t=="AMI"]) - median(z[y.t!="AMI"]))) %>%
  mutate(log.FC.validation=sapply(x.v[genes.elastic], function(z) median(z[y.v=="AMI"]) - median(z[y.v!="AMI"]))) %>%
  mutate(fc.discovery=2^log.FC.discovery) %>%
  mutate(fc.validation=2^log.FC.validation) %>%
  arrange(-Coefficient) %>%

fc.elastic$Previous.qPCR = fc.elastic$gene %in% colnames(qpcr)


save(m.lasso, m.elastic, file="../data/linear_results_filtered.Rdata")

write.table(fc.lasso, file="../figures/5_gene_model.xls", sep="\t",
            quote=F, row.names=F)

write.table(fc.elastic, file="../figures/18_gene_model.xls", sep="\t",
            quote=F, row.names=F)
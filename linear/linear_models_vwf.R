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
ind = sapply(x.t, function(z) median(z[y.t == "AMI"]) > median(z[y.t !="AMI"]))
x.t = x.t[,ind]
x.v = x.v[,ind]

## train initial model
m = glm(y.t ~ VWF, data=x.t, family="binomial")
o.t = predict(m, type="response")
o.v = predict(m, newdata=x.v, type="response")

## train elastic net
net = cv.glmnet(as.matrix(x.t), y.t, family="binomial", alpha=1, offset=o.t)
p.t.net = predict(net, newx=as.matrix(x.t), offset=o.t) # predict on training set
p.v.net = predict(net, newx=as.matrix(x.v), offset=o.v) # predict on validation set
pROC::auc(roc(y.v, p.v.net))

## find selected genes
co = coefficients(net)
gene.names = rownames(co)[co[,1] != 0] # find non-zero coefficients
gene.names = gene.names[grep("Intercept", gene.names, invert=T)] # delete intercept term


## post-hoc fold changes
fc = data.frame(Gene=gene.names) %>%
  mutate(Coefficient=-co[gene.names,1]) %>%
  mutate(log.FC.discovery=sapply(x.t[gene.names], function(z) median(z[y.t=="AMI"]) - median(z[y.t!="AMI"]))) %>%
  mutate(log.FC.validation=sapply(x.v[gene.names], function(z) median(z[y.v=="AMI"]) - median(z[y.v!="AMI"]))) %>%
  mutate(FC.discovery=2^log.FC.discovery) %>%
  mutate(FC.validation=2^log.FC.validation) %>%
  arrange(-Coefficient)

write.table(fc, file="../figures/coefficients.xls", sep="\t",
            quote=F, row.names=F)

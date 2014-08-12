load("~/Projects/CEC/data/cec.Rdata")

set.seed(1)

select = dplyr::select

x.t = genes[genes$Cohort != "VALIDATION", ] %>%
  select(-Cohort, -SN, -Status, -Array)

y.t = genes[genes$Cohort != "VALIDATION", ]$Status

ind = sapply(x.t, function(z) median(z[y == "AMI"]) > median(z[y!="AMI"]))

x.v = genes[genes$Cohort == "VALIDATION" & genes$SN > 100, ] %>%
  select(-Cohort, -SN, -Status, -Array)
y.v = genes[genes$Cohort == "VALIDATION" & genes$SN > 100, ]$Status

x.t = x.t[,ind]
x.v = x.v[,ind]

m = cv.glmnet(as.matrix(x.t), y.t, pmax=10, family="binomial", alpha=1)

p = predict(m, newx=as.matrix(x.v))
pROC::auc(roc(y.v, p))

co = coefficients(m)
co = as.matrix(co)
gene.names = rownames(co)[co != 0]

plot(genes[gene.names[2:length(gene.names)]],
     col=genes$Status)


library("ROCR")
library("RColorBrewer")
library("glmnet")

pdf("../figures/roc_curves.pdf")
par(cex=1.3)

col = brewer.pal(5, "Dark2")

## Discovery
load("~/Projects/CEC/data/cec_filtered2.Rdata")
load("../data/linear_results_filtered.Rdata")

x.t = genes[genes$Cohort != "VALIDATION", ] %>%
  select(-Cohort, -SN, -Status, -Array)
y.t = genes[genes$Cohort != "VALIDATION", ]$Status

x.v = genes[genes$Cohort == "VALIDATION" & genes$SN > 100, ] %>%
  select(-Cohort, -SN, -Status, -Array)
y.v = genes[genes$Cohort == "VALIDATION" & genes$SN > 100, ]$Status

ind = sapply(x.t, function(z) median(z[y.t == "AMI"]) > 1 + median(z[y.t !="AMI"]))
x.t = x.t[,ind]
x.v = x.v[,ind]

p.t = predict(m.elastic, newx=as.matrix(x.t))
p.v = predict(m.elastic, newx=as.matrix(x.v))

roc.m.t = performance(prediction(p.t, y.t), "tpr", "fpr")
roc.m.v = performance(prediction(p.v, y.v), "tpr", "fpr")

plot(roc.m.t, 
     main="Discovery Set",
     col=col[1],
     lwd=2)
legend(0.2, 0.5, c("Microarray: 18 Gene Model"), col=col[1], lwd=2)

plot(roc.m.v,
     main="Validation Set",
     col=col[2],
     lwd=2)
legend(0.2, 0.5, c("Microarray: 18 Gene Model"), col=col[2], lwd=2)


## Validation


# qpcr 
qpcr = read.delim("../data/qpcr.txt")

qpcr.genes = df %>%
  filter(sig) %>%
  filter(previous.qpcr) %>%
  select(gene) 

qpcr2 = filter(qpcr, Used.in.Validation != "NO")

y = qpcr$Condition
x = qpcr[qpcr.genes$gene] %>% as.matrix
x = apply(x, 2, function(y) y - qpcr$GABPB1)

m.qpcr = cv.glmnet(x, 
                   y, 
                   nfolds=nrow(x), 
                   keep=T, 
                   family="binomial",
                   grouped=F,
                   alpha=0.5)

p.qpcr = m.qpcr$fit.preval[ , m.qpcr$lambda == m.qpcr$lambda.1se]

roc.qpcr = performance(prediction(p.qpcr, y), "tpr", "fpr")


# cec data
cec = read.delim("../data/cec_counts.txt")
roc.cec = performance(prediction(cec$Count, cec$Status), "tpr", "fpr")


plot(roc.m.v,
     main="Validation Set",
     col=col[2],
     lwd=2.5)
plot(roc.qpcr, col=col[3], lwd=2, add=T)
plot(roc.cec, col=col[5], lwd=2, add=T)

legend(0.2, 
       0.5, 
       c("Microarray: 18 Gene Model",
                   "qPCR: 4 Gene Model",
                   "CEC count"), 
       col=col[c(2,3,5)], lwd=2.4)
dev.off()

x = -x 

pdf("../figures/qpcr_values.pdf")

par(cex=2)

par(mfrow=c(2,2))


lapply(1:4,
       function(i) stripchart(x[,i] ~ y,
                              vertical=T,
                              method="jitter",
                              pch=20,
                              frame=F,
                              ylab='-ddCT',
                              col=brewer.pal(3, "Set2")[1:2],
                              main=colnames(x)[i],
                              cex=1.3))

lapply(1:4,
       function(i) stripchart(2^x[,i] ~ y,
                              vertical=T,
                              method="jitter",
                              pch=20,
                              frame=F,
                              ylab=parse(text='2^-ddCT'),
                              col=brewer.pal(3, "Set2")[1:2],
                              main=colnames(x)[i],
                              cex=1.3))


dev.off()


x.t.q = x.t[,qpcr.genes$gene]
x.v.q= x.v[qpcr.genes$gene]
mmq = cv.glmnet(as.matrix(x.t.q), y.t, family="binomial", alpha=0.5)

p.x.v.q = predict(mmq, newx=as.matrix(x.v.q))


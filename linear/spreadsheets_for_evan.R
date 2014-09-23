library("dplyr")
library("glmnet")

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

d.t = genes %>%
  filter(Cohort != "VALIDATION") %>%
  select(Array, Status) %>%
  mutate(Cohort="Discovery") %>%
  mutate(logOR=as.numeric(p.t))
  
d.v = genes %>%
  filter(Cohort == "VALIDATION") %>%
  filter(SN > 100) %>%
  select(Array, Status) %>%
  mutate(Cohort="Validation") %>%
  mutate(logOR=as.numeric(p.v))

d.n = genes %>%
  filter(SN < 100) %>%
  select(Array, Status) %>%
  mutate(Cohort="DISCARDED (LOW SN)") %>%
  mutate(logOR=NA)

d = rbind(d.t, d.v, d.n) %>%
  mutate(Patient=sapply(strsplit(as.character(Array), "_"), function(x) x[2])) %>%
  select(Patient, Array, Status, Cohort, logOR) %>%
  mutate(logOR=-logOR)

## qpcr

qpcr = read.delim("../data/qpcr.txt")

z = qpcr %>%
  select(-Process, -Sample.ID, -Condition, -Used.in.Validation, -Run.ID) %>%
  as.matrix

qpcr.genes = df %>%
  filter(sig) %>%
  filter(previous.qpcr) %>%
  select(gene) 

qpcr2 = qpcr %>%
  filter(!(Sample.ID %in% c("P1-2", "P3-2", "D30-2", "D5-2")))

y = qpcr2$Condition
x = qpcr2[qpcr.genes$gene] %>% as.matrix
x = apply(x, 2, function(y) y - qpcr2$GABPB1)

m.qpcr = cv.glmnet(x, 
                   y, 
                   nfolds=nrow(x), 
                   keep=T, 
                   family="binomial",
                   grouped=F,
                   alpha=0.5)

p.qpcr = logit(m.qpcr$fit.preval[ , m.qpcr$lambda == m.qpcr$lambda.1se])

d.q = qpcr2 %>% 
  select(Sample.ID, Condition) %>%
  mutate(logOR=p.qpcr) %>%
  mutate(Discarded=F) 

d.dq = qpcr %>%
  filter(Sample.ID %in% c("P1-2", "P3-2", "D30-2", "D5-2")) %>%
  select(Sample.ID, Condition) %>%
  mutate(logOR=NA, Discarded=T)
  
d2 = rbind(d.q, d.dq)

write.table(d,
            quote=F,
            row.names=F,
            sep="\t", 
            file="../figures/microarray_samples.xls")

write.table(d2,
            quote=F,
            row.names=F,
            sep="\t",
            file="../figures/qpcr_samples.xls")

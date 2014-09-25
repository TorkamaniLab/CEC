library("limma")
library("plyr")
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

lm.t = lmFit(t(x.t), model.matrix( ~ y.t))
lm.v = lmFit(t(x.v), model.matrix( ~ y.v))

elm.t = eBayes(lm.t)
elm.v = eBayes(lm.v)

df.t = topTable(elm.t, coef=2, number=ncol(x.t))
df.v = topTable(elm.v, coef=2, number=ncol(x.v))

p.vals = inner_join(
  df.t %>% mutate(gene=row.names(df.t)),
  df.v %>% mutate(gene=row.names(df.v)),
  by="gene"
  )
colnames(p.vals) = gsub(".y$", ".validation", colnames(p.vals), perl=T)
colnames(p.vals) = gsub(".x$", ".discovery", colnames(p.vals), perl=T)


coefficients = m.elastic %>%
  coef %>%
  as.matrix %>%
  as.data.frame
coefficients$gene = row.names(coefficients)
colnames(coefficients) = c("coef", "gene")

df = inner_join(p.vals, coefficients) %>%
  filter(coef != 0) %>%
  arrange(coef) %>%
  select(gene, coef, logFC.discovery, P.Value.discovery, adj.P.Val.discovery,
         logFC.validation, P.Value.validation, adj.P.Val.validation) %>%
  mutate(coef=-coef) %>%
  mutate(logFC.validation=-logFC.validation) %>%
  mutate(FC.validation=2^logFC.validation) %>%
  mutate(logFC.discovery=-logFC.discovery) %>%
  mutate(FC.validation=2^logFC.validation) %>%
  mutate(previous.qpcr=gene %in% colnames(qpcr)) %>%
  mutate(sig= adj.P.Val.discovery < 1e-3 & adj.P.Val.validation < 1e-3)

write.table(df, quote=F, row.names=F, col.names=T, sep="\t",
            file="~/Projects/CEC/figures/11_gene_model.xls")




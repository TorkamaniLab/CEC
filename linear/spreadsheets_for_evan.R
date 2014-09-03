library("pROC")
library("dplyr")

load("~/Projects/CEC/data/cec_filtered2.Rdata")

evan.data = read.delim("./data/evan_data.txt", stringsAsFactors=F)

## sample xls
samples = genes %>%
  select(Status, Array, Cohort, SN) %>%
  mutate(New.Cohort=ifelse(Cohort=="VALIDATION", "VALIDATION", "DISCOVERY")) %>%
  mutate(keep=SN>100 | is.na(SN)) %>%
  mutate(ID=sapply(strsplit(as.character(genes$Array), "_"), function(x) x[2])) %>%
  select(ID, Array, Status, Original.Cohort=Cohort, New.Cohort, Keep=keep)
write.table(samples, sep="\t", quote=F, row.names=F, file="./data/sample_info.xls")

## gene xls
d = genes[genes$Cohort != "VALIDATION", ]
d1 = genes[genes$Cohort == "TRAINING", ]
d2 = genes[genes$Cohort == "TESTING", ]

y = d$Status
y1 = d1$Status
y2 = d2$Statuss

genes = gsub(" (.*)", "", evan.data$Gene)
genes = genes[genes %in% colnames(d)]

gene.info = data.frame(Gene=genes,
                       AUC1=sapply(genes, function(g) auc(roc(y1, d1[[g]]))),
                       AUC2=sapply(genes, function(g) auc(roc(y2, d2[[g]]))),
                       AUC.overall=sapply(genes, function(g) auc(roc(y, d[[g]]))))

write.table(gene.info, sep="\t", quote=F, row.names=F, file="./data/gene_info.xls")


## correlation matrix
s = cor(d[c("HBEGF", "NR4A2", "NR4A3", "SYTL3", "EDN1")])
s = as.data.frame(s)
s$Gene = row.names(s)
s = s[c("Gene", "HBEGF", "NR4A2", "NR4A3", "SYTL3", "EDN1")]

write.table(s, sep="\t", quote=F, row.names=F, file="./data/correlations.xls")

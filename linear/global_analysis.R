setwd("/gpfs/home/ekramer/Projects/CEC/linear/")

load("~/Projects/CEC/data/cec.Rdata")

## fold change
genes = genes[genes$SN > 100 | is.na(genes$SN), ]

y = genes$Status
x = select(genes, -Status, -Array, -Cohort, -SN)

fc = sapply(x, function(z) median(z[y=="AMI"]) - median(z[y!="AMI"]))

## pca
p = prcomp(as.matrix(x), center=T, scale.=T)
s = p$rotation[,3]

fc.rnk = data.frame(Gene=names(fc), fc=fc)
s.rnk = data.frame(Gene=names(s), s=s)


write.rnk = function(x, f) write.table(x, quote=F, sep="\t", row.names=F, col.names=F, file=f)
write.rnk(fc.rnk, "../data/GSEA/fc.rnk")
write.rnk(s.rnk, "../data/GSEA/s.rnk")

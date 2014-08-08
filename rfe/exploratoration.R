library("dplyr")

load("/gpfs/home/ekramer/Projects/CEC/data/exprs.Rdata")

d = merge(gene.exprs, pheno[c("PTID", "status")])

x = as.matrix(select(d, -status, -PTID))

p = prcomp(x, center=T, scale.=T)
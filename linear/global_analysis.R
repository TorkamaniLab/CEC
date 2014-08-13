setwd("/gpfs/home/ekramer/Projects/CEC/linear/")

load("~/Projects/CEC/data/cec.Rdata")

## fold change
genes = genes[genes$SN > 100 | is.na(genes$SN), ]

y = genes$Status
x = select(genes, -Status, -Array, -Cohort, -SN)

fc = sapply(x, function(z) median(z[y=="AMI"]) - median(z[y!="AMI"]))

## pca
p = prcomp(as.matrix(x), center=T, scale.=T)
s = p$rotation[,1:3]

fc.rnk = data.frame(Gene=names(fc), fc=fc)



write.rnk = function(x, f) write.table(x, quote=F, sep="\t", row.names=F, col.names=F, file=f)
write.rnk(fc.rnk, "../data/GSEA/fc.rnk")
lapply(1:3, function(i){
  s.rnk = data.frame(Gene=row.names(s), s=s[,1])
  write.rnk(s.rnk, paste("../data/GSEA/s", i, ".rnk", sep=""))
})

load("~/Projects/CEC/data/cec.Rdata")
ind = laply(genes[gene.names], 
      function(g) which(laply(probes, function(p) all(g == p))), 
      .progress="text")
probe.names = colnames(probes)[ind]

library("gplots")
library("RColorBrewer")

hm.simple <- function(d, main, f.names, cexCol=1){
  x = as.matrix(d[f.names])
  x = t(scale(x))
  
  heatmap.2(x,
            Rowv=T,
            col=redgreen(75), 
            scale="row", 
            trace="none", 
            key=F, 
            lwid=c(.1,1), 
            cexCol=cexCol, 
            cexRow=1,
            main=main,
            labCol=as.character(d$Status))
  
}

load("../data/rfe_data.Rdata")
load("../data/cec.Rdata")

d = genes[genes$SN > 100 | is.na(genes$SN), ]

validation = d[d$Cohort == "VALIDATION", ]
training = d[d$Cohort != "VALIDATION", ]
f.names = gene.rfe$denovo$optVariables


pdf("../figures/heatmaps.pdf", width=10)

hm.simple(training, "Training Set", f.names)
hm.simple(validation, "Validation Set", f.names)

dev.off()

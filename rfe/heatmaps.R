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

load("../data/linear_results.Rdata")
load("../data/cec.Rdata")

d = genes[genes$SN > 100 | is.na(genes$SN), ]

validation = d[d$Cohort == "VALIDATION", ]
training = d[d$Cohort != "VALIDATION", ]

pdf("../figures/heatmaps.pdf", width=10)

hm.simple(training, "Discovery Set", gene.names)
hm.simple(validation, "Validation Set", gene.names)

dev.off()

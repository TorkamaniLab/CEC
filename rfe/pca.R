library("pcaMethods")
library("plyr")
library("dplyr")
library("RColorBrewer")

setwd("/gpfs/home/ekramer/Projects/CEC/rfe/")

load("/gpfs/home/ekramer/Projects/CEC/data/cec.Rdata")

tmp = genes[genes$SN > 100 | is.na(genes$SN), ]

x = select(tmp, -Status, -Array, -Cohort, -SN)
p = prcomp(as.matrix(x), center=T, scale.=T)

## wilcox tests
s = as.data.frame(p$x[,1:4])
lapply(s, function(z) wilcox.test(z ~ tmp$Status))

## plotting
pdf("../figures/pca.pdf", width=12)

par(mar=c(5.1, 5.1, 4.1, 2.1))
par(mfrow=c(1,3))

col = brewer.pal(3, "Set2")
names(col) = levels(tmp$Status)

make.stripchart <- function(i){
  z = s[,i]
  stripchart(z ~ tmp$Status, 
             vertical=T,
             pch=20,
             frame=F,
             method="jitter",
             cex=1.7,
             col=col,
             main=paste("Principle Component", i),
             ylab="Score",
             cex.axis=1.5,
             cex.main=1.7,
             cex.lab=1.5)
  
}

lapply(1:3, make.stripchart)


par(mfrow=c(1,1))

barplot(summary(p)$importance[2, 1:10], space=0.1,
        main="Variance Explained by Principle component",
        col=col[3])

plot(s[,1:3], col=col[as.numeric(tmp$Status)], pch=20, main="Principle Component Scatter Plots")

dev.off()




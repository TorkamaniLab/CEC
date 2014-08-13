load("../data/cec.Rdata")
load("../data/linear_results.Rdata")

d = genes[genes$SN > 100 | is.na(genes$SN), ]
x = d[gene.names]
y = d$Status

r = as.data.frame(lapply(x, function(z) residuals(lm(z ~ y))))


lapply(x, function(z) wilcox.test(z ~ y))

par(mfrow=c(1,1))

pdf("../figures/scatterplots.pdf")

par(cex=1.3)

plot(x[c("HBEGF", "SYTL3")], 
     bty="n",
     pch=20,
     main="Expression Values",
     col=brewer.pal(3, "Set2")[1:2][as.numeric(y)])

plot(r[c("HBEGF", "SYTL3")], 
     bty="n",
     pch=20,
     main="Controlling for Patient's Status",
     col=brewer.pal(3, "Set2")[1:2][as.numeric(y)])

plot(x[c("CXCL16", "CD55")], 
     bty="n",
     pch=20,
     main="Expression Values",
     col=brewer.pal(3, "Set2")[1:2][as.numeric(y)])


plot(r[c("CXCL16", "CD55")], 
     bty="n",
     pch=20,
     main="Controlling for Patient's Status",
     col=brewer.pal(3, "Set2")[1:2][as.numeric(y)])

plot(x, 
     pch=20,
     col=brewer.pal(3, "Set2")[1:2][as.numeric(y)],
     main="Expression Values")

plot(r,
     pch=20,
     col=brewer.pal(3, "Set2")[1:2][as.numeric(y)],
     main="Controlling for Patient's Status")

dev.off()


library("gplots")

hm.simple = function(x, y, main="", Rowv=NULL){
  heatmap.2(t(as.matrix(x)),
            scale="row",
            col=redgreen(75),
            Rowv=Rowv,
            trace="none",
            labCol=y,
            key=F,
            dist=function(z) dist(z, method="manhattan"),
            main=main)
  
}

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

co = m.elastic %>%
  coef %>%
  as.matrix
co = co[co!=0]

m.genes = names(co)

x.t = scale(x.t)
x.v = scale(x.v)

x.t = x.t[ , colnames(x.t) %in% m.genes]
x.v = x.v[ , colnames(x.v) %in% m.genes]

x = scale(rbind(x.t, x.v))

dend = x %>%
  t %>%
  dist %>%
  hclust %>%
  as.dendrogram

# order the genes by their coefficients
o = rank(co)[colnames(x)]

pdf("../figures/heatmaps.pdf")
par(cex=1.3)
hm.simple(x.t, y.t, main="Discovery", Rowv=F)
hm.simple(x.v, y.v, main="Validation", Rowv=F)
dev.off()


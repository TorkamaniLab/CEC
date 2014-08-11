
hm.simple <- function(d, main, f.names, cexCol=1){
  x = as.matrix(d[f.names])
  x = t(scale(x))
  
  heatmap.2(x,
            Rowv=T,
            col=col, 
            scale="row", 
            trace="none", 
            key=F, 
            lwid=c(.1,1), 
            cexCol=cexCol, 
            cexRow=1,
            main=main,
            labCol=as.character(d$Status),
            dist=function(x) dist(x, method="manhattan"))
  
}

load("../data/rfe_data.Rdata")

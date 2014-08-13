library("ROCR")
library("caret")
library("doMC")
library("dplyr")
library("pROC")

roc.plot <- function(p1, p2, y, ...){
  
    ## Plots ROC curves for two sets of predictions
    ## p1 -- "Linear model"
    ## p2 -- "Gaussian Process"
    ## y -- common set of labels for predictions
  
    require("ROCR")
    require("pROC")
    require("RColorBrewer")
    
    roc1 = p1 %>% prediction(y) %>% performance("tpr", "fpr")
    roc2 = p2 %>% prediction(y) %>% performance("tpr", "fpr")
    
    auc1 = y %>% roc(p1) %>% auc
    auc2 = y %>% roc(p2) %>% auc
    
    col = brewer.pal(3, "Accent")[1:2]
    
    plot(roc1, 
         col=col[1], 
         lwd=2, ...)
    
    plot(roc2, 
         col=col[2],
         lwd=2,
         add=T)
    
    legend.text = paste(c("Gaussian Process", "Linear Model"),
                         " (",
                         round(c(auc1, auc2), 2),
                         ")", 
                        sep="")
    
    legend(0.2, 
           0.5, 
          legend.text, 
          lwd=2, 
          col=col)
}



registerDoMC(10)
setwd("/gpfs/home/ekramer/Projects/CEC/linear")
load("../data/linear_results.Rdata")
load("../data/cec.Rdata")

y.train = genes$Status[genes$Cohort != "VALIDATION"]
y.validation = genes$Status[genes$Cohort == "VALIDATION" & genes$SN > 100]

pdf("../figures/roc_curves.pdf")
par(cex=1.3)
roc.plot(p.t.gp[,2], p.t.net, y.train, main="Discovery Set")
roc.plot(p.v.gp[,2], p.v.net, y.validation, main="Validation Set")
dev.off()


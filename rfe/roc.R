library("ROCR")
library("caret")
library("doMC")
library("dplyr")
library("pROC")

roc.plot <- function(p1, p2, y, ...){
  
    ## Plots ROC curves for two sets of predictions
    ## p1 -- "All probes"
    ## p2 -- "10 probe model"
    ## y -- common set of labels for predictions
  
    require("ROCR")
    require("pROC")
    require("RColorBrewer")
    
    roc1 = p1[,2] %>% prediction(y) %>% performance("tpr", "fpr")
    roc2 = p2[,2] %>% prediction(y) %>% performance("tpr", "fpr")
    
    auc1 = y %>% roc(p1[,2]) %>% auc
    auc2 = y %>% roc(p2[,2]) %>% auc
    
    col = brewer.pal(3, "Accent")[1:2]
    
    plot(roc1, 
         col=col[1], 
         lwd=2, ...)
    
    plot(roc2, 
         col=col[2],
         lwd=2,
         add=T)
    
    legend.text = paste(c("Linear Model", "Gaussian Process"),
                         " (",
                         round(c(auc1, auc2), 2),
                         ")", 
                        sep="")
    
    legend(0.5, 
           0.5, 
          legend.text, 
          lwd=2, 
          col=col)
}



registerDoMC(10)
setwd("/gpfs/home/ekramer/Projects/CEC/sbt")
load("../data/rfe_data.Rdata")
load("../data/cec.Rdata")

genes = genes[genes$SN > 100 | is.na(genes$SN), ]
probes = probes[probes$SN > 100 | is.na(probes$SN), ]

y = genes$Status
ten.probes = genes[c("Cohort", gene.rfe$denovo$optVariables)]
all.probes= select(probes, -Status, -Array, -SN)
  
## RUN LOOCV AND PREDICT ON VALIDATION SET
p.ten = get.predictions(ten.probes, y, ntree=1000)
p.all = get.predictions(all.probes, y, ntree=1000)

## PLOT ROC CURVES
y.train = y[genes$Cohort != "VALIDATION"]
y.validation = y[genes$Cohort == "VALIDATION"]

pdf("roc_curves.pdf")
roc.plot(p.all$p.train, p.ten$p.train, y.train, main="Discovery Set")
roc.plot(p.all$p.validation, p.ten$p.validation, y.validation, main="Validation Set")
dev.off()


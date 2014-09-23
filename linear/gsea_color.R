library("RColorBrewer")
library("plyr")
library("dplyr")

col = c("#E0AEBB",
       "#C2F195",
       "#88D9D8",
       "#DDBC6A",
       "#E7E2B7",
       "#A5CDF0",
       "#ECA982",
       "#89E0A7",
       "#BDC1BD",
       "#DABEE5",
       "#B4D59A",
       "#EBEC95",
       "#D0EBDE",
       "#91D3B3",
       "#D7C689",
       "#C1BB97",
       "#EDCFB5",
       "#C1D173",
       "#C1D2E2",
       "#DCB1A8",
       "#CAF2BF",
       "#F1A3B8",
       "#EECBD7",
       "#A0C681")

gsea.dir = "/gpfs/home/ekramer/Projects/CEC/data/GSEA/sep03/fc.rnk_ReactomePathways.gmt.GseaPreranked.1409788773303/"

pos = read.delim(paste(gsea.dir, 
            "gsea_report_for_na_pos_1409788773303.xls", 
            sep=""))

neg = read.delim(paste(gsea.dir, 
                       "gsea_report_for_na_neg_1409788773303.xls", 
                       sep=""))




d = rbind(pos, neg)

load("~/Projects/CEC_old/CEC2/figures/general/reactome.types.Rdata")
reactome.types$NAME = toupper(as.character(reactome.types$displayName))

d2 = merge(d, reactome.types, all.x=T)
d2 = d2[order(d2$NES, decreasing=T), ]
d2 = d2[!is.na(d2$parent), ]
d2 = d2[d2$SIZE < 1000, ]

pdf("../figures/GSEA.pdf")

make.plot <- function(i){

set.seed(i)
x = runif(nrow(d2))
r = sqrt(d2$SIZE)/400

## without labels
# symbols(x, 
#         d2$NES, 
#         circles=r,
#         bty="n",
#         inches=F,
#         axes=F,
#         xlab="",
#         ylab="Normalized Enrichment Score",
#         main="Gene Enrichment Analysis",
#         bg=col[as.numeric(factor(d2$parent))],
#         fg="white",
#         xlim=c(0,1.5))
# axis(2)
# 
# n = unique(d2$parent)
# n = gsub(" by Scavenger Receptors", "", n)
# n = gsub(" of small molecules", "", n)
# n = gsub(" and maintenance", "", n)
# 
# legend(1.05, 5.7, 
#        n, 
#        col=col[as.numeric(factor(unique(d2$parent)))],
#        pch=20,
#        cex=0.65,
#        ncol=1,
#        y.intersp=1,
#        pt.cex=1,
#        bty="n")


## with labels
symbols(x, 
        d2$NES, 
        circles=r,
        bty="n",
        inches=F,
        axes=F,
        xlab="",
        ylab="Normalized Enrichment Score",
        main="Gene Enrichment Analysis",
        bg=col[as.numeric(factor(d2$parent))],
        fg="white",
        xlim=c(0,1.5))
axis(2)

n = unique(d2$parent)
n = gsub(" by Scavenger Receptors", "", n)
n = gsub(" of small molecules", "", n)
n = gsub(" and maintenance", "", n)

legend(1.05, 4.7, 
       n, 
       col=col[as.numeric(factor(unique(d2$parent)))],
       pch=20,
       cex=0.65,
       ncol=1,
       y.intersp=1,
       pt.cex=1,
       bty="n")

labels = c("GPCR ligand\nbinding",
           "A1 Rhodopsin like\nreceptors",
           "Neuronal System",
           "Hemostasis",
           "Platelet\nActivation")
text(x[1:5], d2$NES[1:5], labels, cex=0.6)

ind = nrow(d2)
labels2 = c("Generic Transcription\nPathway")

text(x[ind], d2$NES[ind], labels2, cex=0.6)

}

#lapply(1:5, make.plot)
make.plot(5)

dev.off()
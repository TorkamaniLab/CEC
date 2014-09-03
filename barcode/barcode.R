library("plyr")
library("reshape2")
library("dplyr")
library("affy")
library("frma")
library("data.table")

### Barcoding our samples

celfile.path = "/gpfs/home/ekramer/Projects/CEC/data/CEL"

bc.cec = ReadAffy(celfile.path=celfile.path) %>%
  frma %>%
  barcode %>%
  as.data.frame
bc.cec$AffyID = row.names(bc.cec)

### Processing Catalog

cat = get(load("catalogGPL570-celltype_v3.rda"))

bc.tissue = cat %>% 
  group_by(AffyID) %>%
  do(data.frame(AffyID=.$AffyID,
                CellTypes=unlist(strsplit(.$CellTypes, ";", fixed=T)),
                stringsAsFactors=F)) %>%
  ungroup %>%
  filter(!(CellTypes %in% c("Many", "None"))) %>%
  as.data.table %>% 
  dcast(AffyID ~ CellTypes, function(...) 1, fill=0) %>%
  as.data.frame
row.names(bc.tissue) = bc.tissue$AffyID

### Joining and Clustering

# bc= bc.cec %>%
#   merge(bc.tissue) %>%
#   select(-AffyID)
# 
# d = dist(bc %>% as.matrix %>% t, method="binary")
# 
# closest = d %>%
#   as.matrix %>%
#   melt %>%
#   filter(value != 0) %>%
#   filter(grepl("CEL", Var1)) %>%
#   group_by(Var1) %>%
#   filter(value %in% sort(value)[1:3]) %>%
#   ungroup %>%
#   arrange(Var1, value)
# 
# pdf("clustering.pdf", width=10)
# plot(hclust(d), cex=0.2)
# dev.off()
# 
# counts = colSums(bc)
# pdf("total_gene_count.pdf")
# stripchart(counts ~ grepl("CEL", colnames(bc)), vertical=T, method="jitter",
#            pch=20, frame=F, ylab="Number of Genes Expressed",
#            main="Total Gene Counts",
#            group.names=c("Gene Barcode Database", "CEC samples"))
# dev.off()

tissue.pairs = bc.tissue %>% 
  select(-AffyID) %>% 
  colnames %>% 
  combn(2) %>%
  t %>%
  as.data.frame(stringsAsFactors=F)

tmp = mlply(tissue.pairs, 
            function(V1, V2) bc.tissue[[V1]] | bc.tissue[[V2]],
            .progress="text")

bc.tissue.pairs = as.data.frame(tmp) %>%
  mutate(AffyID=bc.tissue$AffyID) %>%
  arrange(AffyID)

p = bc.tissue.pairs[,1]

bc.cec = merge(bc.cec, select(bc.tissue, AffyID)) %>%
  arrange(AffyID)

bc.tissue.pairs = bc.tissue.pairs %>% arrange(AffyID)

library("affy")
library("reshape2")
library("data.table")
library("hgu133plus2.db")
library("plyr")
library("dplyr")
library("entropy")

select = dplyr::select

get.pheno <- function(){
  pheno.dir = "/gpfs/home/ekramer/Projects/CEC/data/CEC_MI_TRANSFER/"
  
  training.file = paste(pheno.dir, "Scripps CEC sample info_cofact.txt", sep="")
  testing.file = paste(pheno.dir, "Scripps CEC sample info_testing.txt", sep="")
  validation.file = paste(pheno.dir, "Scripps CEC sample info_verification_cofact.txt", sep="")
  
  training = read.table(training.file, sep="\t", header=T, comment.char="#")
  testing = read.delim(testing.file, sep="\t", header=T, comment.char="#")
  validation = read.delim(validation.file, sep="\t", header=T, comment.char="#")
  
  training$Cohort = "TRAINING"
  testing$Cohort = "TESTING"
  validation$Cohort = "VALIDATION"
  
  training %>%
    merge(testing, all=T) %>%
    merge(validation, all=T) %>%
    mutate(Status=factor(ifelse(group=="Donor", "Ctrl", "AMI"))) %>%
    mutate(Array=gsub("PLUS2", "PLUS_2", id)) %>%
    mutate(Cohort=factor(Cohort)) %>% 
    select(Array, PID=pid, Cohort, Status, SN=sn)
}

probes2genes <- function(probes, filter.func=IQR){
  
  iqrs = data.frame(probe_id=colnames(probes)) %>% # creates a data.frame of probe_ids
    mutate(iqr=apply(probes, 2, filter.func)) %>% # calculates function for each probe
    merge(toTable(hgu133plus2SYMBOL)) %>% # maps probe ids to symbols
    group_by(symbol) %>%
    filter(iqr == max(iqr)) %>% # filters for top probe for each gene
    ungroup
  
  x = probes[ ,as.character(iqrs$probe_id)] # select best probes
  colnames(x) = as.character(iqrs$symbol) # convert column names to genes
  return(x)
}

ent = function(x) x %>%
    cut(5) %>%
    table %>%
    entropy(method="shrink", verbose=F)


## get phenotype information 
pheno = get.pheno()

## set up base directories
celfile.path = "/gpfs/home/ekramer/Projects/CEC/data/CEL"

## read affy file and normalize with RMA
affy = ReadAffy(celfile.path=celfile.path)
eset = rma(affy, bgversion=1)
probes = t(exprs(eset))

## keep probes that could to be unique to endothelial cells
cec.possible <-read.table("~/Projects/CEC/data/CEC_MI_TRANSFER/absent_blood_ngt-igt-t2d.txt",
                          stringsAsFactors=F)$V1
probes = probes[ , colnames(probes) %in% cec.possible]

## map probes to genes
# genes = probes2genes(probes, filter.func=ent)
genes = probes2genes(probes)

## merge with phenotype data
merge.pheno = function(x) x %>%
  as.data.frame %>%
  mutate(Array=gsub(".CEL$", "", row.names(x), perl=T)) %>%
  merge(pheno[c("Array", "Status", "Cohort", "SN")])

probes = merge.pheno(probes)
genes = merge.pheno(genes)

save(probes, genes, pheno, file="/gpfs/home/ekramer/Projects/CEC/data/cec_filtered2.Rdata")

library("affy")
library("reshape2")
library("data.table")
library("hgu133plus2.db")
library("plyr")
library("dplyr")

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

probes2genes <- function(probes){
  
  iqrs = data.frame(probe_id=colnames(probes)) %>% # creates a data.frame of probe_ids
    mutate(iqr=apply(probes, 2, IQR)) %>% # calculates IQR for each probe
    merge(toTable(hgu133plus2SYMBOL)) %>% # maps probe ids to symbols
    group_by(symbol) %>%
    filter(iqr == max(iqr)) %>% # filters for top probe for each gene
    ungroup
  
  x = probes[ ,as.character(iqrs$probe_id)] # select best probes
  colnames(x) = as.character(iqrs$symbol) # convert column names to genes
  return(x)
}

## get phenotype information 
pheno = get.pheno()

## set up base directories
celfile.path = "/gpfs/home/ekramer/Projects/CEC/data/CEL"

## read affy file and normalize with RMA
affy = ReadAffy(celfile.path=celfile.path)
eset = rma(affy, normalize=F)
probes = t(exprs(eset))

## reduce probes to genes
iqrs = data.frame(probe_id=colnames(probes)) %>%
  mutate(iqr=apply(probes, 2, IQR)) %>%
  merge(toTable(hgu133plus2SYMBOL))



library("affy")
library("reshape2")
library("data.table")
library("hgu133plus2.db")
library("plyr")
library("dplyr")

### FUNCTIONS ###

get.probe.exprs <- function(celfile.path){
  
  ## reads in CEL files from celfile.path
  ## returns data.frame of expression values
  
  probes = ReadAffy(celfile.path=celfile.path)
  
  eset = expresso(probes,
                  bgcorrect.method="rma",
                  normalize.method="constant",
                  pmcorrect.method="pmonly",
                  summary.method="avgdiff")
  
  probe.exprs = as.data.frame(t(exprs(eset)))
  probe.exprs$PTID = gsub(".CEL", "", row.names(probe.exprs), perl=T)
  
  return(probe.exprs)
}

get.gene.exprs <- function(probe.exprs, fun.aggregate=median){
  
  ## takes probe-set expression
  ## returns gene-level expression 
  ## agg.func aggregates across probe sets
  
  probe2gene = toTable(hgu133plus2SYMBOL)
  
  tmp = probe.exprs %>%
    as.data.table %>% # converting to data.table makes this go much faster
    melt(id.vars="PTID", variable.name="probe_id") # see reshape2 manual for details of melting/casting
    
  tmp2 = merge(tmp, probe2gene, by="probe_id") %>% # merge with mapping to genes
    select(-probe_id) %>% # delete column with probe_set id
    dcast(PTID~symbol, fun.aggregate=fun.aggregate) %>% # aggregate across probe sets
    as.data.frame
}

get.pheno <- function(){
  pheno.dir = "/gpfs/home/ekramer/Projects/CEC/data/CEC_MI_TRANSFER/"
  
  training.file = paste(pheno.dir, "Scripps CEC sample info_cofact.txt", sep="")
  testing.file = paste(pheno.dir, "Scripps CEC sample info_testing.txt", sep="")
  validation.file = paste(pheno.dir, "Scripps CEC sample info_verification_cofact.txt", sep="")
  
  training = read.table(training.file, sep="\t", header=T, comment.char="#")
  testing = read.delim(testing.file, sep="\t", header=T, comment.char="#")
  validation = read.delim(validation.file, sep="\t", header=T, comment.char="#")
  
  training$cohort = "TRAINING"
  testing$cohort = "TESTING"
  validation$cohort = "VALIDATION"
  
  training %>%
    merge(testing, all=T) %>%
    merge(validation, all=T) %>%
    mutate(status=ifelse(group=="Donor", "Ctrl", "AMI")) %>%
    mutate(PTID=gsub("PLUS2", "PLUS_2", id))
    
}

### SCRIPT ###

probe.exprs = get.probe.exprs(celfile.path)
gene.exprs = get.gene.exprs(probe.exprs)
pheno = get.pheno()

save(probe.exprs, gene.exprs, pheno, file="/gpfs/home/ekramer/Projects/CEC/data/exprs.Rdata")

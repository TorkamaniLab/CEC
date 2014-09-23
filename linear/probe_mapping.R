library("glmnet")
library('dplyr')
library("hgu133plus2.db")

select = dplyr::select

load("~/Projects/CEC/data/cec_filtered2.Rdata")
load("../data/linear_results_filtered.Rdata")

genes = genes %>%
  as.matrix %>%
  t %>%
  as.data.frame
genes$Symbol = row.names(genes)

probes = probes %>%
  as.matrix %>%
  t %>%
  as.data.frame
probes$Probe_set = row.names(probes)

gene2probe = inner_join(genes, probes) %>%
  select(Symbol, Probe_set) %>%
  filter(Symbol != "Array")

gene2ref_seq = hgu133plus2REFSEQ %>%
  as.data.frame %>%
  select(Probe_set=probe_id, Ref_seq=accession) %>%
  inner_join(gene2probe)

m.ref_seq = gene2ref_seq %>%
  filter(Symbol %in% m.genes) %>%
  filter(grepl("^NM", Ref_seq, perl=T)) %>%
  select(Symbol, Probe_set, Ref_seq) %>%
  arrange(Symbol)


write.table(m.ref_seq,
            row.names=F,
            col.names=T,
            sep="\t",
            quote=F,
            file="../figures/RefSeq_IDs.xls")
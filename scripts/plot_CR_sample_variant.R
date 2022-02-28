library(tidyverse)

variant_cr_before=read_tsv("data/HAIL_GWAS_vcf/var_CR_initial.tsv.bgz",skip=1,
                           col_names=c("locus","allele","CR"))

ggplot(variant_cr_before)+
  geom_histogram(aes(x=CR), bins=100)

sample_cr_before=read_tsv("data/HAIL_GWAS_vcf/sample_CR_initial.tsv.bgz",skip=1,
                          col_names=c("sample","CR"))

ggplot(sample_cr_before)+
  geom_histogram(aes(x=CR), bins=100)

var_cr_afterSample=read_tsv("data/HAIL_GWAS_vcf/var_CR_afterSampleCR.tsv.bgz",skip=1,
                            col_names=c("locus","allele","CR"))
ggplot(var_cr_afterSample)+
  geom_histogram(aes(x=CR), bins=100)


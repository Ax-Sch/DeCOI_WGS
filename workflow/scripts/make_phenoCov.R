library(tidyverse)

#cov<-read_tsv("/media/axel/Dateien/Arbeit_Gen/COVID_WGS/WGS_GWAS_pipeline/results/regenie_pheno/cov_A2_EUR")
#pheno<-read_tsv("/media/axel/Dateien/Arbeit_Gen/COVID_WGS/WGS_GWAS_pipeline/results/regenie_pheno/pheno_EUR_A2")

cov<-read_tsv(snakemake@input[[2]])%>%
  select(-FID)
pheno<-read_tsv(snakemake@input[[1]])%>%
  select(-FID)

pheno_cov<-pheno %>%
 left_join(cov, by = "IID")

write_tsv(x=pheno_cov,file=snakemake@output[[1]])

males<-pheno %>% 
  filter(sex_for_fam==1) %>%
  select(IID)

write_tsv(x=males, file=snakemake@output[[2]], col_names=FALSE)



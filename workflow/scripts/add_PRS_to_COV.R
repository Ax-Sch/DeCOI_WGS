library(tidyverse)

cov_path<-"results/var_sets/EURs_unrel/B2/cov.cov"
prs_path<-"results/PRS/scores.profile"
cov_out_path<-"results/var_sets/EURs_unrel/B2/cov.cov"

cov_path<-snakemake@input[[1]]
prs_path<-snakemake@input[[2]]
cov_out_path<-snakemake@output[[1]]

COV<-read_tsv(file=cov_path, col_names=TRUE)
PRS<-read.table(file = prs_path, header = TRUE)

print("Percentage of individuals in the COV file that also got a PRS calculated:")
mean(COV$IID %in% PRS$IID)

PRS_red<-PRS %>%
  select(IID, SCORESUM)

COV<-COV %>% 
  left_join(PRS_red, by="IID")

write_tsv(x=COV, file=cov_out_path, col_names = TRUE)
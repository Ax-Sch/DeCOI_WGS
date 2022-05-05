library(tidyverse)
library(readxl)

# Params:
# 1. input pheno-excel
# 2. input PC file
# 3. output file

args=c("config/post_df3_HGI_sample_QC_summary.xlsx","data/PCAcovar/A2_EUR_PCA.eigenvec")

# first: excel_file, second: fam-file, third: Phenotype column in excel

args = commandArgs(trailingOnly=TRUE)
print(args)

# File formats:
#Line 1 : Header with FID, IID and C covariate names.
#Followed by lines of C+2 values. Space/tab separated.

samplesheet<-read_tsv(args[1])%>%
  mutate(sex_for_fam=ifelse(is_female==TRUE, 2, ifelse(is_female==FALSE,1,0)))




col_names_eigenv<-c("FID","IID", paste0("PC",1:10))

# sex coding: 1=male, 2=female
# pheno coding: 1=control, 2=case
eigenv<-read_delim(args[2], col_names=col_names_eigenv, delim=" ")

cov_info<-eigenv %>% left_join(samplesheet, by=c("IID"="s"))%>%
  mutate(Age=as.integer(Age)) %>%
  mutate(age_sex=Age*sex_for_fam)%>%
  mutate(Age2=Age*Age)%>%
  select(all_of(col_names_eigenv), Age, Age2, age_sex, sex_for_fam)

write_tsv(x = cov_info, file = args[3], col_names = TRUE) # pheno file for regenie
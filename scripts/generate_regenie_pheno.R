library(tidyverse)
library(readxl)
# Params:
# 1. input pheno-excel
# 2. input fam file
# 3. pheno column name in excel
# 4. Output: fam file (Plink)
# 5. Output: pheno file (Regenie)

#testing
#setwd("/media/axel/Dateien/Arbeit_Gen/COVID_WGS/WGS_GWAS_pipeline")
args=c("config/post_df3_HGI_sample_QC_summary.xlsx","data/anc_vcf/EUR_vcf_not_rel.vcf.bgz.fam", "A2", 
       "data/regenie_pheno/fam_file.fam", "data/regenie_pheno/phenox")
# first: excel_file, second: fam-file, third: Phenotype column in excel

args = commandArgs(trailingOnly=TRUE)
print(args)

# File formats:
#FAM file: 
# A text file with no header line, and one line per sample with the following six fields:
#Family ID ('FID')
#Within-family ID ('IID'; cannot be '0')
#Within-family ID of father ('0' if father isn't in dataset)
#    Within-family ID of mother ('0' if mother isn't in dataset)
#Sex code ('1' = male, '2' = female, '0' = unknown)
#Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)

# Regnie pheno file
# FID IID Y1 Y2
#  0=control, 1=case, missing values must be coded as NA

samplesheet<-read_excel(args[1], sheet = "qc_table")%>% 
        rename(sex_to_use=`Sex (f=fem, m=m)`)%>%
  mutate(sex_for_fam=ifelse(sex_to_use=="f", 2, ifelse(sex_to_use=="m",1,0)))

pheno_col_num<-which(colnames(samplesheet) == args[3])
pheno_PLINK<-as.integer(unlist(samplesheet[,pheno_col_num]))+1
samplesheet$pheno_P=replace_na(pheno_PLINK, -9)

pheno_REGENIE<-as.integer(unlist(samplesheet[,pheno_col_num]))
samplesheet$pheno_R=pheno_REGENIE


col_names_fam<-c("FID","IID","father","mother","sex","pheno")

# sex coding: 1=male, 2=female
# pheno coding: 1=control, 2=case
fam<-read_tsv(args[2], 
col_names=col_names_fam,
col_types=c(
  FID = col_character(),
  IID = col_character(),
  father = col_character(),
  mother = col_character(),
  sex = col_integer(),
  pheno = col_integer()
)
)

all_info_joined<-fam %>% left_join(samplesheet, by=c("IID"="s"))%>%
  mutate(sex=sex_for_fam)

fam_new<-all_info_joined %>% 
  select(all_of(col_names_fam[1:4]), sex_for_fam, pheno_P)
write_tsv(x = fam_new, file = args[4], col_names = FALSE) # fam_file for plink

pheno_file_regenie <- all_info_joined %>%
  select(FID, IID, pheno_R)
colnames(pheno_file_regenie)<-c("FID","IID",args[3])
write_tsv(x = pheno_file_regenie, file = args[5], col_names = TRUE) # pheno file for regenie
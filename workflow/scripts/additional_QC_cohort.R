library(tidyverse)
#1d
# .lmiss files: variant based missingnness
missingness_case<-read_table("Geno05_CR_sex_snp_qc_snpqcCAS.lmiss")
missingness_ctrl<-read_table("Geno05_CR_sex_snp_qc_snpqcCON.lmiss")
missingness_comb<-missingness_case %>% 
  left_join(missingness_ctrl, by=c("CHR"="CHR","SNP"="SNP"))%>%
  mutate(miss_diff=abs(F_MISS.x-F_MISS.y))

length((missingness_comb %>% filter(miss_diff>=0.02))$SNP)
head(missingness_comb %>% filter(miss_diff>=0.02))

snps_to_filter_missingness<-missingness_comb %>% filter(miss_diff>=0.02) %>% select(SNP)

# hardy
hardy<-read_table("Geno05_CR_sex_snp_qc_snpqcALL.hwe")
head(hardy)

print("exclude by hardy weinberg")
hardy_weinberg_filtered<-hardy %>% filter(CHR!=23 & ((TEST=="AFF" & P<1e-10) | (TEST=="UNAFF" & P<1e-6))  ) %>% distinct(SNP)
head(hardy_weinberg_filtered)
hardy_females<-read_table("Geno05_CR_sex_snp_qc_snpqcALL_female.hwe")
## just ChrX - very few females in cases >> taking the ALL TEST
hardy_weinberg_filtered_xchr<-hardy_females %>% filter((TEST=="ALL" & P<1e-6)) %>% distinct(SNP)

filter_missingness_hardy<-rbind(snps_to_filter_missingness, hardy_weinberg_filtered, hardy_weinberg_filtered_xchr) %>% distinct()

write_tsv(filter_missingness_hardy, "missingness_hardy_weinberg_filter.txt", col_names = FALSE)
## used to do 1e final QC

library(tidyverse)
library(plotly)

wgs_plink<-read_tsv("../results/regenie/GWAS_EUR_B2.PHENO1.glm.logistic")%>%
  rename(CHROM=`#CHROM`)%>%
  mutate(CHROM=str_replace(CHROM, fixed("X"), "23")) %>%
  mutate(CHROM=as.numeric(CHROM))


wgs_regenie<-read_tsv("../results/regenie/GWAS_EUR_B2.regenie.gz")

gwas_covid19HGI<-read_tsv("/mnt/int1/COVID_HGI_GWAS/COVID19_HGI_B2_ALL_leave_23andme_20210607.txt.gz")
gwas_covid19HGI_log<-gwas_covid19HGI %>% 
  mutate(LOG10P=-log10(all_inv_var_meta_p)) %>%
  rename(CHROM=`#CHR`)

wgs_genomicc<-read_tsv("/mnt/int1/COVID_HGI_GWAS/GenOMICC_EUR.txt.gz")

wgs_genomicc_log<-wgs_genomicc%>%
  mutate(LOG10P=-log10(Pval)) %>%
  rename(CHROM=`CHR`)

wgs_plink_filtered <- wgs_plink %>% 
  filter(LOG10_P>2)

wgs_regenie_filtered <- wgs_regenie %>% 
  rename(CHROM=Chr)%>%
  mutate(LOG10P=-log10(Pval),
         GENPOS=Pos)%>%
  filter(LOG10P>2)

gwas_covid19HGI_log_filtered <- gwas_covid19HGI_log %>%
  filter(LOG10P>3)

wgs_genomicc_log_filtered <- wgs_genomicc_log %>%
  filter(LOG10P>3)


make_man_plot_chrom<-function(chrom){
lala<-ggplot()+
  geom_point(data=gwas_covid19HGI_log_filtered %>%filter(CHROM==chrom), aes(x=POS, y=LOG10P), color="green")+
  geom_point(data=wgs_genomicc_log_filtered %>%filter(CHROM==chrom), aes(x=POS, y=LOG10P), color="brown")+
  geom_point(data=wgs_regenie_filtered %>%filter(CHROM==chrom), aes(x=GENPOS, y=LOG10P), color="blue")+
  geom_point(data=wgs_plink_filtered %>%filter(CHROM==chrom), aes(x=POS, y=LOG10_P))

ggplotly(lala)
}

make_man_plot_chrom(3)

abc<-gwas_covid19HGI_log %>%
  left_join(wgs_genomicc_log, by=c("CHROM"="CHROM", "POS"="POS", "REF"="Non.Effect.Allele"))

cor.test(abc$LOG10P.x, abc$LOG10P.y)


abc<-wgs_plink %>%
  left_join(wgs_regenie, by=c("CHROM"="CHROM", "POS"="GENPOS", "REF"="ALLELE0"))

cor.test(abc$LOG10_P, abc$LOG10P)


abc<-wgs_plink %>%
  left_join(wgs_genomicc_log, by=c("CHROM"="CHROM", "POS"="POS", "REF"="Non.Effect.Allele"))

cor.test(abc$LOG10_P, abc$LOG10P)

abc<-wgs_regenie %>%
  left_join(wgs_genomicc_log, by=c("CHROM"="CHROM", "GENPOS"="POS", "ALLELE0"="Non.Effect.Allele"))

cor.test(abc$LOG10P.x, abc$LOG10P.y)




abc<-wgs_plink %>% filter(LOG10_P>1)%>%
  left_join(wgs_genomicc_log, by=c("CHROM"="CHROM", "POS"="POS", "REF"="Non.Effect.Allele"))

cor.test(abc$LOG10_P, abc$LOG10P)


abc<-wgs_regenie %>% filter(LOG10P>1)%>%
  left_join(wgs_genomicc_log, by=c("CHROM"="CHROM", "GENPOS"="POS", "ALLELE0"="Non.Effect.Allele"))

cor.test(abc$LOG10P.x, abc$LOG10P.y)


abc<-wgs_plink %>% filter(LOG10_P>2)%>%
  left_join(wgs_genomicc_log, by=c("CHROM"="CHROM", "POS"="POS", "REF"="Non.Effect.Allele"))

cor.test(abc$LOG10_P, abc$LOG10P)


abc<-wgs_regenie %>% filter(LOG10P>2)%>%
  left_join(wgs_genomicc_log, by=c("CHROM"="CHROM", "GENPOS"="POS", "ALLELE0"="Non.Effect.Allele"))

cor.test(abc$LOG10P.x, abc$LOG10P.y)







abc<-wgs_plink %>% filter(LOG10_P>2)%>%
  left_join(gwas_covid19HGI_log, by=c("CHROM"="CHROM", "POS"="POS", "REF"="REF"))

cor.test(abc$LOG10_P, abc$LOG10P)


abc<-wgs_regenie %>% filter(LOG10P>2)%>%
  left_join(gwas_covid19HGI_log, by=c("CHROM"="CHROM", "GENPOS"="POS", "ALLELE0"="REF"))

cor.test(abc$LOG10P.x, abc$LOG10P.y)

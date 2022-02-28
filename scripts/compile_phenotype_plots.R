library(tidyverse)
library(readxl)

args=c("config/post_df3_HGI_sample_QC_summary.xlsx", "data/hail_gather_data/for_sample_QC.tsv", "data/plotPhenoInfo/")

args = commandArgs(trailingOnly=TRUE)
print(args)
PhenoFileLocation=args[1]
qc_file<-args[2]
outFolder=args[3]

samplesheet<-read_excel(PhenoFileLocation, sheet = "qc_table")

qc_table<-read_tsv(qc_file)

dir.create(outFolder,recursive = TRUE, showWarnings = FALSE)
setwd(outFolder)

samplesheet<-samplesheet%>% 
  rename(sex_to_use=`Sex (f=fem, m=m)`)%>%
  mutate(sex_for_fam=ifelse(sex_to_use=="f", 2, ifelse(sex_to_use=="m",1,0)))

sample_sex<-samplesheet %>% select(sex_to_use, s)

qc_table<-qc_table %>%
  left_join(sample_sex, by=(c("s"="s")))

qc_table<-qc_table %>% 
  mutate(Age=as.integer(Age),
         A2=as.integer(A2),
         B2=as.integer(B2),
         C2=as.integer(C2),
         is_case=as.integer(is_case)
         )



# Check Age vs. Sex

SexPlot<-ggplot(qc_table)+
  geom_histogram(aes(x=Age, fill=sex_to_use), bins=20)+
  theme_bw()

print(SexPlot)
ggsave(filename="SexPlot.pdf",
       plot=SexPlot,
       width=7,
       height=3)

SexInfo<-qc_table %>% 
  summarise(mean(sex_to_use=="m"), num=sum(sex_to_use!="noSex"))
write_tsv(x=SexInfo, file="SexInfo.tsv")

# Check Age vs. Severity
qc_table <- qc_table %>% 
  mutate(is_case_F=as.factor(is_case))

AgePlot<-ggplot(qc_table)+
  geom_histogram(aes(x=Age, fill=is_case_F), bins=20)+
  geom_vline(data=qc_table %>% filter(is_case==0), aes(xintercept = mean(Age)), col="red")+
  geom_vline(data=qc_table %>% filter(is_case==1), aes(xintercept = mean(Age)), col="green")+
  geom_vline(data=qc_table %>% filter(is_case==2), aes(xintercept = mean(Age)), col="blue")+
  theme_bw()

print(AgePlot)
ggsave(filename="AgePlot.pdf",
       plot=AgePlot,
       width=7,
       height=3)

AgeInfo<-qc_table %>% 
  group_by(is_case)%>%
  summarise(MEAN=mean(Age), SD=sd(Age), N=length(Age))%>%
  mutate(SEM=SD/sqrt(N))
  

write_tsv(x=AgeInfo, file="AgeInfo.tsv")
  
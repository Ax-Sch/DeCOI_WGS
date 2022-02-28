library(tidyverse)
library(readxl)
sampleQC<-read_tsv("dataset_split_cols.tsv")
n_cols<-ncol(sampleQC)
samplesheet<-read_excel("../../config/post_df3_HGI_sample_QC_summary.xlsx", sheet = "qc_table")[,1:13] %>%
  filter(is.na(`missing/not to use df3`))%>%
  filter(!is.na(A2))%>%
  distinct(s, .keep_all=TRUE) 

samplesheet_w_sampleQC <- samplesheet %>% 
  left_join(sampleQC, by=c("s"="value"))

# correct sex info
samplesheet_w_sampleQC <- samplesheet_w_sampleQC %>%
  mutate(impute_sex2_is_female=as.logical(impute_sex2_is_female))%>%
  rename(sex_samplesheet=`Sex (f=fem, m=m)`)%>%
  mutate(sex_samplesheet_bool=ifelse(sex_samplesheet=="f", TRUE, 
                          ifelse(sex_samplesheet=="m", FALSE, NA)))%>%
  mutate(sex_problem=(is.na(sex_samplesheet_bool) | 
                      is.na(impute_sex2_is_female) | 
                      impute_sex2_is_female!=sex_samplesheet_bool))


# QC criteria
samplesheet_w_sampleQC<-samplesheet_w_sampleQC %>% 
  mutate(filtered_cause=ifelse(sample_qc6_call_rate < 0.90, "CR", ""))%>%
  mutate(filtered_cause=ifelse(sample_qc_dp_stats_mean < 20, paste(filtered_cause, "DP", sep="_"), filtered_cause))%>%
  mutate(filtered_cause=ifelse(sample_qc_dp_stats_stdev > 20, paste(filtered_cause, "DPSD", sep="_"), filtered_cause))%>%
  mutate(filtered_cause=ifelse(sample_qc_gq_stats_mean < 45, paste(filtered_cause, "GQ", sep="_"), filtered_cause))%>%
  mutate(filtered_cause=ifelse(sample_qc_gq_stats_stdev > 35, paste(filtered_cause, "GQSD", sep="_"), filtered_cause))%>%
  mutate(filtered_cause=ifelse(sex_problem, paste(filtered_cause, "SexProb", sep="_"), filtered_cause))%>%
  mutate(filtered=filtered_cause!="")
  
knitr::kable(table(samplesheet_w_sampleQC$filtered_cause))

  
# add column describing if the QC is passed
samplesheet_w_sampleQC<-samplesheet_w_sampleQC %>% 
  mutate(qc=ifelse(!(filtered | sex_problem), 
            ifelse(is.na(impute_sex2_is_female) | is.na(sex_samplesheet_bool), "1sexMIS", "0ok"), 
            ifelse(filtered & sex_problem, "4filtered_plus_sexmm",
            ifelse(filtered, "2filtered", "3sexmm"))))

knitr::kable(table(samplesheet_w_sampleQC$qc))

samplesheet_w_sampleQC$batch<-strtrim(samplesheet_w_sampleQC$s,2)
d_types=sapply(samplesheet_w_sampleQC, class)

dir.create("filtered", showWarnings = FALSE)

for (col in 1:n_cols){
  if (d_types[col] %in% c("integer","numeric")){
    vals_as_vect<-unlist(samplesheet_w_sampleQC[,col])
    median_val=median(vals_as_vect)
    IQR_val=quantile(vals_as_vect, 3/4, na.rm=T) - quantile(vals_as_vect, 1/4, na.rm=T)
    upper=quantile(vals_as_vect, 3/4, na.rm=T)+ 1.5 * IQR_val
    lower=quantile(vals_as_vect, 1/4, na.rm=T)- 1.5 * IQR_val
    
    p<-ggplot(samplesheet_w_sampleQC)+
      geom_histogram(aes_string(x=colnames(samplesheet_w_sampleQC)[col], fill="filtered_cause"), bins=100)+
      geom_vline(xintercept=median_val)+
      geom_vline(xintercept=upper)+
      geom_vline(xintercept=lower)+
      annotate("text", label = round(upper,3), x = upper, y = 5) +
      annotate("text", label = round(lower,3), x = lower, y = 5) 
    ggsave(filename=paste0("filtered/", colnames(samplesheet_w_sampleQC)[col],"_plotted_filtered.png"), plot=p, height=7, width=7)
  }
}

dir.create("filtered_sex", showWarnings = FALSE)

for (col in 1:n_cols){
  if (d_types[col] %in% c("integer","numeric")){
    p<-ggplot(samplesheet_w_sampleQC)+
      geom_histogram(aes_string(x=colnames(samplesheet_w_sampleQC)[col], fill="qc"), bins=100)
    ggsave(filename=paste0("filtered_sex/", colnames(samplesheet_w_sampleQC)[col],"_plotted_all_qc.png"), plot=p)
  }
}

dir.create("center", showWarnings = FALSE)

for (col in 1:n_cols){
  if (d_types[col] %in% c("integer","numeric")){
    vals_as_vect<-unlist(samplesheet_w_sampleQC[,col])
    median_val=median(vals_as_vect)
    IQR_val=quantile(vals_as_vect, 3/4, na.rm=T) - quantile(vals_as_vect, 1/4, na.rm=T)
    upper=quantile(vals_as_vect, 3/4, na.rm=T)+ 1.5 * IQR_val
    lower=quantile(vals_as_vect, 1/4, na.rm=T)- 1.5 * IQR_val
    
    p<-ggplot(samplesheet_w_sampleQC)+
      geom_histogram(aes_string(x=colnames(samplesheet_w_sampleQC)[col], fill="batch"), bins=100)+
      geom_vline(xintercept=median_val)+
      geom_vline(xintercept=upper)+
      geom_vline(xintercept=lower)+
      annotate("text", label = round(upper,3), x = upper, y = 5) +
      annotate("text", label = round(lower,3), x = lower, y = 5) 
    ggsave(filename=paste0("center/", colnames(samplesheet_w_sampleQC)[col],"_plotted_colored.png"), plot=p)
  }
}

hets<-qplot(samplesheet_w_sampleQC$sample_qc6_n_het/samplesheet_w_sampleQC$sample_qc6_n_called)
ggsave(filename="heterozygosity.png")

export_file<-samplesheet_w_sampleQC %>% filter(filtered==FALSE, sex_problem==FALSE, is_case!=-9) %>% select(s, filtered, sex_problem, is_case, A2, B2, C2, Age)
write_tsv(x=export_file, file="for_sample_QC.tsv")

export_file<-samplesheet_w_sampleQC %>% filter(sex_problem==FALSE, is_case!=-9) %>% select(s, filtered, filtered_cause,sample_qc_dp_stats_mean,  sample_qc6_call_rate, sex_problem, is_case, A2, B2, C2, Age)
write_tsv(x=export_file, file="for_sample_QC_w_cause.tsv")

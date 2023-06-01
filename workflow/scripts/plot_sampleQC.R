library(tidyverse)
library(readxl)
library(plotly)

sampleQC<-read_tsv("sampleqc.tsv")
n_cols<-ncol(sampleQC)
samplesheet<-read_excel("post_df3_HGI_sample_QC_summary.xlsx", sheet = "qc_table")[,1:13] %>%
  filter(is.na(`missing/not to use df3`))%>%
  filter(!is.na(A2))%>%
  distinct(s, .keep_all=TRUE) 

  samplesheet_w_sampleQC <- samplesheet %>% 
  left_join(sampleQC, by="s")





# correct sex info
samplesheet_w_sampleQC <- samplesheet_w_sampleQC %>%
  rename(sex_samplesheet=`Sex (f=fem, m=m)`)%>%
  mutate(sex_samplesheet_bool=ifelse(sex_samplesheet=="f", TRUE, 
                          ifelse(sex_samplesheet=="m", FALSE, NA)))


# QC criteria
samplesheet_w_sampleQC<-samplesheet_w_sampleQC %>% 
  mutate(filtered_cause=ifelse(sample_qc.call_rate < 0.90, "CR", ""))%>%
  mutate(filtered_cause=ifelse(sample_qc.dp_stats.mean < 20, paste(filtered_cause, "DP", sep="_"), filtered_cause))%>%
  mutate(filtered_cause=ifelse(sample_qc.dp_stats.stdev > 20, paste(filtered_cause, "DPSD", sep="_"), filtered_cause))%>%
  mutate(filtered_cause=ifelse(sample_qc.gq_stats.mean < 45, paste(filtered_cause, "GQ", sep="_"), filtered_cause))%>%
  mutate(filtered_cause=ifelse(sample_qc.gq_stats.stdev > 35, paste(filtered_cause, "GQSD", sep="_"), filtered_cause))%>%
  mutate(filtered=filtered_cause!="")
  
knitr::kable(table(samplesheet_w_sampleQC$filtered_cause))

  
# add column describing if the QC is passed


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


#dir.create("center", showWarnings = FALSE)

#for (col in 1:n_cols){
#  if (d_types[col] %in% c("integer","numeric")){
#    vals_as_vect<-unlist(samplesheet_w_sampleQC[,col])
#    median_val=median(vals_as_vect)
#    IQR_val=quantile(vals_as_vect, 3/4, na.rm=T) - quantile(vals_as_vect, 1/4, na.rm=T)
#    upper=quantile(vals_as_vect, 3/4, na.rm=T)+ 1.5 * IQR_val
#    lower=quantile(vals_as_vect, 1/4, na.rm=T)- 1.5 * IQR_val
    
#    p<-ggplot(samplesheet_w_sampleQC)+
#      geom_histogram(aes_string(x=colnames(samplesheet_w_sampleQC)[col], fill="batch"), bins=100)+
#      geom_vline(xintercept=median_val)+
#      geom_vline(xintercept=upper)+
#      geom_vline(xintercept=lower)+
#      annotate("text", label = round(upper,3), x = upper, y = 5) +
#      annotate("text", label = round(lower,3), x = lower, y = 5) 
#    ggsave(filename=paste0("center/", colnames(samplesheet_w_sampleQC)[col],"_plotted_colored.png"), plot=p)
#  }
#}

hets<-qplot(samplesheet_w_sampleQC$sample_qc.n_het/samplesheet_w_sampleQC$sample_qc.n_called)
ggsave(filename="heterozygosity.png")

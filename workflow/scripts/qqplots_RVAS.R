if (!require("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel", repos='http://cran.us.r-project.org')
}



library(data.table)
library(tidyverse)
library(ggrepel)

in_file="results/RVAS/rRVAS_max_EUR_POP_A2_PRS/0.001_A2.regenie.gz"
in_file<-"../GWAS/EUR_A2.regenie.gz"
in_file<-"../plink_RVAS/EUR_POP_A2_A2.model"
in_file=snakemake@input[[1]]

ENSG_to_name_path="resources/ENSG_to_name.txt"
ENSG_to_name_path=snakemake@input[[2]]

BURDEN_TEST=TRUE
BURDEN_TEST=snakemake@params[[1]]

con <- file(in_file,"r")
first_line <- readLines(con,n=1)
line11<-readLines(con,n=10)[10]
close(con)

print(paste("First line:", first_line))
print(paste("Line 11:", line11))

# multi-space delimited?
if (grepl(pattern = "  ",x = line11, )){
  if (grepl(pattern="CHR", x=first_line)){
    HEADER=TRUE
  }else{
    HEADER=FALSE
  }
  print("Using read.table as double space was detected in line 11")
  assoc_data<-read.table(in_file, header=HEADER)
}else if(grepl(pattern = "##",x = first_line)){
  #comments at the beginning?
  print("Using fread and skipping lines containing #")
  assoc_data<-fread(file=in_file, skip ="#", blank.lines.skip=TRUE)
}else{
  print("Using fread")
  assoc_data<-fread(file=in_file, blank.lines.skip=TRUE)
}

head(assoc_data)

assoc_data<-as_tibble(assoc_data)

# add p-value column
# check for LOG10
clns<-colnames(assoc_data)

if ("LOG10P" %in% clns){
  assoc_data <- assoc_data %>%
    mutate(P=10^(-LOG10P))
}

# rename
PVAL_NAMES<-c("Pval","PVAL", "P")
for (p_name in PVAL_NAMES){
  if (p_name %in% clns){
    assoc_data$P <- unlist(assoc_data[,p_name])
  }
}

SNP_NAMES<-c("ID")
for (snp_name in SNP_NAMES){
  if (snp_name %in% clns){
    assoc_data$SNP <- unlist(assoc_data[,snp_name])
  }
}

# add gene Name/Symbol
if (BURDEN_TEST){   #grepl(pattern="ENSG", assoc_data$SNP[1])){
  #BURDEN_TEST=TRUE
  ENSG_to_name<-read_tsv(ENSG_to_name_path, col_names = c("ENSG", "HGNC", "ENS_NAME"), skip=1)%>%
    filter(!is.na(ENS_NAME))
  SNP_info<-str_split(assoc_data$SNP, pattern = fixed("."), simplify = TRUE, n=3)
  colnames(SNP_info)<- c("ENSG", "MASK", "AF_GROUP")
  assoc_data<-cbind(assoc_data, SNP_info)
  
  assoc_data<-assoc_data%>%
    left_join(ENSG_to_name, by="ENSG")%>%
    mutate(gene_ensg=ifelse(is.na(HGNC), ENSG, HGNC))%>%
    mutate(SNP=paste(gene_ensg, MASK, AF_GROUP, sep="."))%>%
    filter(AF_GROUP!="singleton")
}

assoc_data <- assoc_data%>%
  filter(!is.na(P))%>%
  arrange(P)


makeQQ<-function(ASSOC_DATA, suff, in_file){
  ci <- 0.95
  nSNPs <- length(ASSOC_DATA$P)
  
  plotdata <- data.frame(
    observed = -log10((ASSOC_DATA$P)),
    SNP=ASSOC_DATA$SNP,
    expected = -log10(ppoints(nSNPs)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs))))
  )
  

  plotdata_sub <- plotdata %>%
    filter(expected <= 2)
  
  # uncomment, if subsampling should be activated
  #n_r<-nrow(plotdata_sub)
  #sze<-ifelse(2000<n_r, 2000, n_r)
  
  #print(sze)
  
  #plotdata_sub <- plotdata_sub %>%
  #  sample_n(size=sze)
  
  plotdata_sup <- plotdata %>%
    filter(expected > 2)
  
  plotdata_small <- rbind(plotdata_sub, plotdata_sup)
  
  plotdata_small<-plotdata_small %>%
    arrange(-observed)%>%
    mutate(ind=1:n())%>%
    mutate(top10=ind<=5)
  
  qqplot <- ggplot(plotdata_small, aes(x = expected, y = observed)) +
    geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
    geom_point(color = "black", size = 1.1, direction = "vh") +
    geom_segment(data = . %>% filter(expected == max(expected)), 
                 aes(x = 0, xend = expected, y = 0, yend = expected),
                 size = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
    labs(x = expression(paste("Expected -log"[10],"(", plain(P),")")),
         y = expression(paste("Observed -log"[10],"(", plain(P),")"))) +
    theme_minimal()
    ggsave(filename = paste0(in_file, suff, ".pdf"), plot = qqplot, height=5, width=5)
  
    qqplot_w_label<- qqplot+ 
      geom_text_repel(data=plotdata_small %>% filter(top10==TRUE),aes(label=SNP))
  
    ggsave(filename = paste0(in_file, suff, "Label.pdf"), plot = qqplot_w_label, height=5, width=5)
  
  return(qqplot)
}

# plot with all masks/af cut offs
qqplot<-makeQQ(assoc_data, "QQ", in_file)


#print(qqplot)

if (BURDEN_TEST){
# plot for each mask/af cut off combination
for (mask in unique(assoc_data$MASK)){
  for (af_group in unique(assoc_data$AF_GROUP)){
    temp_assoc_data<-assoc_data%>%
      filter(MASK==mask, AF_GROUP==af_group)
    suffix<-paste(mask, af_group,"QQ", sep="_")
    qqplot<-makeQQ(temp_assoc_data,suffix , in_file)
  }
}
}

# export
write_tsv(x = assoc_data[1:min(5000, nrow(assoc_data)),], file = paste0(in_file, "_head.tsv"))





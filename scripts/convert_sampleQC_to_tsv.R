library(tidyverse)

base_dataset<-read_tsv("sampleqc.tsv")
final_dataset<-as_tibble(base_dataset$s)

for (col in 2:ncol(base_dataset)){
dataset_first_split<- lapply(base_dataset[,col], function(x) {str_split(x, "\\{|\\}", simplify=TRUE)} ) %>%
  as.data.frame(stringsAsFactors=FALSE) %>% type.convert(as.is=TRUE)
prefix<-gsub("[^A-Za-z0-9_ ]","",colnames(base_dataset)[col])

col_names_temp<-colnames(dataset_first_split)
data_cols=c()
for (i in 1:ncol(dataset_first_split)){
  n_char_first=nchar(dataset_first_split[1,i], keepNA=FALSE)
  if (n_char_first < 30  & n_char_first > 2 & i!=ncol(dataset_first_split)){
    col_names_temp[i+1]<-paste(prefix, unique(dataset_first_split[,i])[1],sep="_")
    data_cols<-c(data_cols, FALSE)
  }else{
    data_cols<-c(data_cols, TRUE)
  }
}
colnames(dataset_first_split)<-col_names_temp
dataset_first_split<-dataset_first_split[, data_cols]
d_types=sapply(dataset_first_split, class)
dataset_first_split<-dataset_first_split[d_types!="logical"]

for (i in 1:ncol(dataset_first_split)){
  col_second_split<-as_tibble(str_split(dataset_first_split[,i], ",|:", simplify=TRUE)) %>% 
    type.convert(as.is=TRUE)
  name_of_col<-gsub("[^A-Za-z0-9_ ]","",colnames(dataset_first_split)[i])
  d_types_loop=sapply(col_second_split, class)
  k=1 # remove logical cols at the beginning
  while (d_types_loop[k]=="logical"){
    k=k+1
  }
  col_second_split=col_second_split[,k:ncol(col_second_split)]
  d_types_loop=sapply(col_second_split, class)
  col_names<-colnames(col_second_split)
 for (h in 1:ncol(col_second_split)){
    if (h %% 2 == 1){
      colname_append<- gsub("[^A-Za-z0-9_]","",col_second_split[1,h])
      col_names[h+1]<-paste(name_of_col,colname_append, sep="_")
      print(paste(h, colname_append))
    }
  }
  colnames(col_second_split)<-col_names
  col_second_split<-col_second_split[,(1:ncol(col_second_split)) %% 2 == 0]
  final_dataset<-cbind(final_dataset,col_second_split)
}
}

write_tsv(x=final_dataset, file="dataset_split_cols.tsv")


library(randomForest)
library(tidyverse)
library(plotly)

# read in the eigenvectors, produced in PLINK


args=c('MergeFullForPCA.eigenvec', '20130606_g1k.ped', '../../resources/covid_individuals.txt')
# first: excel_file, second: fam-file, third: Phenotype column in excel

#args = commandArgs(trailingOnly=TRUE)
#print(args)

eigenvec_path<-args[1]
ped_path<-args[2]
phenotype_path<-args[3]


eigenvec <- read.table(eigenvec_path, header = FALSE, skip=0, sep = ' ')
eigenvec <- eigenvec[,2:ncol(eigenvec)]
head(eigenvec)
colnames(eigenvec) <- c("Individual.ID",paste('PC', c(1:20), sep = ''))

# read in the PED data
PED <- read.table(ped_path, header = TRUE, skip = 0, sep = '\t')

pheno<-read_tsv(phenotype_path)

eigenvec<-eigenvec %>% left_join(pheno, by=c("Individual.ID"="s"))

#build data frame for random forest classifier
dataRF <- merge(eigenvec, PED[, c("Individual.ID", "Population")], all.x=TRUE)

#build plot
dataRF$Population <- factor(dataRF$Population, levels=c(
  "ACB","ASW","ESN","GWD","LWK","MSL","YRI",
  "CLM","MXL","PEL","PUR",
  "CDX","CHB","CHS","JPT","KHV",
  "CEU","FIN","GBR","IBS","TSI",
  "BEB","GIH","ITU","PJL","STU"))

dataRF$Continental <- rep(NA_character_, nrow(dataRF))
dataRF$Continental[which(dataRF$Population %in% c("ACB","ASW","ESN","GWD","LWK","MSL","YRI"))]<-"AFR"
dataRF$Continental[which(dataRF$Population %in% c("CLM","MXL","PEL","PUR"))]<-"AMR"
dataRF$Continental[which(dataRF$Population %in% c("CDX","CHB","CHS","JPT","KHV"))]<-"EAS"
dataRF$Continental[which(dataRF$Population %in% c("CEU","FIN","GBR","IBS","TSI"))]<-"EUR"
dataRF$Continental[which(dataRF$Population %in% c("BEB","GIH","ITU","PJL","STU"))]<-"SAS"
dataRF$Continental<-as.factor(dataRF$Continental)

rf_classifier = randomForest(Continental ~ ., 
                             data=dataRF[which(!is.na(dataRF$Continental)), c("PC1","PC2","PC3","PC4","PC5","PC6", "Continental")], 
                             ntree=3000, importance=TRUE)

#predict population in your cohort
dataPred<-dataRF[which(is.na(dataRF$Continental)),]
dataRF<-dataRF[which(!is.na(dataRF$Continental)),]
dataPred$Prediction<-rep(NA, nrow(dataPred))
dataPred$Prediction<-predict(rf_classifier,dataPred[,c("PC1","PC2","PC3","PC4","PC5","PC6")])

#outliers<-c()
# manuell outlier entfernt:
#dataPred


#outliers <- c("FO13863x01_02","FO14026x01_02","DE17BOSUKDD100080","DE39BOSUKDD100072","DE93BOSUKDD100070","DE06BOSUKDD100084","FO14016x01_02","DE38BOSUKED100009","DE76BOSUKDD100085","DE66BOSUKDD100071","DE88BOSUKDD100063","FO14344x02_02","DE37BOSUKDD100011","FO14432x01_02")

#dataPred <- dataPred %>% 
#  mutate(Prediction=as.factor(ifelse(Prediction=="EUR" & (PC1>-0.0067 | PC2<1.5e-2 | PC3<2e-3), "UNK", as.character(Prediction))))


PC1_2_plot<-ggplot()+
  geom_point(data=dataRF, aes(x=PC1,y=PC2,shape=Continental), alpha=0.5,size=3)+ 
  geom_point(data=dataPred, aes(x=PC1,y=PC2,color=Prediction, text=Individual.ID),shape=1, alpha=0.8,size=3 )+
  theme_bw()+
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))

ggsave(filename="PC1_2_plot.pdf",
       plot=PC1_2_plot,
       width=4,
       height=3.5)
ggplotly(PC1_2_plot)
  

ggplotly(ggplot()+
           geom_point(data=dataRF, aes(x=PC3,y=PC4,shape=Continental), alpha=0.5,size=3)+ 
           geom_point(data=dataPred, aes(x=PC3,y=PC4,color=Prediction, text=Individual.ID),shape=1, alpha=0.8,size=2 ) )

ggplotly(ggplot()+
           geom_point(data=dataRF, aes(x=PC5,y=PC6,shape=Continental), alpha=0.5,size=3)+ 
           geom_point(data=dataPred, aes(x=PC5,y=PC6,color=Prediction, text=Individual.ID),shape=1, alpha=0.8,size=2 ) )

ggplotly(ggplot()+
           geom_point(data=dataRF, aes(x=PC7,y=PC8,shape=Continental), alpha=0.5,size=3)+ 
           geom_point(data=dataPred, aes(x=PC7,y=PC8,color=Prediction, text=Individual.ID),shape=1, alpha=0.8,size=2 ) )

# lastly write IDs per population in files to be used for further analyses
write_tsv(dataPred %>% select(Individual.ID, Prediction),file="populations.txt")

### by seq center
dataPred_center<-dataPred %>%
  mutate(seq_center=ifelse(grepl("^FO|^DX",Individual.ID), "Tuebingen", NA))%>%
  mutate(seq_center=ifelse(grepl("^IT|^DE",Individual.ID), "Bonn", seq_center))%>%
  mutate(seq_center=ifelse(is.na(seq_center), "Koeln", seq_center))

#case control

ggplotly(ggplot()+
           geom_point(data=dataPred_center, aes(x=PC1,y=PC2,shape=Prediction, color=as.factor(is_case)), alpha=0.8,size=2)+
           theme_bw())

ggplotly(ggplot()+
           geom_point(data=dataPred_center, aes(x=PC3,y=PC4,shape=Prediction, color=as.factor(is_case)), alpha=0.8,size=2)+
           theme_bw())

ggplotly(ggplot()+
           geom_point(data=dataPred_center, aes(x=PC5,y=PC6,shape=Prediction, color=as.factor(is_case)), alpha=0.8,size=2)+
           theme_bw())

ggplotly(ggplot()+
           geom_point(data=dataPred_center, aes(x=PC7,y=PC8,shape=Prediction, color=as.factor(is_case)), alpha=0.8,size=2)+
           theme_bw())

ggplotly(ggplot()+
           geom_point(data=dataPred_center, aes(x=PC9,y=PC10,shape=Prediction, color=as.factor(is_case)), alpha=0.8,size=2)+
           theme_bw())






ggplotly(ggplot()+
           geom_point(data=dataPred_center, aes(x=PC1,y=PC2,shape=Prediction, color=seq_center), alpha=0.8,size=2)+
           theme_bw())

ggplotly(ggplot()+
           geom_point(data=dataPred_center, aes(x=PC3,y=PC4,shape=Prediction, color=seq_center), alpha=0.8,size=2)+
           theme_bw())

ggplotly(ggplot()+
           geom_point(data=dataPred_center, aes(x=PC5,y=PC6,shape=Prediction, color=seq_center), alpha=0.8,size=2)+
           theme_bw())

ggplotly(ggplot()+
           geom_point(data=dataPred_center, aes(x=PC7,y=PC8,shape=Prediction, color=seq_center), alpha=0.8,size=2)+
           theme_bw())

ggplotly(ggplot()+
           geom_point(data=dataPred_center, aes(x=PC9,y=PC10,shape=Prediction, color=seq_center), alpha=0.8,size=2)+
           theme_bw())




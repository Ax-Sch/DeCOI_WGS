library(tidyverse)

create_AF_plots<-function(vars_private) {

ggplot(vars_private %>% distinct(ID, .keep_all=T))+
  geom_histogram(aes(x=pmax(0,gnomAD_ge_AC, na.rm=T), fill=as.factor(AF_category)), bins=200)+
  xlim(c(-10,190))

ggplot(vars_private %>% distinct(ID, .keep_all=T))+
  geom_histogram(aes(x=pmax(0,gnomAD_ge_AC, na.rm=T), fill=as.factor(AF_category)), bins=30)+
  xlim(c(-2,27))

ggplot(vars_private %>% distinct(ID, .keep_all=T))+
  geom_histogram(aes(x=pmax(0,gnomAD_ge_nhomalt, na.rm=T), fill=as.factor(AF_category)), bins=30)+
  xlim(c(-1,29))

ggplot(vars_private, aes(x=gnomAD_ge_AC, y=(`C(HET)`+2*`C(HOM A1)`)))+
  geom_jitter(alpha=0.1)+
  xlim(c(-0.1,1000))+
  geom_smooth()

ggplot(vars_private, aes(x=gnomAD_AMR_AF, y=(`C(HET)`+2*`C(HOM A1)`)))+
  geom_point(alpha=0.1)+
  xlim(c(-0.001,0.01))+
  geom_smooth()

}

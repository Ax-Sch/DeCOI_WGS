library(tidyverse)

check_S_path<-snakemake@input[[1]]
plot_out<-snakemake@output[[1]]


check_S<-read.table(check_S_path, header = T)
F_plot<-ggplot(check_S, aes(x=F, fill=PEDSEX==1 | PEDSEX==0))+
  geom_histogram(bins=100)

ggsave(filename = plot_out, plot = F_plot)
ggsave(filename = paste0(plot_out, "zoomed.jpg"), plot = F_plot + ylim(c(0,100)) )

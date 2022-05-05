# QQ plot incl. CI intervals

library(data.table)
library(dplyr)
library(janitor)   # only available under R 3.5.1!
library(ggplot2)

setwd("C:/Users/Axel/Downloads/")

meta <- fread("GWAS_EUR_B2.regenie") %>% clean_names()
meta$p_value<-10^-meta$log10p
print(paste("preparing qq plot data, starts at ", date()))
ci <- 0.95
nSNPs <- nrow(meta)
lambda_calculated <- median(qchisq(meta$p_value, 1, lower.tail = F), na.rm = T) / qchisq(0.5, 1)
print(lambda_calculated) 

plotdata <- data.frame(
  observed = -log10(sort(meta$p_value)), 
  expected = -log10(ppoints(nSNPs)), 
  ci_lower = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))), 
  ci_upper = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs))))
)

print(paste("making qq plot, starts at ",date()))
qqplot <- ggplot(plotdata, aes(x = expected, y = observed)) + 
  geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), fill = "grey40", alpha = 0.5) + 
  geom_point(color = "black", size = 1.1) + 
  geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               size = 1.25, alpha = 0.5, color = "grey40", lineend = "round") +
  labs(x = expression(paste("Expected -log"[10],"(", plain(P),")")),
       y = expression(paste("Observed -log"[10],"(", plain(P),")"))) +
  theme_minimal() +
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=14)) +
  annotate("text", x=1.75, y=8, 
           label= paste("lambda = ",round_half_up(lambda_calculated, 4)),  
           color="grey40", size = 7, fontface="bold") 


options(bitmapType='cairo')
jpeg("B2.jpeg", width = 10, height = 10, units = 'in', res = 300)
print(qqplot)
dev.off()


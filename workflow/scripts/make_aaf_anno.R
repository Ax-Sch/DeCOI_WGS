library(data.table)
vars<-fread(snakemake@input[[1]], fill=TRUE, header = TRUE, na.strings = ".")
# vars<-fread("results/for_RVAS/annotated_split_vep_1.tsv", fill=TRUE, header = TRUE, na.strings = ".")

vars<-vars[, c("chrom","pos","ref","alt"):=tstrsplit(ID, ":", fixed=TRUE)]

###### MASKS ######
vars<-vars[BIOTYPE == "protein_coding"]

vars<-vars[,REVEL_score:=as.numeric(REVEL_score)]
                   

vars<-vars[, ':=' (revel03=!is.na(REVEL_score) & REVEL_score>0.3, 
                   revel05=!is.na(REVEL_score) & REVEL_score>0.5, 
                   revel07=!is.na(REVEL_score) & REVEL_score>0.7,
                   missense=like(Consequence, "missense_variant"),
                   syn=like(Consequence, "synonymous_variant"),
                   pLOF=(IMPACT=="HIGH"),
                   UTR5=like(Consequence, "5_prime_UTR_variant"),
                   UTR3=like(Consequence, "3_prime_UTR_variant"),
		   promoter=like(Consequence,"upstream_gene_variant") & TSSDistance <= 1000,
                   enhancer=like(Consequence, "upstream_gene_variant")& TSSDistance <= 50000 & !is.na(DHS),
                   CADD10= !(is.na(CADD_PHRED) | CADD_PHRED<10) ) ]

vars<-vars[, ':=' (revel03_mis=(revel03 & missense),
                   revel05_mis=(revel05 & missense), 
                   revel07_mis=(revel07 & missense), 
                   moderate_non_missense=(!missense & (IMPACT=="MODERATE")))]

vars<-vars[, by=c("ID", "Gene"), mask_anno:=ifelse(sum(pLOF)>0,"pLoF", 
                                            ifelse(sum(revel07_mis)>0, "missense.revel07",
                                            ifelse(sum(revel05_mis)>0, "missense.revel05",
                                            ifelse(sum(revel03_mis)>0, "missense.revel03",
                                            ifelse(sum(missense | moderate_non_missense)>0, "moderate",
                                            ifelse(sum(syn)>0, "synonymous", 
                                            ifelse(sum(UTR5)>0 & sum(CADD10)>0, "UTR5_CADD", 
                                            ifelse(sum(UTR5)>0, "UTR5",
                                            ifelse(sum(UTR3)>0 & sum(CADD10)>0, "UTR3_CADD",
                                            ifelse(sum(UTR3)>0, "UTR3", 
                                            ifelse(sum(promoter)>0 & sum(CADD10)>0, "promoter_CADD", 
                                            ifelse(sum(promoter)>0, "promoter", 
                                            ifelse(sum(enhancer)>0 & sum(CADD10)>0, "enhancer_CADD", 
                                            ifelse(sum(enhancer)>0, "enhancer",
                                                   "not_relevant")
                                                                        ) ) ) ) ) ) ) ) ) ) ) ) )
]

anno_list<-rbindlist(list(vars[, c("ID", "Gene", "mask_anno")]), use.names = FALSE)

###### AF ######
vars=vars[, maxAF:=pmax(as.numeric(gnomAD_ex_AF), as.numeric(gnomAD_ge_AF), 0,
                        na.rm = TRUE)]

vars=vars[, minAF:=pmin(as.numeric(gnomAD_ex_AF), as.numeric(gnomAD_ge_AF), 1,
                        na.rm = TRUE)]

vars=vars[, maxAC:=pmax(as.numeric(gnomAD_ex_AC), as.numeric(gnomAD_ge_AC), 0, na.rm = TRUE)]

vars=vars[, maxHom:=pmax(as.numeric(gnomAD_ex_nhomalt), as.numeric(gnomAD_ge_nhomalt), 0, na.rm = TRUE)]


#vars=vars[,AF_category:=ifelse( (maxAF < 0.005) & (maxAC > 9 | maxHom > 0 ), 0.005, maxAF) ]
vars=vars[,AF_category:=maxAF]

###### SET LIST #######
collapsed_vars<-vars[,concats:=paste(ID, collapse=","), by=Gene]
collapsed_vars<-unique(collapsed_vars[,c("Gene","chrom","concats")])
collapsed_vars<-collapsed_vars[,pos:=1:nrow(collapsed_vars)]
collapsed_vars<-collapsed_vars[Gene!="" & !is.na(Gene), ]


fwrite(unique(anno_list), file = snakemake@output[[1]], col.names = FALSE, sep=" ")
fwrite(unique(vars[, c("ID", "AF_category")]), file = snakemake@output[[2]], col.names = FALSE, sep=" ")
fwrite((collapsed_vars[,c("Gene","chrom","pos","concats")]), file = snakemake@output[[3]], col.names = FALSE, sep=" ")

fwrite(vars[, -"concats"], file = snakemake@output[[4]], col.names = TRUE, sep="\t")



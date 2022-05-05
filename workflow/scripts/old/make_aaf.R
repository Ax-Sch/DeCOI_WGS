library(data.table)
path_annot<-"/media/axel/Dateien/Arbeit_Gen/COVID_WGS/WGS_GWAS_pipeline/data/rare_chrom_annotated/variants_chr18.tsv.gz"
path_aaf<-"/media/axel/Dateien/Arbeit_Gen/COVID_WGS/WGS_GWAS_pipeline/data/rare_chrom_annotated/variants_chr18.frqx.gz"

path_annot<-snakemake@input[["pathAnnot_tsv"]]
path_aaf<-snakemake@input[["pathAAF_tsv"]]

aafs<-fread(path_aaf)
annots<-as.data.table(read.table(path_annot, header = TRUE, comment.char="#", na.strings = "-"))
setnames(annots, "Uploaded_variation", "ID")

annots[, c("CHROM", "POS") := tstrsplit(Location, ":" , fixed=TRUE)]

annots_aafs<-merge(annots, aafs, by.x = "ID", by.y = "SNP")

###### MASKS ######
vars<-annots_aafs[, ':=' (pph1=like(Polyphen2_HDIV_pred,"D"), 
                  pph2=like(Polyphen2_HVAR_pred,"D"), 
                  lrt=like(LRT_pred,"D"), 
                  mt=like(MutationTaster_pred,"D") | like(MutationTaster_pred,"A"), 
                  sif=like(SIFT_pred,"D"),
                  missense=like(Consequence, "missense_variant"),
                  syn=like(Consequence, "synonymous_variant"),
                  pLOF=(IMPACT=="HIGH"))]

vars<-vars[, ':=' (five_path=(pph1 & pph2 & lrt & mt & sif & missense), 
             one_path=((pph1 | pph2 | lrt | mt | sif)& missense),
           moderate_non_missense=(!missense & (IMPACT=="MODERATE")))]

vars<-vars[, by=c("ID", "Gene"), mask_anno:=ifelse(sum(pLOF)>0,"pLoF", 
             ifelse(sum(five_path)>0, "missense.5in5",
             ifelse(sum(one_path)>0, "missense.1in5",
             ifelse(sum(moderate_non_missense)>0, "moderate.non.missense",
             ifelse(sum(syn)>0, "synonymous", "not_relevant")
             )  ) ) )
            ]

vars<-vars[ mask_anno!="not_relevant"]
vars<-vars[, M0:=syn]
vars<-vars[, M1:=pLOF]
vars<-vars[, M3:=(pLOF| moderate_non_missense | five_path)]
vars<-vars[, M4:=(M3 | one_path)]

###### AF ######

vars=vars[, maxAF:=pmax(AFR_AF, AMR_AF, EAS_AF, EUR_AF, SAS_AF, AA_AF, EA_AF,
                        gnomAD_AFR_AF, gnomAD_AMR_AF, gnomAD_ASJ_AF, gnomAD_EAS_AF, gnomAD_FIN_AF, gnomAD_NFE_AF, gnomAD_OTH_AF, gnomAD_SAS_AF, gnomAD_ge_AF,
                        na.rm = TRUE)]

vars=vars[, minAF:=pmin(AFR_AF, AMR_AF, EAS_AF, EUR_AF, SAS_AF, AA_AF, EA_AF,
                        gnomAD_AFR_AF, gnomAD_AMR_AF, gnomAD_ASJ_AF, gnomAD_EAS_AF, gnomAD_FIN_AF, gnomAD_NFE_AF, gnomAD_OTH_AF, gnomAD_SAS_AF, gnomAD_ge_AF,maxAF,
                        na.rm = TRUE)]


vars=vars[,AF_category:=ifelse( (maxAF<0.001 | minAF>0.999 | is.na(maxAF)), 0.0005,
                        ifelse( (maxAF<0.01 | minAF>0.99), 0.005, 10)) ]

vars[, max_ac:=pmax(0,gnomAD_ge_AC,na.rm=T)]

vars=vars[,AF_category:=ifelse(max_ac>=8, 0.005, AF_category) ]

####### ANNO LIST #######
anno_list<-rbindlist(list(vars[, c("ID", "Gene", "mask_anno")]), use.names = FALSE)

###### SET LIST #######
collapsed_vars<-vars[,concats:=paste(ID, collapse=","), by=Gene]
collapsed_vars<-unique(collapsed_vars[,c("Gene","CHROM","concats")])
collapsed_vars<-collapsed_vars[,pseudo_pos:=1:nrow(collapsed_vars)]


fwrite(unique(anno_list), file = snakemake@output[["regenie_anno_file"]], col.names = FALSE, sep=" ")
fwrite(unique(vars[, c("ID", "AF_category")]), file = snakemake@output[["regenie_aaf_file"]], col.names = FALSE, sep=" ")
fwrite((collapsed_vars[,c("Gene","CHROM","pseudo_pos","concats")]), file = snakemake@output[["regenie_set_list"]], col.names = FALSE, sep=" ")
fwrite(vars, file = snakemake@output[["variables_categorized"]], col.names = TRUE, sep="\t")

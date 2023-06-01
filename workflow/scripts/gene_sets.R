if (!require("GenomicRanges", quietly = TRUE) & (!require("rtracklayer", quietly = TRUE)) ){
  install.packages("BiocManager", repos='http://cran.us.r-project.org')
  BiocManager::install("GenomicRanges")
  BiocManager::install("rtracklayer")
}

library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)


gene_set_path<-"resources/gene_sets.tsv"
gene_set_path<-snakemake@input[[1]]
bed_path<-"resources/genomic_ranges.bed"
bed_path<-snakemake@input[[2]]
aaf_path<-"results/all_vars_for_RVAS/data.anno.aaf.file.txt"
aaf_path<-snakemake@input[[3]]
anno_path<-"/home/aschmidt/Arbeit_Gen/COVID_WGS/COVID_annotation/results/for_RVAS/all_contigs_anno.file.txt"
anno_path<-snakemake@input[[4]]
bim_path<-"results/all_vars_for_RVAS/eurs.bim"
bim_path<-snakemake@input[[5]]
#order_csq<-"resources/genomic_ranges.bed"
#order_csq<-snakemake@input[[6]]

new_anno_path<-"set.annos.tsv"
new_anno_path<-snakemake@output[[1]]
new_aafs_path<-"aafs.tsv"
new_aafs_path<-snakemake@output[[2]]
new_sets_path<-"sets.tsv"
new_sets_path<-snakemake@output[[3]]
relevant_variants_path<-"relevant_variants.tsv"
relevant_variants_path<-snakemake@output[[4]]
new_bim_path<-"new_bim.bim"
new_bim_path<-snakemake@output[[5]]

aafs<-read_delim(aaf_path, delim=" ", col_names = c("ID", "AAF"))
print("aafs:")
head(aafs)
annos<-read_delim(anno_path, delim=" ", col_names = c("ID", "gene_id", "annot"))
print("annos")
head(annos)

bim<-fread(bim_path, sep ="\t", col.names = c("CHR", "ID", "cM","BP","A1","A2"))

print("bim")
head(bim)

#bim_test<-head(bim, n=5000000)

bim_new_id<-bim %>%
  mutate(IDnew=paste("1", 1:n(), A1, A2, sep=":"))%>%
  mutate(BPnew=1:n())%>%
  mutate(CHRnew=1)

comb_sources<-bim_new_id %>%
  left_join(aafs, by="ID")%>%
  left_join(annos, by="ID")

print("combsource")
head(comb_sources)
head(comb_sources %>% filter(!is.na(annot)))

# make sets based on gene definitions
gene_sets<-read_tsv(gene_set_path)
new_sets<-tibble()
new_annos<-tibble()
new_aaf<-tibble()

score_conseq<-function(conseq){
  conseq_tbl<-c("pLoF",
                "missense.revel07",
                "missense.revel05",
                "missense.revel03",
                "moderate",
                "synonymous",
                "UTR5_CADD",
                "UTR5",
		            "UTR3_CADD",
                "UTR3",
                "promoter_CADD",
                "promoter",
                "enhancer_CADD",
                "enhancer",
                "regionBased")
  
  
conseq_rank<-c()
for (single_conseq in conseq){
  single_rank<-which(conseq_tbl %in% single_conseq)
  single_rank<-ifelse(is.na(single_rank), 99, single_rank)
  conseq_rank<-c(conseq_rank, single_rank)
}
conseq_rank<-
return(conseq_rank)
}


i=1
for (single_gene_set in unique(gene_sets$set)){
  rel_genes<-gene_sets %>%
    filter(set==single_gene_set)
  
  tmp_annos<-comb_sources %>%
    filter(gene_id %in% rel_genes$gene_id)%>%
    filter(annot!="not_relevant")%>%
    mutate(gene_id=single_gene_set)
  if (nrow(tmp_annos)==0){
    next()	
  }
  
  # check for variants that are added to the gene set more than once, take the most severe conseq
  mltpl_occuring_tmp<-tmp_annos %>% 
    add_count(IDnew)%>%
    filter(n>1)%>%
    select(-n)
  
  if (nrow(mltpl_occuring_tmp)>0){
    tmp_annos<-tmp_annos %>%
      filter(!IDnew %in% mltpl_occuring_tmp$IDnew)
    
    for (IDnew_tmp in unique(mltpl_occuring_tmp$IDnew)){
      mltpl_occuring_single<-mltpl_occuring_tmp %>%
        filter(IDnew==IDnew_tmp)%>%
        mutate(eff_score=score_conseq(annot))%>%
        arrange(eff_score) 
      
      print(IDnew_tmp)
      print(mltpl_occuring_single)
      
      mltpl_occuring_single<-mltpl_occuring_single %>%
        select(-eff_score)
      
      tmp_annos<-rbind(tmp_annos, mltpl_occuring_single[1,])
    }
  }
  
  new_annos<-rbind(new_annos, tmp_annos) 
  tmp_IDs<-paste(tmp_annos$IDnew, collapse=",")
  tmp_set<-tibble(gene_id=single_gene_set, chr="1", pos=i, IDs=tmp_IDs)
  new_sets<-rbind(new_sets, tmp_set)
  i=i+1 
}



#### 

# make sets based on bed file

region_sets<-import(bed_path)
head(region_sets)
region_sets_names<-str_split(region_sets$name, ";", simplify=TRUE)
region_sets_names<-as_tibble(region_sets_names)
colnames(region_sets_names)<- c("set", "gene")
set_names<-unique(region_sets_names$set)

sources_w_AF<-comb_sources %>%
  filter(!is.na(AAF))%>%
  mutate(pos_start=BP,
         pos_end=BP+nchar(A1))

head(sources_w_AF)
chroms<-paste0("chr",sources_w_AF$CHR)
chroms<-str_replace(chroms, "chr23", "chrX")

annos_granges <- GRanges(
  seqnames = chroms,
  ranges = IRanges(start=sources_w_AF$pos_start, end = sources_w_AF$pos_end) )
head(annos_granges)
for (set_name in set_names){
  
  tmp_region<-region_sets[region_sets_names$set==set_name]
  tmp_intersect<-GenomicRanges::intersect(annos_granges, tmp_region)
  
  tmp_chr_pos_inReg <-sources_w_AF[annos_granges %in% tmp_intersect,]%>%
    distinct(ID, .keep_all=TRUE)%>%
    mutate(gene_id=set_name,
           annot="regionBased")
  if (nrow(tmp_chr_pos_inReg)==0){
    next()
  }

  new_annos<-rbind(new_annos, tmp_chr_pos_inReg, fill=TRUE)
  
  tmp_IDs<-paste(tmp_chr_pos_inReg$IDnew, collapse=",")
  tmp_set<-tibble(gene_id=set_name, chr="1", pos=i, IDs=tmp_IDs)
  new_sets<-rbind(new_sets, tmp_set)
  i=i+1 
}


####

new_annos_export<-new_annos %>%
 distinct(IDnew, gene_id, annot)

new_aaf<-comb_sources %>% distinct(IDnew, AAF) %>%
  filter(!is.na(AAF))

relevant_variants<-unique(new_annos_export$ID)

write_delim(x=new_annos_export,
            file = new_anno_path,
            col_names = FALSE )


write_delim(x=new_aaf,
            file = new_aafs_path,
            col_names = FALSE )

write_delim(x=new_sets,
            file = new_sets_path,
            col_names = FALSE )

write(x = relevant_variants,
      file=relevant_variants_path)

bim_new_id<-bim_new_id%>%
  select(CHRnew,IDnew,cM,BPnew,A1,A2)%>%
  relocate(CHRnew,IDnew,cM,BPnew,A1,A2)

write_tsv(x=bim_new_id,
            file = new_bim_path,
            col_names = FALSE )


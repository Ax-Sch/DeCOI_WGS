import sys
import os
os.environ["LC_ALL"] = "en_US.UTF-8"
import hail as hl

# give the following as arguments: input[0] params["tmp_dir"] output[0]
in_mt=sys.argv[1] # input[0]
tmp_dir=sys.argv[2] # params["tmp_dir"]
out_tsv=sys.argv[3] #output[0]

# initialize hail
hl.init(tmp_dir=tmp_dir)
hl._set_flags(no_whole_stage_codegen='1')

mtAll = hl.read_matrix_table(in_mt)
mtAll=mtAll.drop("PL")
mtAll=hl.sample_qc(mtAll)
mtAll=hl.variant_qc(mtAll)

mt_high_qual = hl.split_multi_hts(mtAll, permit_shuffle = True)
mt_high_qual= mt_high_qual.annotate_entries(AB = (mt_high_qual.AD[1] / hl.sum(mt_high_qual.AD) ))

mt_high_qual = mt_high_qual.filter_entries( (mt_high_qual.GQ>=20) &
	(mt_high_qual.DP >= 8) &
	((mt_high_qual.GT.is_hom_ref() & (mt_high_qual.AB <= 0.1)) |
	(mt_high_qual.GT.is_het() & (mt_high_qual.AB >= 0.25) & (mt_high_qual.AB <= 0.75)) |
	(mt_high_qual.GT.is_hom_var() & (mt_high_qual.AB >= 0.9))))
	
mt_high_qual = mt_high_qual.filter_rows(
	(hl.len(mt_high_qual.alleles) == 2) &
	hl.is_snp(mt_high_qual.alleles[0], mt_high_qual.alleles[1]) &
	(mt_high_qual.variant_qc.call_rate > 0.95) )
	
mt_sex_imputation = mt_high_qual.filter_rows(mt_high_qual.row.locus.in_x_nonpar()==True)
imputed_sex = hl.impute_sex(mt_sex_imputation.GT,aaf_threshold=0.005, female_threshold=0.2, male_threshold=0.8)

mt_high_qual = mt_high_qual.annotate_cols(impute_sex = imputed_sex[mt_high_qual.s])
mt_high_qual = mt_high_qual.filter_rows(mt_high_qual.row.locus.in_autosome()==True)
mt_high_qual = hl.sample_qc(mt_high_qual, name="sample_QC_SNP_AUT_90Call")
mt_high_qual.write(tmp_dir+"hail_gather_data_HC", overwrite=True)
mt_high_qual = hl.read_matrix_table(tmp_dir+"hail_gather_data_HC")
mt_high_qual.col.export(out_tsv)



import sys
import os
os.environ["LC_ALL"] = "en_US.UTF-8"
import hail as hl
import pandas as pd

# give the following as arguments: input[0] tmp_dir output[0]
mt_path=sys.argv[1] # input[0]
tmp_dir=sys.argv[2] # tmp_dir
out_tsv=sys.argv[3]
out_plink=sys.argv[4]
# initialize hail

hl.init(tmp_dir=tmp_dir)
hl._set_flags(no_whole_stage_codegen='1')

mtAll = hl.read_matrix_table(mt_path)
mtAll = hl.variant_qc(mtAll)
mtAll = hl.sample_qc(mtAll)

mtHQ= mtAll.annotate_entries(AB = (mtAll.AD[1] / hl.sum(mtAll.AD) ))

mtHQ = mtHQ.filter_entries( (mtHQ.GQ>=20) &
	(mtHQ.DP >= 8) &
	((mtHQ.GT.is_hom_ref() & (mtHQ.AB <= 0.1)) |
	(mtHQ.GT.is_het() & (mtHQ.AB >= 0.25) & (mtHQ.AB <= 0.75)) |
	(mtHQ.GT.is_hom_var() & (mtHQ.AB >= 0.9))))

mtHQ = hl.variant_qc(mtHQ)

mtHQ = mtHQ.filter_rows(
	(hl.len(mtHQ.alleles) == 2) &
	hl.is_snp(mtHQ.alleles[0], mtHQ.alleles[1]) &
	(mtHQ.variant_qc.call_rate > 0.95) )

mtHQ = mtHQ.filter_rows((mtHQ.variant_qc.AF[1] <= 0.99) & 
	(mtHQ.variant_qc.AF[1] >= 0.01))

mtHQ.write(tmp_dir+"QC_95CR", overwrite=True)
mtHQ = hl.read_matrix_table(tmp_dir+"QC_95CR")

hl.export_plink(mtHQ,out_plink,ind_id = mtHQ.s, fam_id = mtHQ.s, is_female=mtHQ.keep_samples.is_female)

# here the sex imputation was placed

mtHQ_auto = mtHQ.filter_rows(mtHQ.row.locus.in_autosome() == True)
mtHQ_auto.write(tmp_dir+"QC_AUT_95CR", overwrite=True)
mtHQ_auto = hl.read_matrix_table(tmp_dir+"QC_AUT_95CR")
mtHQ_auto = hl.sample_qc(mtHQ_auto, name="QC_AUT_95CR")

#export
mtHQ_Table=mtHQ_auto.cols()
mtHQ_Table_pandas=mtHQ_Table.to_pandas()
mtHQ_Table_pandas.to_csv(out_tsv, sep ='\t')




#mt_sex_imputation = mtHQ.filter_rows(mtHQ.row.locus.in_x_nonpar()==True)
#mt_sex_imputation.write(tmp_dir+"mt_sex", overwrite=True)
#mt_sex_imputation = hl.read_matrix_table(tmp_dir+"mt_sex")
#imputed_sex = hl.impute_sex(mt_sex_imputation.GT, aaf_threshold=0.01, female_threshold=0.2, male_threshold=0.8)
#mtHQ = mtHQ.annotate_cols(impute_sex = imputed_sex[mtHQ.s])

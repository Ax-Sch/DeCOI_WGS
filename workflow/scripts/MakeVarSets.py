import sys
import os
import pandas as pd
os.environ["LC_ALL"] = "en_US.UTF-8"
import hail as hl

tmp_dir=sys.argv[1] # params["tmp_dir"]
in_hail_dir=sys.argv[2] 
input_vcf=sys.argv[3]
individual_list=sys.argv[4]
out_trunk=sys.argv[5] # output[0]
out_Commontrunk=sys.argv[6]

# initialize hail
hl.init(tmp_dir=tmp_dir)
hl._set_flags(no_whole_stage_codegen='1')


mtAll = hl.read_matrix_table(in_hail_dir)

to_keep = hl.import_table(individual_list, 
      impute=True, 
      key='s') # impute=True >> type imputation


mtAll = mtAll.annotate_cols(keep_samplesQC = to_keep[mtAll.s])
mtAll = mtAll.filter_cols(hl.is_defined(mtAll.col.keep_samplesQC))

mtAll = mtAll.filter_rows(hl.len(mtAll.alleles) == 2)

mtAll= mtAll.annotate_entries(AB = (mtAll.AD[1] / hl.sum(mtAll.AD) ))
mtAll = mtAll.filter_entries(
	(mtAll.locus.in_autosome_or_par() | (mtAll.keep_samples.is_female==True)  | (hl.is_missing(mtAll.keep_samples.is_female)) ) & 
	((mtAll.GQ>=20) &
	(mtAll.DP >= 8) &
	((mtAll.GT.is_hom_ref() & (mtAll.AB <= 0.1)) |
	(mtAll.GT.is_het() & (mtAll.AB >= 0.25) & (mtAll.AB <= 0.75)) |
	(mtAll.GT.is_hom_var() & (mtAll.AB >= 0.9))))  |
	(mtAll.locus.in_x_nonpar() & (mtAll.keep_samples.is_female==False))  &
	((mtAll.GQ>=20) &
	(mtAll.DP >= 4) &
	(mtAll.GT.is_hom_ref() & (mtAll.AB <= 0.1)) |
	(mtAll.GT.is_hom_var() & (mtAll.AB >= 0.9))) )

mtAll=hl.variant_qc(mtAll)


mtAll.write(tmp_dir+"pre_call_rate", overwrite=True)
mtAll = hl.read_matrix_table(tmp_dir+"pre_call_rate")

#mtAll.variant_qc.call_rate.export()
mtAll_downsampled=mtAll.sample_rows(0.005)
mtAll_downsampled.write(tmp_dir+"pre_call_rate_downS", overwrite=True)
mtAll_downsampled = hl.read_matrix_table(tmp_dir+"pre_call_rate_downS")

mt_rows_table=mtAll_downsampled.rows()
mt_rows_table_flat=mt_rows_table.flatten()
mt_rows_table_flat.export(out_trunk+"_row_info.tsv.gz")

mtAll = mtAll.filter_rows(mtAll.variant_qc.call_rate > 0.95)

mtAll.write(out_trunk +"_mt", overwrite=True)
mtAll = hl.read_matrix_table(out_trunk + "_mt")

hl.export_plink(mtAll,out_trunk,ind_id = mtAll.s, fam_id = mtAll.s, is_female=mtAll.keep_samples.is_female)
metaData = hl.get_vcf_metadata(input_vcf)
hl.export_vcf(mtAll, out_trunk+".vcf.bgz", metadata=metaData, tabix=True)


# common Vars:
mtAll = mtAll.filter_rows((mtAll.variant_qc.AF[1] <= 0.999) & (mtAll.variant_qc.AF[1] >= 0.001))
mtAll = mtAll.filter_rows(mtAll.variant_qc.call_rate > 0.98)

mtAll.write(out_trunk +"_mtCommon", overwrite=True)
mtAll = hl.read_matrix_table(out_trunk + "_mtCommon")

hl.export_plink(mtAll,out_Commontrunk,ind_id = mtAll.s, fam_id = mtAll.s, is_female=mtAll.keep_samples.is_female)


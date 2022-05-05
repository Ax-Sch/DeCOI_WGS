import sys
import os
os.environ["LC_ALL"] = "en_US.UTF-8"
import hail as hl

tmp_dir=sys.argv[1] # params["tmp_dir"]
in_hail_dir=sys.argv[2] #
individual_list=sys.argv[3]
out_trunk=sys.argv[4] # output[0]

# initialize hail
hl.init(tmp_dir=tmp_dir)
hl._set_flags(no_whole_stage_codegen='1')

mtAll = hl.read_matrix_table(in_hail_dir)
to_keep = hl.import_table(individual_list, 
      impute=True, 
      key='s') # impute=True >> type imputation


mtAll = mtAll.annotate_cols(keep_samples = to_keep[mtAll.s])
mtAll = mtAll.filter_cols(hl.is_defined(mtAll.col.keep_samples))

mtAll = mtAll.filter_rows(hl.len(mtAll.alleles) == 2)

mtAll= mtAll.annotate_entries(AB = (mtAll.AD[1] / hl.sum(mtAll.AD) ))
mtAll = mtAll.filter_entries((mtAll.GQ>=20) &
	(mtAll.DP >= 8) &
	((mtAll.GT.is_hom_ref() & (mtAll.AB <= 0.1)) |
	(mtAll.GT.is_het() & (mtAll.AB >= 0.25) & (mtAll.AB <= 0.75)) |
	(mtAll.GT.is_hom_var() & (mtAll.AB >= 0.9))))

mtAll=hl.variant_qc(mtAll)
mtAll = mtAll.filter_rows((mtAll.variant_qc.AF[1] <= 0.999) & (mtAll.variant_qc.AF[1] >= 0.001))

mtAll.write(tmp_dir+"pre_call_rate", overwrite=True)
mtAll = hl.read_matrix_table(tmp_dir+"pre_call_rate")

mtAll = mtAll.filter_rows(mtAll.variant_qc.call_rate > 0.90)
print(mtAll.count())
mtAll = mtAll.filter_rows(mtAll.variant_qc.call_rate > 0.95)
print(mtAll.count())
mtAll = mtAll.filter_rows(mtAll.variant_qc.call_rate > 0.98)
print(mtAll.count())


mtAll.write(out_trunk +"_mt", overwrite=True)
mtAll = hl.read_matrix_table(out_trunk + "_mt")

hl.export_plink(mtAll,out_trunk,ind_id = mtAll.s, fam_id = mtAll.s)

import sys
import os
import csv
os.environ["LC_ALL"] = "en_US.UTF-8"
import hail as hl

# give the following as arguments: input[0] params["tmp_dir"] output[0]
tmp_dir=sys.argv[1] # params["tmp_dir"]
in_table=sys.argv[2] # input[1]
in_hail_dir=sys.argv[3] # input[2]
out_mt=sys.argv[4] # output[0]

out_dir_list=out_mt.split("/")[0:-1]
out_path="/".join(out_dir_list)


# initialize hail
hl.init(tmp_dir=tmp_dir)
hl._set_flags(no_whole_stage_codegen='1')

sampleQC = hl.import_table(in_table, 
      impute=True, 
      key='s') # impute=True >> type imputation

mtAll = hl.read_matrix_table(in_hail_dir)
mtAll = mtAll.drop("PL")

consecutive_sam_var_count=[]
consecutive_sam_var_count.append("var/samples before filtering: " + str(mtAll.count()))

# Biallelic only, just samples that passed initial QC
mtAll = hl.split_multi_hts(mtAll, permit_shuffle = True)
#mtAll = mtAll.filter_rows(hl.len(mtAll.alleles) == 2)
mtAll = mtAll.annotate_cols(remove_samples = sampleQC[mtAll.s])
mtAll = mtAll.filter_cols(hl.is_defined(mtAll.col.remove_samples))

# Remove bad calls
mtAll= mtAll.annotate_entries(AB = (mtAll.AD[1] / hl.sum(mtAll.AD) ))
mtAll=hl.variant_qc(mtAll)
mtAll = mtAll.filter_entries( (mtAll.GQ>=20) &
	(mtAll.DP >= 8) &
	((mtAll.GT.is_hom_ref() & (mtAll.AB <= 0.1)) |
	(mtAll.GT.is_het() & (mtAll.AB >= 0.25) & (mtAll.AB <= 0.75)) |
	(mtAll.GT.is_hom_var() & (mtAll.AB >= 0.9))))

mtAll=hl.variant_qc(mtAll)

mtAll.write(tmp_dir+"hail_generate_plink_tmp2", overwrite=True)
mtAll = hl.read_matrix_table(tmp_dir+"hail_generate_plink_tmp2")

mtAll.variant_qc.call_rate.export(out_path+ "/var_CR_initial.tsv.bgz")
# 95 % Variant CR
mtAll = mtAll.filter_rows(mtAll.variant_qc.call_rate >= 0.95)

consecutive_sam_var_count.append("var/samples after 95% var CR: " + str(mtAll.count()))

mtAll=hl.sample_qc(mtAll)

mtAll.sample_qc.call_rate.export(out_path+ "/sample_CR_initial.tsv.bgz")

# 97 % Sample CR
mtAll = mtAll.filter_cols(mtAll.sample_qc.call_rate >= 0.97) # 0.97: 15418 variants, 877 samples; 0.98: 15300 variants, 859 samples
mtAll=hl.variant_qc(mtAll)

mtAll.write(out_mt, overwrite=True)
mtAll = hl.read_matrix_table(out_mt)

mtAll.variant_qc.call_rate.export(out_path + "/var_CR_afterSampleCR.tsv.bgz")
mtAll.variant_qc.AF.export(out_path + "/AF_CR_afterSampleCR.tsv.bgz")
consecutive_sam_var_count.append("var/samples after 97% sample CR: " + str(mtAll.count()))

with open(out_path + '/consecutive_sample_count.tsv', 'w', newline='') as f_output:
	tsv_output = csv.writer(f_output, delimiter='\n')
	tsv_output.writerow(consecutive_sam_var_count)




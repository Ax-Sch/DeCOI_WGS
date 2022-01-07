import sys
import os
os.environ["LC_ALL"] = "en_US.UTF-8"
import hail as hl

# give the following as arguments: input[0] params["tmp_dir"] output[0]
input_vcf=sys.argv[1] # input[0]
tmp_dir=sys.argv[2] # params["tmp_dir"]
in_table=sys.argv[3] # input[1]
in_hail_dir=sys.argv[4] # input[2]
out_vcf=sys.argv[5] # output[0]

# initialize hail
hl.init(tmp_dir=tmp_dir)
hl._set_flags(no_whole_stage_codegen='1')

sampleQC = hl.import_table(in_table, 
      impute=True, 
      key='s') # impute=True >> type imputation

metaData = hl.get_vcf_metadata(input_vcf)
mtAll = hl.read_matrix_table(in_hail_dir)
mtAll = mtAll.drop("PL")
print("variants / samples before filtering: " + str(mtAll.count()))
mtAll = hl.split_multi_hts(mtAll, permit_shuffle = True)
mtAll = mtAll.annotate_cols(remove_samples = sampleQC[mtAll.s])

mtAll = mtAll.filter_cols(hl.is_defined(mtAll.col.remove_samples))

# set bad calls to NA:
mtAll= mtAll.annotate_entries(AB = (mtAll.AD[1] / hl.sum(mtAll.AD) ))
mtAll=hl.variant_qc(mtAll)
mtAll = mtAll.filter_entries( (mtAll.GQ>=20) &
	(mtAll.DP >= 8) &
	((mtAll.GT.is_hom_ref() & (mtAll.AB <= 0.1)) |
	(mtAll.GT.is_het() & (mtAll.AB >= 0.25) & (mtAll.AB <= 0.75)) |
	(mtAll.GT.is_hom_var() & (mtAll.AB >= 0.9))))

mtAll=hl.variant_qc(mtAll)
# sapply call rate filter:
mtAll = mtAll.filter_rows(mtAll.variant_qc.call_rate >= 0.95)

mtAll.write(tmp_dir+"hail_generate_plink_tmp2", overwrite=True)
mtAll = hl.read_matrix_table(tmp_dir+"hail_generate_plink_tmp2")
print("variants / samples after first variant QC: " + str(mtAll.count()))

mtAll=hl.sample_qc(mtAll)
mtAll = mtAll.filter_cols(mtAll.sample_qc.call_rate >= 0.97) # 0.97: 15418 variants, 877 samples; 0.98: 15300 variants, 859 samples
mtAll=hl.variant_qc(mtAll)

mtAll.write(tmp_dir+"hail_generate_plink_tmp1", overwrite=True)
mtAll = hl.read_matrix_table(tmp_dir+"hail_generate_plink_tmp1")
print("variants / samples after additional sample QC: " + str(mtAll.count()))

mtAll = mtAll.filter_rows((mtAll.variant_qc.AF[1] <= 0.999) & (mtAll.variant_qc.AF[1] >= 0.001))
mtAll = mtAll.filter_rows(mtAll.variant_qc.call_rate >= 0.98)

mtAll.write(tmp_dir+"hail_generate_plink_tmp2", overwrite=True)
mtAll = hl.read_matrix_table(tmp_dir+"hail_generate_plink_tmp2")
print("variants / samples after second variant QC: " + str(mtAll.count()))

mtAll.write(tmp_dir+"hail_generate_plink_OUT", overwrite=True)
mtAll = hl.read_matrix_table(tmp_dir+"hail_generate_plink_OUT")
print("variants / samples after " + str(mtAll.count()))

hl.export_vcf(mtAll, out_vcf, metadata=metaData, tabix=True)

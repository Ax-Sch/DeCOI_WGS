import sys
import os
import csv
os.environ["LC_ALL"] = "en_US.UTF-8"
import hail as hl

# give the following as arguments: input[0] params["tmp_dir"] output[0]
tmp_dir=sys.argv[1] # params["tmp_dir"]
in_hail_dir=sys.argv[2] # 
input_vcf=sys.argv[3] # 
input_UCSCBed=sys.argv[4]
out_vcf=sys.argv[5] # output[0]

out_dir_list=out_vcf.split("/")[0:-1]
out_path="/".join(out_dir_list)


consecutive_sam_var_count=[]


# initialize hail
hl.init(tmp_dir=tmp_dir)
hl._set_flags(no_whole_stage_codegen='1')

mtAll = hl.read_matrix_table(in_hail_dir)
metaData = hl.get_vcf_metadata(input_vcf)

consecutive_sam_var_count.append("var/samples before filtering: " + str(mtAll.count()))

# 98 % Variant CR, 0.1% AF
mtAll = mtAll.filter_rows((mtAll.variant_qc.AF[1] >= 0.99) | (mtAll.variant_qc.AF[1] <= 0.01))
#mtAll = mtAll.filter_rows(mtAll.variant_qc.call_rate >= 0.98)

mtAll.write(tmp_dir+"hail_rare_var_tmp2", overwrite=True)
mtAll = hl.read_matrix_table(tmp_dir+"hail_rare_var_tmp2")
consecutive_sam_var_count.append("var/samples after MAF < 1% AF Filter: " + str(mtAll.count()))
hl.export_vcf(mtAll, out_vcf, metadata=metaData, tabix=True)

with open(out_path + '/consecutive_sample_count.tsv', 'w', newline='') as f_output:
	tsv_output = csv.writer(f_output, delimiter='\n')
	tsv_output.writerow(consecutive_sam_var_count)


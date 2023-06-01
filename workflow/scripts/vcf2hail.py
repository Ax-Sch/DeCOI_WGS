import sys
import os
os.environ["LC_ALL"] = "en_US.UTF-8"
import hail as hl

# give the following as arguments: input[0] params["tmp_dir"] output[0]
input_vcf=sys.argv[1] # input[0]
out_mt=sys.argv[2]
tmp_dir=sys.argv[3]
individual_list=sys.argv[4]

# initialize hail
hl.init(tmp_dir=tmp_dir)
hl._set_flags(no_whole_stage_codegen='1')
hl.import_vcf(input_vcf, min_partitions=4, reference_genome='GRCh38', force_bgz=True, array_elements_required=False, find_replace=('\t1/.:','\t1/1:') ).write(tmp_dir+"rawMT", overwrite=True )

mtAll = hl.read_matrix_table(tmp_dir+"rawMT")

mtAll = mtAll.drop("PL")

to_keep = hl.import_table(individual_list,
	impute=True, 
	key='s') # impute=True >> type imputation

mtAll = mtAll.annotate_cols(keep_samples = to_keep[mtAll.s])
mtAll = mtAll.filter_cols(hl.is_defined(mtAll.col.keep_samples))

mtAll.write(out_mt, overwrite=True)


import sys
import os
os.environ["LC_ALL"] = "en_US.UTF-8"
import hail as hl

# give the following as arguments: input[0] params["tmp_dir"] output[0]
input_vcf=sys.argv[1] # input[0]
out_mt=sys.argv[2]
tmp_dir=sys.argv[3]

# initialize hail
hl.init(tmp_dir=tmp_dir)
hl._set_flags(no_whole_stage_codegen='1')
hl.import_vcf(input_vcf, min_partitions=4, reference_genome='GRCh38', force_bgz=True, array_elements_required=False).write(out_mt, overwrite=True )

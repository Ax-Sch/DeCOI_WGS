import sys
import os
os.environ["LC_ALL"] = "en_US.UTF-8"
import hail as hl

# give the following as arguments: input[0] params["tmp_dir"] output[0]
input_mt=sys.argv[1] # input[0]
tmp_dir=sys.argv[2] # params["tmp_dir"]
out_vcf=sys.argv[3] #output[0]
prs_score=sys.argv[4]
working_dir=sys.argv[5]

print("working directory: " + working_dir)
os.chdir(working_dir)

hl._set_flags(no_whole_stage_codegen='1')
mtAll = hl.read_matrix_table(input_mt)
mtAll = mtAll.drop("PL")
mtAll = hl.split_multi_hts(mtAll, permit_shuffle = True)
t = hl.import_table(prs_score, 
      impute=True, no_header=True) # impute=True >> type imputation
t=t.annotate(rsid="chr"+hl.str(t.f0)+":"+hl.str(t.f2)+":"+hl.str(t.f3)+":"+hl.str(t.f4))
t=t.key_by(t.rsid)
mtAll = mtAll.annotate_rows(in_score = t[mtAll.rsid])
mtAll.count()
mtAll=mtAll.filter_rows(hl.is_defined(mtAll.in_score.f0))

#mtAll.count()
mtAll.write(tmp_dir + "German_variants", overwrite=True)
mtAll=hl.read_matrix_table(tmp_dir + "/German_variants")
hl.export_vcf(mtAll, out_vcf)

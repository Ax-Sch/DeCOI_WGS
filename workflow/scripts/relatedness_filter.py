import sys
import os
os.environ["LC_ALL"] = "en_US.UTF-8"
import hail as hl

# give the following as arguments: input[0] tmp_dir output[0]
GWAS_MT=sys.argv[1] # input[0]
tmp_dir=sys.argv[2] # tmp_dir
pheno_table=sys.argv[3]
out_kinship=sys.argv[4] # 
out_remove=sys.argv[5] 

# initialize hail

hl.init(tmp_dir=tmp_dir)
hl._set_flags(no_whole_stage_codegen='1')

pheno = hl.import_table(pheno_table, 
      impute=True, 
      key='s')

#mtAll = hl.read_matrix_table(GWAS_MT)
#mtAll = mtAll.annotate_cols(pheno = pheno[mtAll.s])

#mtAll.write(tmp_dir+"hail_related_OUT1", overwrite=True)
mtAll = hl.read_matrix_table(tmp_dir+"hail_related_OUT1")
# https://github.com/broadinstitute/gnomad_qc/blob/master/gnomad_qc/v2/sample_qc/joint_sample_qc.py

#eig, scores, _ = hl.hwe_normalized_pca(mtAll.GT, k=10, compute_loadings=False)
#scores.write(tmp_dir+"hail_scores.ht", overwrite=True )

scores = hl.read_table(tmp_dir+"hail_scores.ht")
# this needs fast ssds on worker nodes
relatedness_ht = hl.pc_relate(mtAll.GT, min_individual_maf=0.01, scores_expr=scores[mtAll.col_key].scores,
                                      block_size=4096, min_kinship=0.01, statistics='all')

relatedness_ht.write(tmp_dir+"rel_matrix", overwrite=True )
relatedness_ht = hl.read_table(tmp_dir+"rel_matrix")

pairs = relatedness_ht.filter(relatedness_ht['kin'] > 0.088)

samples = mtAll.cols()
pairs_with_case = pairs.key_by(
	i=hl.struct(id=pairs.i, is_case=samples[pairs.i].pheno.is_case),
	j=hl.struct(id=pairs.j, is_case=samples[pairs.j].pheno.is_case))

def tie_breaker(l, r):
	return (hl.case()
		.when(hl.is_missing(r.is_case) & hl.is_missing(l.is_case), 0)
		.when(hl.is_missing(r.is_case), -1)
		.when(hl.is_missing(l.is_case), 1)
		.when(l.is_case > r.is_case, -1)
		.when(l.is_case < r.is_case, 1)
		.default(0)
		)
	

related_samples_to_remove = hl.maximal_independent_set(
	pairs_with_case.i, pairs_with_case.j, False, tie_breaker)

relatedness_ht_export=relatedness_ht.annotate(IID1=relatedness_ht.i.s,
	IID2=relatedness_ht.j.s)

relatedness_ht_export.export(out_kinship)

related_samples_to_remove=related_samples_to_remove.annotate(IID=related_samples_to_remove.node.id.s)
related_samples_to_remove.export(out_remove)


import sys
import os
os.environ["LC_ALL"] = "en_US.UTF-8"
import hail as hl

# give the following as arguments: input[0] tmp_dir output[0]
input_vcf=sys.argv[1] # input[0]
tmp_dir=sys.argv[2] # tmp_dir
pop_table=sys.argv[3] # input["populations"]
pheno_table=sys.argv[4]
out_vcf1=sys.argv[5] # output["pop_vcf_rel"]
out_vcf2=sys.argv[6] # 
out_dir=sys.argv[7]
population=sys.argv[8]
# initialize hail
hl.init(tmp_dir=tmp_dir)
hl._set_flags(no_whole_stage_codegen='1')



population_list = hl.import_table(pop_table, 
      impute=True, 
      key='Individual.ID')
pheno = hl.import_table(pheno_table, 
      impute=True, 
      key='s') 

#hl.import_vcf(input_vcf, min_partitions=4, reference_genome='GRCh38', force_bgz=True, array_elements_required=False).write(tmp_dir+"hail_related_IN", overwrite=True ) 
metaData = hl.get_vcf_metadata(input_vcf)
mtAll = hl.read_matrix_table(tmp_dir+"hail_related_IN")

mtAll = mtAll.annotate_cols(pop_samples = population_list[mtAll.s])
mtAll = mtAll.annotate_cols(pheno = pheno[mtAll.s])
mtAll = mtAll.filter_cols((mtAll.col.pop_samples.Prediction == population))

#mtAll.write(tmp_dir+"hail_related_OUT1", overwrite=True)
mtAll = hl.read_matrix_table(tmp_dir+"hail_related_OUT1")
#hl.export_vcf(mtAll, out_vcf1, metadata=metaData, tabix=True)
# https://github.com/broadinstitute/gnomad_qc/blob/master/gnomad_qc/v2/sample_qc/joint_sample_qc.py

#RELATEDNESS FILTER
#eig, scores, _ = hl.hwe_normalized_pca(mtAll.GT, k=10, compute_loadings=False)
#scores.write(tmp_dir+"hail_scores.ht", overwrite=True )

scores = hl.read_table(tmp_dir+"hail_scores.ht")
# this needs fast ssds on worker nodes
relatedness_ht = hl.pc_relate(mtAll.GT, min_individual_maf=0.01, scores_expr=scores[mtAll.col_key].scores,
                                      block_size=4096, min_kinship=0.05, statistics='kin2')

relatedness_ht.write(tmp_dir+"rel_matrix", overwrite=True )
relatedness_ht = hl.read_table(tmp_dir+"rel_matrix")

pairs = relatedness_ht.filter(relatedness_ht['kin'] > 0.088)

samples = mtAll.cols()
pairs_with_case = pairs.key_by(
	i=hl.struct(id=pairs.i, is_case=samples[pairs.i].pheno.is_case),
	j=hl.struct(id=pairs.j, is_case=samples[pairs.j].pheno.is_case))

def tie_breaker(l, r):
	return hl.if_else(l.is_case > r.is_case, -1,
		hl.if_else(l.is_case < r.is_case, 1, 0))

related_samples_to_remove = hl.maximal_independent_set(
	pairs_with_case.i, pairs_with_case.j, False, tie_breaker)

related_samples_to_remove.export(out_dir + "hail_SUC_relatives_to_remove.tsv")

mtAll = mtAll.filter_cols(hl.is_defined(
	related_samples_to_remove.key_by(
	s = related_samples_to_remove.node.id.s)[mtAll.col_key]), keep=False)
	
mtAll.write(tmp_dir+"hail_related_OUT2", overwrite=True)
mtAll = hl.read_matrix_table(tmp_dir+"hail_related_OUT2")
#hl.export_vcf(mtAll, out_vcf2, metadata=metaData, tabix=False)
hl.export_plink(mtAll, out_vcf2, ind_id = mtAll.s, fam_id = mtAll.s)

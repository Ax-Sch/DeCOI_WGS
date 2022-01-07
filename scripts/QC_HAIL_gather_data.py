import sys
import os
os.environ["LC_ALL"] = "en_US.UTF-8"
import hail as hl

# give the following as arguments: input[0] params["tmp_dir"] output[0]
input_vcf=sys.argv[1] # input[0]
tmp_dir=sys.argv[2] # params["tmp_dir"]
out_tsv=sys.argv[3] #output[0]
working_dir=sys.argv[4]
print("working directory: " + working_dir)
os.chdir(working_dir)


# initialize hail
hl.init(tmp_dir=tmp_dir)
hl._set_flags(no_whole_stage_codegen='1')
hl.import_vcf(input_vcf, min_partitions=4, reference_genome='GRCh38', force_bgz=True, array_elements_required=False).write(tmp_dir+"hail_gather_data_IN", overwrite=True )
metaData = hl.get_vcf_metadata(input_vcf)
mtAll = hl.read_matrix_table(tmp_dir+"hail_gather_data_IN")
mtAll=mtAll.drop("PL")
mtAll=hl.sample_qc(mtAll)
mtAll=hl.variant_qc(mtAll)
mt_high_qual = mtAll.filter_rows(
	(hl.len(mtAll.alleles) == 2) &
	hl.is_snp(mtAll.alleles[0], mtAll.alleles[1]) &
	(mtAll.variant_qc.call_rate > 0.90) )
mt_sex_imputation = mt_high_qual.filter_rows(mt_high_qual.row.locus.in_x_nonpar()==True)
imputed_sex = hl.impute_sex(mt_sex_imputation.GT,aaf_threshold=0.005, female_threshold=0.5, male_threshold=0.75)
mt_high_qual = mt_high_qual.annotate_cols(impute_sex = imputed_sex[mt_high_qual.s])
mt_high_qual = mt_high_qual.filter_rows(mt_high_qual.row.locus.in_autosome()==True)
mt_high_qual = hl.sample_qc(mt_high_qual, name="sample_QC_SNP_AUT_90Call")
mt_high_qual.write(tmp_dir+"hail_gather_data_HC", overwrite=True)
mt_high_qual = hl.read_matrix_table(tmp_dir+"hail_gather_data_HC")
mt_high_qual.col.export(out_tsv)


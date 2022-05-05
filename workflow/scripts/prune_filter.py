import sys
import os
os.environ["LC_ALL"] = "en_US.UTF-8"
import hail as hl
import pandas as pd

# give the following as arguments: input[0] tmp_dir output[0]
mt_path=sys.argv[1] # input[0]
tmp_dir=sys.argv[2] # tmp_dir
bed_file=sys.argv[3]
out_mt=sys.argv[4] # 
out_plink=sys.argv[5] 
out_tsv=sys.argv[6]

# initialize hail

hl.init(tmp_dir=tmp_dir)
hl._set_flags(no_whole_stage_codegen='1')

mtAll = hl.read_matrix_table(mt_path)
mtAll = mtAll.drop("PL")

mtAll = hl.variant_qc(mtAll)
mtAll = hl.sample_qc(mtAll)


# sex imputation
mt_high_qual = hl.split_multi_hts(mtAll, permit_shuffle = True)
mt_high_qual= mt_high_qual.annotate_entries(AB = (mt_high_qual.AD[1] / hl.sum(mt_high_qual.AD) ))

mt_high_qual = mt_high_qual.filter_entries( (mt_high_qual.GQ>=20) &
	(mt_high_qual.DP >= 8) &
	((mt_high_qual.GT.is_hom_ref() & (mt_high_qual.AB <= 0.1)) |
	(mt_high_qual.GT.is_het() & (mt_high_qual.AB >= 0.25) & (mt_high_qual.AB <= 0.75)) |
	(mt_high_qual.GT.is_hom_var() & (mt_high_qual.AB >= 0.9))))
	
mt_high_qual = mt_high_qual.filter_rows(
	(hl.len(mt_high_qual.alleles) == 2) &
	hl.is_snp(mt_high_qual.alleles[0], mt_high_qual.alleles[1]) &
	(mt_high_qual.variant_qc.call_rate > 0.95) )
	
mt_sex_imputation = mt_high_qual.filter_rows(mt_high_qual.row.locus.in_x_nonpar()==True)
mt_sex_imputation.write(tmp_dir+"mt_HQ", overwrite=True)
mt_sex_imputation = hl.read_matrix_table(tmp_dir+"mt_HQ")

imputed_sex = hl.impute_sex(mt_sex_imputation.GT, aaf_threshold=0.005, female_threshold=0.2, male_threshold=0.8)
mtAll = mtAll.annotate_cols(impute_sex = imputed_sex[mtAll.s])


# filter for HQ set of variants
mtAll = mtAll.filter_rows(hl.len(mtAll.alleles) == 2)
mtAll = mtAll.filter_rows((mtAll.variant_qc.AF[1] <= 0.99) & 
	(mtAll.variant_qc.AF[1] >= 0.01)&
	hl.is_snp(mtAll.alleles[0], mtAll.alleles[1]))
mtAll.write(tmp_dir+"common_vars", overwrite=True)
mtAll = hl.read_matrix_table(tmp_dir+"common_vars")

mtAll= mtAll.annotate_entries(AB = (mtAll.AD[1] / hl.sum(mtAll.AD) ))
mtAll = mtAll.filter_entries((mtAll.GQ>=20) &
	(mtAll.DP >= 8) &
	((mtAll.GT.is_hom_ref() & (mtAll.AB <= 0.1)) |
	(mtAll.GT.is_het() & (mtAll.AB >= 0.25) & (mtAll.AB <= 0.75)) |
	(mtAll.GT.is_hom_var() & (mtAll.AB >= 0.9))))

mtAll=hl.variant_qc(mtAll)
mtAll = mtAll.filter_rows(mtAll.variant_qc.call_rate > 0.99)

mtAll.write(tmp_dir+"common_vars_QC", overwrite=True)
mtAll = hl.read_matrix_table(tmp_dir+"common_vars_QC")

#https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)
interval_table = hl.import_bed(bed_file, reference_genome='GRCh38')
mtAll = mtAll.filter_rows(hl.is_missing(interval_table[mtAll.locus]))
mtAll = mtAll.key_rows_by(locus=mtAll['locus'],
	alleles=mtAll['alleles'])

# LD Pruning does not work: hail.utils.java.HailUserError: Error summary: HailException: locus_windows: 'locus_expr' global position must be in ascending order. 34104982 was not less then or equal to 15885723

#mtAll = mtAll.repartition(1000)
#mtAll.write(tmp_dir+"common_vars_QC_BED", overwrite=True)
#mtAll = hl.read_matrix_table(tmp_dir+"common_vars_QC_BED")

#pruned_variant_table = hl.ld_prune(mtAll.GT, r2=0.2, bp_window_size=500000)
#mtAll = mtAll.filter_rows(
#	hl.is_defined(pruned_variant_table[mtAll.row_key]))

mtAll.write(out_mt, overwrite=True)
mtAll = hl.read_matrix_table(out_mt)

hl.export_plink(dataset=mtAll, output=out_plink, ind_id = mtAll.s, fam_id = mtAll.s)

mtAll = mtAll.filter_rows(mtAll.row.locus.in_autosome() == True)
mtAll = hl.sample_qc(mtAll, name="sample_QC_SNP_AUT_90Call")
mtAll_Table=mtAll.cols()
mtAll_Table_pandas=mtAll_Table.to_pandas()
mtAll_Table_pandas.to_csv(out_tsv, sep ='\t')



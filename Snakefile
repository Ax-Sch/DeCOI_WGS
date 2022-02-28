# Snakemake file to QC WGS-data and to perform a GWAS

#mamba install hail libgomp1 bgenix plink plink2

configfile: "config/config_xcat_test.yaml"

include: "rules/rare_variant_annotation.smk"


rule all:
	input:
		"data/hail_gather_data/for_sample_QC.tsv",
		"data/hail_gather_data/sampleqc.tsv",
		"data/hail_gather_data/dataset_split_cols.tsv",
		"data/hail_gather_data/plot_sampleQC.html",
		"data/normalized/df3_norm.vcf.gz",
		"data/HAIL_GWAS_vcf/HAIL_GWAS_vcf.vcf.bgz",
		expand("data/regenie_pheno/fam_{population}_{phenotypes}.fam", population=config["populations"], phenotypes=config["phenotypes"]),
		expand("data/additional_QC/{population}_{phenotypes}/not_rel_QCed.bim", population=config["populations"], phenotypes=config["phenotypes"]),
		expand("data/PCAcovar/{phenotypes}_{population}_PCA.eigenvec", population=config["populations"], phenotypes=config["phenotypes"]),
		expand("data/regenie_pheno/cov_{phenotypes}_{population}", population=config["populations"], phenotypes=config["phenotypes"]),
		expand("data/regenie/GWAS_{population}_{phenotypes}.regenie", population=config["populations"], phenotypes=config["phenotypes"]),
		expand("data/regenie/GWAS_{population}_{phenotypes}.regenie_LOG10P_manhattan.png", population=config["populations"], phenotypes=config["phenotypes"]),
		"data/PRS/variants_for_prs.vcf.bgz",
		"data/PCA/populations.txt",
		multiext("data/ROH/bedForROHAnalysis_ROH", ".hom", ".hom.indiv", ".hom.summary", ".log"),
		"data/plotPhenoInfo/AgePlot.pdf",
		"data/ROH/test_roh.tsv.gz",
		"data/rare_variant_vcf/rare_variant_vcf.vcf.bgz",
		directory("data/hail_gather_data/hail_gather_data_IN"),
		"data/rare_chrom_annotated/variants_chr21.vcf.gz"
		

		
# normalize the vcf file and annotate it with the ID: CHROM:POS:REF:ALT
rule normalize_vcf_annotate_bcftools:
	input:
		bcf=config["input_vcf"],
		fasta=config["fasta"]
	output:
		vcf_gz="data/normalized/df3_norm.vcf.gz",
		vcf_gz_index="data/normalized/df3_norm.vcf.gz.tbi"
	resources: cpus=1, mem_mb=8000, time_job=10080
	params:
		partition='long'
	shell:
		"""
		#wget -nc ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
		#gzip -d GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
		bcftools filter -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX \
		  {input.bcf} -Ou | \
		  bcftools norm --check-ref w -f {input.fasta} -Ou | \
		  bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz > {output.vcf_gz}
		tabix -p vcf {output.vcf_gz}
		"""

# Just calculate sample QC-data using hail
rule GatherSampleQCData:
	input:
		"data/normalized/df3_norm.vcf.gz"
	output:
		sample_qc_file="data/hail_gather_data/sampleqc.tsv",
		hail_MT=directory("data/hail_gather_data/hail_gather_data_IN")
	resources: cpus=4, mem_mb=10000, time_job=10080, additional=" -x " + config["master_nodes_excluded"]
	params:
		partition='long',
		tmp_dir=config["tmp_folder"],
		path_output="data/hail_gather_data"
	shell:
		"""

		mkdir -p {params.path_output}
		mkdir -p {params.tmp_dir}
		
		hail_python_script="scripts/QC_HAIL_gather_data.py {input} {params.tmp_dir} {output.sample_qc_file} {output.hail_MT}"
		
		if [ {config[cluster]} = "yes" ]; then
		worker_nodes_excluded={config[worker_nodes_excluded]}
		num_workers=5
		source scripts/spark_submit_command.sh
		$spark_submit_command $hail_python_script
		
		else
		python $hail_python_script
		fi
		"""

# reformat the hail sample QC file, generate histograms of all values, generate a list of samples passing the sample QC thresholds 
rule ReformatPlotSampleQC:
	input:
		"data/hail_gather_data/sampleqc.tsv",
	output:
		"data/hail_gather_data/dataset_split_cols.tsv",
		"data/hail_gather_data/plot_sampleQC.html",
		"data/hail_gather_data/for_sample_QC.tsv"
	resources: cpus=1, mem_mb=4000, time_job=720
	params:
		partition='batch',
		path_output="data/hail_gather_data/"
	shell:
		"""
		cp scripts/convert_sampleQC_to_tsv.R {params.path_output}
		cp scripts/plot_sampleQC.R {params.path_output}
		cd {params.path_output}
		
		Rscript convert_sampleQC_to_tsv.R
		Rscript -e 'library(rmarkdown); rmarkdown::render("plot_sampleQC.R", "html_document")'
		rm *.R
		"""


# do sample QC based on the list generated from above, do variant QC, filter for variants >0.1% AF (=AC>2)
rule DoFirstSampleVariantQC:
	input:
		table="data/hail_gather_data/for_sample_QC.tsv",
		hail_MT=directory("data/hail_gather_data/hail_gather_data_IN")
	output:
		QCed_MT=directory("data/first_QC/first_QCed_MT")
	resources: cpus=4, mem_mb=10000, time_job=10080, additional=" -x " + config["master_nodes_excluded"]
	params:
		partition='long',
		tmp_dir=config["tmp_folder"],
		path_output="data/first_QC"
	shell:
		"""
		mkdir -p {params.path_output}
		mkdir -p {params.tmp_dir}
		
		hail_python_script="scripts/DoFirstSampleVariantQC.py {params.tmp_dir} {input.table} {input.hail_MT} {output}"
		
		if [ {config[cluster]} = "yes" ]; then
		worker_nodes_excluded={config[worker_nodes_excluded]}
		num_workers=5
		source scripts/spark_submit_command.sh
		$spark_submit_command $hail_python_script
		
		else
		python $hail_python_script
		fi
		"""

rule get_GWAS_vcf:
	input:
		QCed_MT=directory("data/first_QC/first_QCed_MT"),
		input_vcf="data/normalized/df3_norm.vcf.gz"
	output:
		vcf_for_GWAS="data/HAIL_GWAS_vcf/HAIL_GWAS_vcf.vcf.bgz"
	resources: cpus=4, mem_mb=10000, time_job=10080, additional=" -x " + config["master_nodes_excluded"]
	params:
		partition='long',
		tmp_dir=config["tmp_folder"],
		path_output="data/HAIL_GWAS_vcf"
	shell:
		"""
		mkdir -p {params.path_output}
		mkdir -p {params.tmp_dir}
		
		hail_python_script="scripts/MakeGwasVCF.py {params.tmp_dir} {input.QCed_MT} {input.input_vcf} {output}"
		
		if [ {config[cluster]} = "yes" ]; then
		worker_nodes_excluded={config[worker_nodes_excluded]}
		num_workers=5
		source scripts/spark_submit_command.sh
		$spark_submit_command $hail_python_script
		
		else
		python $hail_python_script
		fi
		"""



rule get_rare_variant_VCF:
	input:
		QCed_MT=directory("data/first_QC/first_QCed_MT"),
		input_vcf="data/normalized/df3_norm.vcf.gz",
		input_UCSC_bed="data/bedfilesForPrefilter/ExonsClinVarSpliceAI.bed"
	output:
		vcf_for_GWAS="data/rare_variant_vcf/rare_variant_vcf.vcf.bgz"
	resources: cpus=4, mem_mb=10000, time_job=10080, additional=" -x " + config["master_nodes_excluded"]
	params:
		partition='long',
		tmp_dir=config["tmp_folder"],
		path_output="data/rare_variant_vcf"
	shell:
		"""
		mkdir -p {params.path_output}
		mkdir -p {params.tmp_dir}
		
		hail_python_script="scripts/makeRareVariantVCF.py {params.tmp_dir} {input.QCed_MT} {input.input_vcf} {input.input_UCSC_bed} {output}"
		
		if [ {config[cluster]} = "yes" ]; then
		worker_nodes_excluded={config[worker_nodes_excluded]}
		num_workers=5
		source scripts/spark_submit_command.sh
		$spark_submit_command $hail_python_script
		
		else
		python $hail_python_script
		fi
		"""


#### Get populations based on the 1000G project


rule download_1000G_sample_info:
	output:
		#"GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
		"data/PCA/20130606_g1k.ped" # contig = chr1, ...
	resources: cpus=1, mem_mb=3000, time_job=720
	params:
		partition='batch'
	shell:
		"""
		cd data/PCA/
		wget -nc ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped
		
		#reference genome  (GRCh38)
		#wget -nc ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
		#gzip -d GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
		"""

rule download_1000G_genotypes:
	output:
		config["location_1000G"]+"/ALL.chr{contig}.shapeit2_integrated_v1a.GRCh38.20181129.phased.bcf",
		config["location_1000G"]+"/ALL.chr{contig}.shapeit2_integrated_v1a.GRCh38.20181129.phased.bcf.tbi"#
	resources: cpus=1, mem_mb=3000, time_job=720
	params:
		partition='batch',
		folder_cont=config["location_1000G"]
	shell:
		"""
		mkdir -p {params.folder_cont}
		cd {params.folder_cont}
		prefix="ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr"
		suffix=".shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"
		wget -nc "$prefix""{wildcards.contig}""$suffix" "$prefix""{wildcards.contig}""$suffix".tbi
		"""


rule prepare_1000G_for_ancestry_PCA:
	input:
		vcf=config["location_1000G"]+"/ALL.chr{contig}.shapeit2_integrated_v1a.GRCh38.20181129.phased.bcf",
		fasta=config["fasta"], #/ceph01/projects/bioinformatics_resources/Genomes/Human/GATK/hg38/Homo_sapiens_assembly38.fasta
		hailvcf="data/HAIL_GWAS_vcf/HAIL_GWAS_vcf.vcf.bgz",
	output:
		vcf="data/1000G/1000G_chr{contig}.bcf",
		variant_list="data/1000G/1000G_chr{contig}.vars",
		bed2="data/1000G/1000G_chr{contig}_pruned.bed"
	resources: cpus=1, mem_mb=18000, time_job=720
#	conda:
#		"envs/bcftools_plink_R.yaml"
	params:
		partition='batch',
		bed1="data/1000G/1000G_chr{contig}",
		bed2="data/1000G/1000G_chr{contig}_pruned"
	shell:
		"""
		if bcftools view -q 0.05:minor {input.vcf} | \
		bcftools norm -m-any --check-ref w -f "{input.fasta}" | \
		bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
		bcftools norm -Ob --rm-dup both \
		> {output.vcf} ; then
		echo "no error"
		fi
		
		bcftools index {output.vcf}
		
		bcftools view -r chr{wildcards.contig} {input.hailvcf} -Ou | \
		bcftools query -f "%ID\\n" > {output.variant_list}
		
		plink --noweb \
		--bcf {output.vcf} \
		--keep-allele-order \
		--vcf-idspace-to _ \
		--allow-extra-chr 0 \
		--split-x b38 no-fail \
		--make-bed \
		--out {params.bed1}
		
		plink --noweb \
		--bfile {params.bed1} \
		--extract {output.variant_list} \
		--maf 0.10 --indep 50 5 1.5 \
		--make-bed \
		--out {params.bed2}
		"""
		

rule merge_data_w_1000G_run_PCA:
	input:
		_1000G_data=expand("data/1000G/1000G_chr{contig}_pruned.bed", contig=config["contigs_wo_X"]), #### !!!!! 
		hailvcf="data/HAIL_GWAS_vcf/HAIL_GWAS_vcf.vcf.bgz",
		ped_file_1000G="data/PCA/20130606_g1k.ped"
	output:
		merge_list="data/1000G/merge_list.txt",
		bim_pca="data/PCA/MergeFullForPCA.bim",
		populations="data/PCA/populations.txt"
	resources: cpus=1, mem_mb=18000, time_job=720
	params:
		partition='batch',
#	conda:
#		"envs/bcftools_plink_R.yaml"
	shell:
		"""
		echo {input._1000G_data} | tr " " "\\n" | sed 's/.bed//g' > {output.merge_list}
		plink --merge-list {output.merge_list} --out data/1000G/Merged
		awk '{{ print $2 }}' data/1000G/Merged.bim > data/1000G/MergeVariants.txt
		
		plink --vcf {input.hailvcf} \
		 --double-id \
		 --noweb \
		 --extract data/1000G/MergeVariants.txt \
		 --make-bed \
		 --out data/1000G/hail_for_ancestry
		 
		printf "data/1000G/Merged\\ndata/1000G/hail_for_ancestry" > data/1000G/ForMergeFull.list
		
		mkdir -p data/PCA
		plink --merge-list data/1000G/ForMergeFull.list --out data/PCA/MergeFullForPCA
		 
		awk '{{ print $1,$2 }}' data/1000G/Merged.fam | awk '$(NF+1) = "1000G"' > data/PCA/clusters.txt
		awk '{{ print $1,$2 }}' data/1000G/hail_for_ancestry.fam | awk '$(NF+1) = "Cohort"' >> data/PCA/clusters.txt
		
		plink --bfile data/PCA/MergeFullForPCA \
		 --pca-cluster-names 1000G \
		 --pca \
		 --within data/PCA/clusters.txt
		
		mv plink.* data/PCA/
		cp -f scripts/populations_PCA.R data/PCA/populations_PCA.R
		cd data/PCA/
		Rscript -e 'library(rmarkdown); rmarkdown::render("populations_PCA.R", "html_document")'

		"""
		
rule ApplyRelatednessFilter:
	input:
		vcf_for_GWAS="data/HAIL_GWAS_vcf/HAIL_GWAS_vcf.vcf.bgz",
		populations="resources/populations.txt",
		pheno_file="data/hail_gather_data/for_sample_QC.tsv"
	output:		
		PLINK_fam="data/pop_vcf/{population}_vcf_not_rel.vcf.bgz.fam",
		bed="data/pop_vcf/{population}_vcf_not_rel.vcf.bgz.bed",
		bim="data/pop_vcf/{population}_vcf_not_rel.vcf.bgz.bim"
	resources: cpus=4, mem_mb=32000, time_job=2880, additional=" --gres localtmp:300G --ntasks=1 -x " + config["master_nodes_excluded"]
	params:
		partition='long',
		out_dir="data/pop_vcf/",
		tmp_dir=config["tmp_folder"] +"{population}/",
		pop_vcf_rel="data/pop_vcf/{population}_vcf_rel.vcf.bgz",
		pop_vcf_not_rel="data/pop_vcf/{population}_vcf_not_rel.vcf.bgz",
	shell:
		"""		

		mkdir -p {params.out_dir}
		mkdir -p {params.tmp_dir}
		
		hail_python_script="scripts/relatedness_filter.py $(pwd)/{input.vcf_for_GWAS} $(pwd)/{params.tmp_dir} $(pwd)/{input.populations} $(pwd)/{input.pheno_file} $(pwd)/{params.pop_vcf_rel} $(pwd)/{params.pop_vcf_not_rel} $(pwd)/{params.out_dir} {wildcards.population}"
		
		if [ {config[cluster]} = "yes" ]; then
		queue="long"
		hours_to_run=48
		worker_nodes_excluded={config[worker_nodes_excluded]}
		num_workers=6
		source scripts/spark_submit_command.sh
		$spark_submit_command $hail_python_script
		
		else
		python $hail_python_script
		fi		
		"""

rule generate_pheno_files_for_GWAS:
	input:
		fam="data/pop_vcf/{population}_vcf_not_rel.vcf.bgz.fam",
		xl_file="config/post_df3_HGI_sample_QC_summary.xlsx"
	output:
		fam="data/regenie_pheno/fam_{population}_{pheno}.fam",
		pheno="data/regenie_pheno/pheno_{population}_{pheno}"
	resources: cpus=1, mem_mb=3000, time_job=720
	params:
		partition='batch',
		path="data/regenie_pheno/"
	shell:
		"""
		mkdir -p {params.path} 
		Rscript scripts/generate_regenie_pheno.R {input.xl_file} {input.fam} {wildcards.pheno} {output.fam} {output.pheno}
		"""


rule additional_GWAS_Specifc_QC:
	input:
		bed="data/pop_vcf/{population}_vcf_not_rel.vcf.bgz.bed",
		bim="data/pop_vcf/{population}_vcf_not_rel.vcf.bgz.bim",
		fam="data/regenie_pheno/fam_{population}_{pheno}.fam",
	output:
		bed="data/additional_QC/{population}_{pheno}/not_rel_QCed.bed",
		bim="data/additional_QC/{population}_{pheno}/not_rel_QCed.bim",
		fam="data/additional_QC/{population}_{pheno}/not_rel_QCed.fam"
	resources: cpus=6, mem_mb=90000, time_job=720, additional=" -w b310eth0 --ntasks=1 "
	params:
		partition='batch',
		pathQC="data/additional_QC/{population}_{pheno}"
	shell:
		"""
		mkdir -p {params.pathQC}
		plink \
		--bed {input.bed} \
		--bim {input.bim} \
		--fam {input.fam} \
		--filter-cases \
		--missing \
		--out {params.pathQC}/Geno05_CR_sex_snp_qc_snpqcCAS
		
		plink \
		--bed {input.bed} \
		--bim {input.bim} \
		--fam {input.fam} \
		--filter-controls \
		--missing \
		--out {params.pathQC}/Geno05_CR_sex_snp_qc_snpqcCON
		
		plink \
		--bed {input.bed} \
		--bim {input.bim} \
		--fam {input.fam} \
		--hardy \
		--out {params.pathQC}/Geno05_CR_sex_snp_qc_snpqcALL
				
		plink \
		--bed {input.bed} \
		--bim {input.bim} \
		--fam {input.fam} \
		--hardy \
		--chr 23 \
		--filter-females \
		--out {params.pathQC}/Geno05_CR_sex_snp_qc_snpqcALL_female

		### R session
		wd=$(pwd)
		cp scripts/additional_QC_cohort.R {params.pathQC}/
		cd {params.pathQC}
		Rscript -e 'library(rmarkdown); rmarkdown::render("additional_QC_cohort.R", "html_document")'
		cd $wd
		
		# 1e. final QC (4675116)
		plink \
		--bed {input.bed} \
		--bim {input.bim} \
		--fam {input.fam} \
		--exclude {params.pathQC}/missingness_hardy_weinberg_filter.txt \
		--make-bed \
		--out {params.pathQC}/not_rel_QCed

		"""

rule generate_PCA_covar_for_GWAS:
	input:
		bim="data/additional_QC/{population}_{pheno}/not_rel_QCed.bim",
		bed="data/additional_QC/{population}_{pheno}/not_rel_QCed.bed",
		fam="data/additional_QC/{population}_{pheno}/not_rel_QCed.fam"
	output:
		"data/PCAcovar/{pheno}_{population}_PCA.eigenvec"
	resources: cpus=1, mem_mb=18000, time_job=720
	params:
		partition='batch',
		folderPCA="data/PCAcovar",
		pathinterim="data/PCAcovar/{pheno}_{population}__interimfiles",
		pathPCA="data/PCAcovar/{pheno}_{population}_PCA"
	shell:
		"""
		mkdir -p {params.folderPCA}
		
		#common
		plink \
		--bed {input.bed} \
		--bim {input.bim} \
		--fam {input.fam} \
		--biallelic-only strict \
		--chr 1-22 \
		--geno 0.05 \
		--snps-only 'just-acgt' \
		--hwe 1E-6 midp \
		--indep-pairwise 50 5 0.05 \
		--keep-allele-order \
		--mac 5 \
		--maf 0.01 \
		--out "{params.pathinterim}"

		plink \
		--bed {input.bed} \
		--bim {input.bim} \
		--fam {input.fam} \
		--extract "{params.pathinterim}.prune.in" \
		--pca 10 \
		--out "{params.pathPCA}"
		rm {params.pathinterim}*
		"""


rule generate_covar_files_GWAS:
	input:
		PCA_cov="data/PCAcovar/{pheno}_{population}_PCA.eigenvec",
		xl_file="config/post_df3_HGI_sample_QC_summary.xlsx"
	output:
		"data/regenie_pheno/cov_{pheno}_{population}"
	resources: cpus=1, mem_mb=3000, time_job=720
	params:
		partition='batch',
	shell:
		"""
		Rscript scripts/generate_covar_file.R {input.xl_file} {input.PCA_cov} {output}
		"""
 
rule Do_GWAS_with_regenie:
	input:
		pathCov="data/regenie_pheno/cov_{pheno}_{population}",
		pathPheno="data/regenie_pheno/pheno_{population}_{pheno}",
		bed="data/additional_QC/{population}_{pheno}/not_rel_QCed.bed",
		bim="data/additional_QC/{population}_{pheno}/not_rel_QCed.bim",
		fam="data/additional_QC/{population}_{pheno}/not_rel_QCed.fam"
	output:
		path_reg="data/regenie/GWAS_{population}_{pheno}.regenie"
	resources: cpus=1, mem_mb=18000, time_job=720
	params:
		partition='batch',
		tmp_dir=config["tmp_folder"] + "/{population}_{pheno}",
		pruned_plink="data/regenie/pheno_{population}_{pheno}_pruned_plink",
		regenie_step1="data/regenie/pheno_{population}_{pheno}_step1",
		regenie_step2="data/regenie/GWAS_{population}"
	shell:
		"""
		infile={input.bed}
		mkdir -p {params.tmp_dir}
		
		# LD pruned variants for regenie step 1
		plink \
		--bed {input.bed} \
		--bim {input.bim} \
		--fam {input.fam} \
		--maf 0.01 \
		--chr 1-22 \
		--geno 0.01 \
		--snps-only just-acgt \
		--hwe 1e-10 \
		--exclude range resources/longrange_LD_regions_hg38_GRCh38.txt \
		--indep-pairwise 50kb 5 0.05 \
		--out {params.pruned_plink}
		
		export LD_LIBRARY_PATH={config[LD_LIBRARY_PATH]}
		
		tools/regenie_v2.2.4.gz_x86_64_Linux_mkl \
		  --step 1 \
		  --bed ${{infile::-4}} \
		  --covarFile {input.pathCov} \
		  --phenoFile {input.pathPheno} \
		  --bt \
		  --loocv \
		  --extract {params.pruned_plink}.prune.in \
		  --bsize 1000 \
		  --out {params.regenie_step1}
		  #--lowmem \
		  #--lowmem-prefix {params.tmp_dir} \
		
		#step 2

		tools/regenie_v2.2.4.gz_x86_64_Linux_mkl \
		  --step 2 \
		  --minMAC 16 \
		  --covarFile {input.pathCov} \
		  --phenoFile {input.pathPheno} \
		  --bed ${{infile::-4}} \
		  --bt \
		  --spa \
		  --write-samples \
		  --pred {params.regenie_step1}_pred.list \
		  --bsize 1000 \
		  --out {params.regenie_step2}
		  #--htp DeCOI \
		  #--firth --approx \
		  #--firth-se \
		  #--maxstep-null 5 \
		  #--maxiter-null 10000 \
		  
		  #rm {params.pruned_plink}* {params.regenie_step1}*
		  """
		  

rule generate_qq_plots:
	input:
		path_reg="data/regenie/GWAS_{population}_{pheno}.regenie"
	output:
		"data/regenie/GWAS_{population}_{pheno}.regenie_LOG10P_manhattan.png" #GWAS_EUR_A2.regenie_plot_Pval_manhattan.png
	resources: cpus=1, mem_mb=3000, time_job=720
	params:
		partition='batch',
		output_string="data/regenie/GWAS_{population}_{pheno}.regenie"
	shell:
		"""
		Rscript scripts/qqplot.R -f {input} -o {params.output_string} -c CHROM -p LOG10P -b GENPOS
		"""



rule Extract_PRS_variants:
	input:
		mt_table=directory(config["tmp_folder"] + "hail_gather_data_IN"),
		prs_file="resources/hg38.PRS.tsv"
	output:		
		vcf_for_PRS="data/PRS/variants_for_prs.vcf.bgz"
	resources: cpus=12, mem_mb=48000, time_job=720, additional=" --gres localtmp:200G --ntasks=1 -x " + config["master_nodes_excluded"]
	params:
		partition='batch',
		out_dir="data/PRS/",
		tmp_dir=config["tmp_folder"] +"PRS/",
	shell:
		"""

		mkdir -p {params.out_dir}
		mkdir -p {params.tmp_dir}
		
		hail_python_script="scripts/hail_get_PRS_variants.py $(pwd)/{input.mt_table} $(pwd)/{params.tmp_dir} $(pwd)/{output} $(pwd)/{input.prs_file} $(pwd)"
		
		if [ {config[cluster]} = "yes" ]; then
		worker_nodes_excluded={config[worker_nodes_excluded]}
		num_workers=5
		source scripts/spark_submit_command.sh
		$spark_submit_command $hail_python_script
		
		else
		python $hail_python_script
		fi
		"""


rule Extract_Rohs_plink:
	input:
		vcf_for_Plink="data/HAIL_GWAS_vcf/HAIL_GWAS_vcf.vcf.bgz",
	output:
		bed_file="data/ROH/bedForROHAnalysis.bed",
		ROH_file=multiext("data/ROH/bedForROHAnalysis_ROH", ".hom", ".hom.indiv", ".hom.summary", ".log")
	resources: cpus=12, mem_mb=48000, time_job=720, additional=" --gres localtmp:200G --ntasks=1 -x " + config["master_nodes_excluded"]
	params:
		partition='batch',
		out_dir="data/ROH/",
		bed_out="data/ROH/bedForROHAnalysis",
		ROH_out="data/ROH/bedForROHAnalysis_ROH"
	shell:
		"""
		mkdir -p {params.out_dir}
		
		plink --vcf {input.vcf_for_Plink} --out {params.bed_out} --make-bed
		plink --bfile {params.bed_out} --out {params.ROH_out} --homozyg
		"""

rule Extract_Rohs_bcftools:
	input:
		vcf_for_GWAS="data/HAIL_GWAS_vcf/HAIL_GWAS_vcf.vcf.bgz",
	output:		
		ROH_file="data/ROH/test_roh.tsv.gz",
	resources: cpus=6, mem_mb=24000, time_job=720, additional=" "
	params:
		partition='batch',
		out_dir="data/ROH/",
	shell:
		"""
		mkdir -p {params.out_dir}


		# the code here could e.g. look like this:
		bcftools roh {input.vcf_for_GWAS} --AF-tag AF -G30 --threads 6 | \
		 gzip > {output.ROH_file}
		"""

rule CompilePhenoInfo:
	input:
		xl_file="config/post_df3_HGI_sample_QC_summary.xlsx",
		qc_file="data/hail_gather_data/for_sample_QC.tsv"
	output:
		AgePlot="data/plotPhenoInfo/AgePlot.pdf",

	resources: cpus=1, mem_mb=6000, time_job=720, additional=" --ntasks=1 "
	params:
		partition='batch',
		out_dir="data/plotPhenoInfo/",
	shell:
		"Rscript scripts/compile_phenotype_plots.R {input} {params.out_dir}"

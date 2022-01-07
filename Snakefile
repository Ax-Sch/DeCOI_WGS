# Snakemake file to QC WGS-data and to perform a GWAS

#mamba install hail libgomp1 bgenix plink plink2

configfile: "config/config_xcat_test.yaml"

rule all:
	input:
		"data/hail_gather_data/for_sample_QC.tsv",
		"data/hail_gather_data/sampleqc.tsv",
		"data/hail_gather_data/dataset_split_cols.tsv",
		"data/hail_gather_data/plot_sampleQC.html",
		"data/normalized/df3_norm.vcf.gz",
		"data/HAIL_GWAS_vcf/HAIL_GWAS_vcf.vcf.bgz",
		#expand("data/pop_vcf/{population}_vcf_rel.vcf.bgz", population=config["populations"]),
		expand("data/regenie_pheno/fam_{population}_{phenotypes}.fam", population=config["populations"], phenotypes=config["phenotypes"]),
		expand("data/additional_QC/{population}_{phenotypes}/not_rel_QCed.bim", population=config["populations"], phenotypes=config["phenotypes"]),
		expand("data/PCAcovar/{phenotypes}_{population}_PCA.eigenvec", population=config["populations"], phenotypes=config["phenotypes"]),
		expand("data/regenie_pheno/cov_{phenotypes}_{population}", population=config["populations"], phenotypes=config["phenotypes"]),
		expand("data/regenie/GWAS_{population}_{phenotypes}.regenie", population=config["populations"], phenotypes=config["phenotypes"]),
		expand("data/regenie/GWAS_{population}_{phenotypes}.regenie_Pval_manhattan.png", population=config["populations"], phenotypes=config["phenotypes"]),
		"data/PRS/variants_for_prs.vcf.gz",

		
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
		  bcftools norm -m -any --check-ref w -f {input.fasta} -Ou | \
		  bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz > {output.vcf_gz}
		tabix -p vcf {output.vcf_gz}
		"""

# Just calculate sample QC-data using hail
rule HAIL_gather_QC_data:
	input:
		"data/normalized/df3_norm.vcf.gz"
	output:
		sample_qc_file="data/hail_gather_data/sampleqc.tsv",
		out_dir=directory(config["tmp_folder"] + "hail_gather_data_IN")
	resources: cpus=4, mem_mb=10000, time_job=10080, additional=" -w b308eth0"
	params:
		partition='long',
		tmp_dir=config["tmp_folder"],
		path_output="data/hail_gather_data"
	shell:
		"""
		# set up spark cluster
		mkdir -p {params.path_output}
		mkdir -p {params.tmp_dir}
		export LC_ALL="en_US.UTF-8"
		
		if [ {config[cluster]} = "yes" ]; then
		
		echo "Setting up Spark Cluster."
		export SPARK_LOCAL_DIRS=$(pwd)/{params.tmp_dir}
		export PATH={config[PATH]}
		export PYTHONPATH={config[PYTHONPATH]}
		HAIL_HOME=$(pip3 show hail | grep Location | awk -F' ' '{{print $2 "/hail"}}')
		export SPARK_CLASSPATH=$HAIL_HOME"/backend/hail-all-spark.jar"		
		echo "Env variables set."
		# start spark master in this job, add one worker that will definitely start
		{config[spark_master_path]}
		{config[spark_worker_path]} "spark://${{SLURMD_NODENAME}}.cluster.imbie:7077" --cores 2
		
		# add addtitional workers, save job-ids and add a bash trap to cancel them, when this process stops
		slurm_id=""
		for i in {{1..20}}
		do
		   slurm_id="$slurm_id "$(sbatch --time=12:00:00 -J spark_worker --ntasks=1 --partition=batch \
		   --mem=84G --cpus-per-task=24 \
		   --chdir=$(pwd) \
		   --wrap "{config[spark_worker_path]} spark://${{SLURMD_NODENAME}}.cluster.imbie:7077; sleep 12h" | awk '{{print $NF}}')
		done
		
		echo "started workers as jobs:"
		echo $slurm_id
		trap "scancel $slurm_id" EXIT
		
		# submit the hail script
		
		spark-submit \
		--jars $HAIL_HOME/backend/hail-all-spark.jar \
		--conf spark.driver.extraClassPath=$HAIL_HOME/backend/hail-all-spark.jar \
		--conf spark.executor.extraClassPath=./backend/hail-all-spark.jar \
		--conf spark.serializer=org.apache.spark.serializer.KryoSerializer \
		--conf spark.kryo.registrator=is.hail.kryo.HailKryoRegistrator \
		--conf spark.local.dir=$(pwd)/{params.tmp_dir} \
		--master spark://${{SLURMD_NODENAME}}.cluster.imbie:7077 \
		scripts/QC_HAIL_gather_data.py $(pwd)/{input} $(pwd)/{params.tmp_dir} $(pwd)/{output.sample_qc_file} $(pwd)
		
		else
    		echo "Running locally."
    		python scripts/QC_HAIL_gather_data.py $(pwd)/{input} $(pwd)/{params.tmp_dir} $(pwd)/{output.sample_qc_file} $(pwd)
		fi
		"""

# reformat the hail sample QC file, generate histograms of all values, generate a list of samples passing the sample QC thresholds 
rule QC_HAIL_convert_plot_data:
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
rule QC_HAIL_generate_plink_files_forGWAS:
	input:
		vcf="data/normalized/df3_norm.vcf.gz",
		table="data/hail_gather_data/for_sample_QC.tsv",
		hail_dir=directory(config["tmp_folder"] + "hail_gather_data_IN")
	output:
		vcf_for_GWAS="data/HAIL_GWAS_vcf/HAIL_GWAS_vcf.vcf.bgz"
	resources: cpus=4, mem_mb=10000, time_job=10080, additional=" -w b308eth0"
	params:
		partition='long',
		tmp_dir=config["tmp_folder"],
		path_output="data/HAIL_GWAS_vcf"
	shell:
		"""
		# set up spark cluster
		mkdir -p {params.path_output}
		mkdir -p {params.tmp_dir}
		export LC_ALL="en_US.UTF-8"
		
		if [ {config[cluster]} = "yes" ]; then
		
		echo "Setting up Spark Cluster."
		export SPARK_LOCAL_DIRS=$(pwd)/{params.tmp_dir}
		export PATH={config[PATH]}
		export PYTHONPATH={config[PYTHONPATH]}
		HAIL_HOME=$(pip3 show hail | grep Location | awk -F' ' '{{print $2 "/hail"}}')
		export SPARK_CLASSPATH=$HAIL_HOME"/backend/hail-all-spark.jar"		
		echo "Env variables set."
		# start spark master in this job, add one worker that will definitely start
		{config[spark_master_path]}
		{config[spark_worker_path]} "spark://${{SLURMD_NODENAME}}.cluster.imbie:7077" --cores 2
		
		# add addtitional workers, save job-ids and add a bash trap to cancel them, when this process stops
		slurm_id=""
		for i in b203eth0 b204eth0 b205eth0 b206eth0 b207eth0 b208eth0 b209eth0 b210eth0 b211eth0 b212eth0 b213eth0 b214eth0 b203eth0 b204eth0 b205eth0 b206eth0 b207eth0 b208eth0 b209eth0 b210eth0 b211eth0 b212eth0 b213eth0 b214eth0 b203eth0 b204eth0 b205eth0 b206eth0 b207eth0 b208eth0 b209eth0 b210eth0 b211eth0 b212eth0 b213eth0 b214eth0
		do
		   slurm_id="$slurm_id "$(sbatch --time=12:00:00 -J spark_worker --ntasks=1 --partition=batch \
		   --mem=28G --cpus-per-task=8 \
		   --chdir=$(pwd) \
		   -w $i \
		   --wrap "{config[spark_worker_path]} spark://${{SLURMD_NODENAME}}.cluster.imbie:7077; sleep 12h" | awk '{{print $NF}}')
		done
		
		echo "started workers as jobs:"
		echo $slurm_id
		trap "scancel $slurm_id" EXIT
		
		# submit the hail script
		
		spark-submit \
		--jars $HAIL_HOME/backend/hail-all-spark.jar \
		--conf spark.driver.extraClassPath=$HAIL_HOME/backend/hail-all-spark.jar \
		--conf spark.executor.extraClassPath=./backend/hail-all-spark.jar \
		--conf spark.serializer=org.apache.spark.serializer.KryoSerializer \
		--conf spark.kryo.registrator=is.hail.kryo.HailKryoRegistrator \
		--conf spark.local.dir=$(pwd)/{params.tmp_dir} \
		--master spark://${{SLURMD_NODENAME}}.cluster.imbie:7077 \
		scripts/QC_HAIL_generate_plink_files_forGWAS.py $(pwd)/{input.vcf} $(pwd)/{params.tmp_dir} $(pwd)/{input.table} $(pwd)/{input.hail_dir} $(pwd)/{output}
		
		else
    		echo "Running locally."
    		python scripts/QC_HAIL_generate_plink_files_forGWAS.py $(pwd)/{input.vcf} $(pwd)/{params.tmp_dir} $(pwd)/{input.table} $(pwd)/{input.hail_dir} $(pwd)/{output}
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
	conda:
		"envs/bcftools_plink_R.yaml"
	params:
		partition='batch',
		bed1="data/1000G/1000G_chr{contig}",
		bed2="data/1000G/1000G_chr{contig}_pruned"
	shell:
		"""
		if bcftools view -q 0.05:minor {input.vcf} | \
		bcftools norm -m-any --check-ref w -f {input.fasta} | \
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
	conda:
		"envs/bcftools_plink_R.yaml"
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
		
rule relatedness_filter:
	input:
		vcf_for_GWAS="data/HAIL_GWAS_vcf/HAIL_GWAS_vcf.vcf.bgz",
		populations="data/PCA/populations.txt",
		pheno_file="data/hail_gather_data/for_sample_QC.tsv"
	output:		
		PLINK_fam="data/pop_vcf/{population}_vcf_not_rel.vcf.bgz.fam",
		bed="data/pop_vcf/{population}_vcf_not_rel.vcf.bgz.bed",
		bim="data/pop_vcf/{population}_vcf_not_rel.vcf.bgz.bim"
	resources: cpus=12, mem_mb=90000, time_job=10080, additional=" -w b310eth0 --gres localtmp:400G --ntasks=1 "
	params:
		partition='long',
		out_dir="data/pop_vcf/",
		tmp_dir=config["tmp_folder"] +"{population}/",
		pop_vcf_rel="data/pop_vcf/{population}_vcf_rel.vcf.bgz",
		pop_vcf_not_rel="data/pop_vcf/{population}_vcf_not_rel.vcf.bgz",
	shell:
		"""
		# set up spark cluster
		mkdir -p {params.out_dir}
		mkdir -p {params.tmp_dir}
		export LC_ALL="en_US.UTF-8"
		
		if [ {config[cluster]} = "yes" ]; then
		export SPARK_LOCAL_DIRS=$SCRATCH_DIR
		
		echo "Setting up Spark Cluster."
		export PATH={config[PATH]}
		export PYTHONPATH={config[PYTHONPATH]}
		HAIL_HOME=$(pip3 show hail | grep Location | awk -F' ' '{{print $2 "/hail"}}')
		export SPARK_CLASSPATH=$HAIL_HOME"/backend/hail-all-spark.jar"		
		echo "Env variables set."
		# start spark master in this job, add one worker that will definitely start
		{config[spark_master_path]}
		#{config[spark_worker_path]} "spark://${{SLURMD_NODENAME}}.cluster.imbie:7077" --cores 2
		
		# add addtitional workers, save job-ids and add a bash trap to cancel them, when this process stops
		slurm_id=""
		for i in b311eth0 b312eth0 b313eth0 b314eth0
		do
		   slurm_id="$slurm_id "$(sbatch --time=7-00:00:00 -J spark_worker --ntasks=1 --partition=long \
		   --mem=180G --cpus-per-task=22 \
		   --chdir=$(pwd) \
		   --gres localtmp:400G \
		   -w $i \
		   --wrap "source $(pwd)/scripts/LOCAL_DIR.sh; {config[spark_worker_path]} spark://${{SLURMD_NODENAME}}.cluster.imbie:7077; sleep 12h" | awk '{{print $NF}}')
		done
		
		echo "started workers as jobs:"
		echo $slurm_id
		trap "scancel $slurm_id" EXIT
		
		# submit the hail script
		
		spark-submit \
		--jars $HAIL_HOME/backend/hail-all-spark.jar \
		--conf spark.driver.extraClassPath=$HAIL_HOME/backend/hail-all-spark.jar \
		--conf spark.executor.extraClassPath=./backend/hail-all-spark.jar \
		--conf spark.serializer=org.apache.spark.serializer.KryoSerializer \
		--conf spark.kryo.registrator=is.hail.kryo.HailKryoRegistrator \
		--conf spark.driver.memory=80G \
		--conf spark.executor.memory=170G \
		--master spark://${{SLURMD_NODENAME}}.cluster.imbie:7077 \
		scripts/relatedness_filter.py $(pwd)/{input.vcf_for_GWAS} $(pwd)/{params.tmp_dir} $(pwd)/{input.populations} $(pwd)/{input.pheno_file} $(pwd)/{params.pop_vcf_rel} $(pwd)/{params.pop_vcf_not_rel} $(pwd)/{params.out_dir} {wildcards.population}
		else
    		echo "Running locally."
		python scripts/relatedness_filter.py $(pwd)/{input.vcf_for_GWAS} $(pwd)/{params.tmp_dir} $(pwd)/{input.populations} $(pwd)/{input.pheno_file} $(pwd)/{params.pop_vcf_rel} $(pwd)/{params.pop_vcf_not_rel} $(pwd)/{params.out_dir} $(pwd) {wildcards.population}
		fi
		#tabix -p vcf {params.pop_vcf_rel}
		#tabix -p vcf {params.pop_vcf_not_rel}
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


rule additional_QC_prior_to_each_analysis:
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

rule generate_PCA_covar_for_each_analysis:
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


rule generate_covar_files_for_each_analysis:
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
 
rule GWAS_with_regenie:
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
		 --hwe 1E-15 midp \
		 --maf 0.01 \
		 --geno 0.1 \
		 --indep-pairwise 50 5 0.05 \
		 --out {params.pruned_plink}
		
		export LD_LIBRARY_PATH={config[LD_LIBRARY_PATH]}
		
		tools/regenie_v2.2.4.gz_x86_64_Linux_mkl \
		  --step 1 \
		  --bed ${{infile::-4}} \
		  --covarFile {input.pathCov} \
		  --phenoFile {input.pathPheno} \
		  --bt \
		  --lowmem \
		  --lowmem-prefix {params.tmp_dir} \
		  --extract {params.pruned_plink}.prune.in \
		  --bsize 1000 \
		  --out {params.regenie_step1}
		
		#step 2

		tools/regenie_v2.2.4.gz_x86_64_Linux_mkl \
		  --step 2 \
		  --minMAC 6 \
		  --covarFile {input.pathCov} \
		  --phenoFile {input.pathPheno} \
		  --bed ${{infile::-4}} \
		  --bt \
		  --htp DeCOI \
		  --firth --approx \
		  --firth-se \
		  --maxstep-null 5 \
		  --maxiter-null 10000 \
		  --pred {params.regenie_step1}_pred.list \
		  --bsize 200 \
		  --out {params.regenie_step2}
		  
		  rm {params.pruned_plink}* {params.regenie_step1}*
		  """
		  

rule generate_qq_plots:
	input:
		path_reg="data/regenie/GWAS_{population}_{pheno}.regenie"
	output:
		"data/regenie/GWAS_{population}_{pheno}.regenie_Pval_manhattan.png" #GWAS_EUR_A2.regenie_plot_Pval_manhattan.png
	resources: cpus=1, mem_mb=3000, time_job=720
	params:
		partition='batch',
		output_string="data/regenie/GWAS_{population}_{pheno}.regenie"
	shell:
		"""
		Rscript scripts/qqplot.R -f {input} -o {params.output_string} -c Chr -p Pval -b Pos
		"""



rule get_PRS_variants:
	input:
		mt_table=directory(config["tmp_folder"] + "hail_gather_data_IN"),
		prs_file="resources/hg38.PRS.tsv"
	output:		
		vcf_for_PRS="data/PRS/variants_for_prs.vcf.gz"
	resources: cpus=12, mem_mb=48000, time_job=720, additional=" -w g101eth0 --gres localtmp:200G --ntasks=1 "
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
		worker_nodes="g101eth0 g101eth0 g101eth0 g101eth0 g101eth0 g101eth0 g101eth0 g101eth0 g101eth0"
		source scripts/spark_submit_command.sh
		$spark_submit_command \
		$hail_python_script
		
		else
		python $hail_python_script
		fi
		"""


rule get_rohs:
	input:
		vcf_for_GWAS="data/HAIL_GWAS_vcf/HAIL_GWAS_vcf.vcf.bgz",
	output:		
		ROH_file="data/ROH/test_roh.tsv.gz",
		outfile2="data/ROH/another_file.txt" # you can add many more files like this
	resources: cpus=12, mem_mb=48000, time_job=720, additional=" -w g101eth0 --gres localtmp:200G --ntasks=1 "
	params:
		partition='batch',
		out_dir="data/ROH/",

		tmp_dir=config["tmp_folder"] +"PRS/",
	shell:
		"""
		mkdir -p {params.out_dir}
		#mkdir -p {params.tmp_dir}
		# local SSD Scratch would be available in: $SCRATCH_DIR
		
		# the code here could e.g. look like this:
		bcftools roh {input.vcf_for_GWAS} --AF-tag AF -G30 --threads 12 | \
		grep "RG" | grep 'chr[1-22]' | gzip > {output.ROH_file}
		
		# this would translate to
		bcftools roh data/HAIL_GWAS_vcf/HAIL_GWAS_vcf.vcf.bgz --AF-tag AF -G30 --threads 12 | \
		grep "RG" | grep 'chr[1-22]' | gzip > data/ROH/test_roh.tsv.gz
		
		# the second filename that is defined as output could be accessed like this:
		touch {output.outfile2}
		
		# Have fun! You can also add more rules.
		"""

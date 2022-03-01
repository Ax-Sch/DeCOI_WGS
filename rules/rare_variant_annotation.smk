
samples=""

rule get_bedfile:
	output:
		"data/bedfilesForPrefilter/ExonsClinVarSpliceAI.bed"
	params: 
		partition='batch',
		out_folder="data/bedfilesForPrefilter"
	resources: cpus=1, mem_mb=9000, time_job=720
	conda: "spark"
	shell:
		"""
		# software required: bedops, bedtools

		mkdir -p {params.out_folder}
		cd {params.out_folder}

		# get bedtools reference file
		wget -Nc https://raw.githubusercontent.com/arq5x/bedtools2/master/genomes/human.hg19.genome -O bedtools_hg19_ref_file.genome

		# GENCODE
		wget -Nc ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gff3.gz

		zcat gencode.v39.annotation.gff3.gz | \
		 awk '$3 == "exon"' - | \
		 sort -k1,1 -k4,4n -V | \
		 bedtools merge -i stdin | \
		 cut -f1-3 | gzip > gencode39Exons.bed.gz

		# add padding to gencode file
		bedtools slop -b 20 \
		-i gencode39Exons.bed.gz \
		-g bedtools_hg19_ref_file.genome | \
		gzip > gencode39ExonsPadded.bed.gz

		# ClinVar
		wget -Nc ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
		zcat clinvar.vcf.gz | grep "^#"  > clinvar_pathogenic.vcf
		zcat clinvar.vcf.gz | grep -v "^#" | egrep "pathogenic|Pathogenic" | grep -v "uncertain" >> clinvar_pathogenic.vcf
		bedtools merge -i clinvar_pathogenic.vcf | cut -f1-3 | gzip > clinvar_pathogenic.bed.gz

		# SpliceAI
		# filter splice AI for scores >0.5 and convert to BED
		zcat '/media/axel/Seagate Backup Plus Drive/annotation_sources/spliceai_scores.hg38.vcf.gz' | \
		egrep "0\.[5-9]|1\.0" | \
		awk 'BEGIN {{OFS="\t"}} {{ print $1, $2-1, $2 }}' | \
		gzip > \
		splice_ai_positions.bed.gz


		# merge and sort individual bed files
		zcat gencode39ExonsPadded.bed.gz clinvar_pathogenic.bed.gz splice_ai_positions.bed.gz | cut -f1-3 | bedtools sort | bedtools merge > "ExonsClinVarSpliceAI.bed"
		
		"""




rule split:
	input:
		vcf="data/rare_variant_vcf/rare_variant_vcf.vcf.bgz",
	output:
		"data/rare_chrom_split/all_samples_{contig}.vcf.gz"
	params:
		partition='batch'
	resources: cpus=1, mem_mb=4000, time_job=720, additional="",
	conda: "spark"
	shell:
		"""
		bcftools view {input} -r {wildcards.contig} -Oz -o {output}	
		tabix -pvcf {output}
		"""

rule vep:
	input:
		vcf="data/rare_chrom_split/all_samples_{contig}.vcf.gz"
	output:
		"data/rare_chrom_annotated/variants_{contig}.tsv.gz"
	params:
		partition='long',
		anno_sources=config["database_dir"],
		fasta_file=config["fasta"],
	conda:
		"covid_qc",
		#"envs/vep_101.yml"
	resources: cpus=4, mem_mb=40000, time_job=10040, additional=" --gres localtmp:160G -N1-1" #--exclude=b101eth0,b102eth0,b103eth0,b104eth0,b105eth0,b106eth0,b107eth0,b201eth0,b203eth0,b204eth0,b205eth0,b206eth0,b207eth0,b208eth0,b209eth0,b210eth0,b211eth0,b212eth0,b213eth0,b214eth0,b301eth0,b302eth0,b303eth0,b305eth0,b306eth0,b307eth0,b308eth0,b309eth0,b310eth0,b311eth0,b312eth0,b313eth0,b314eth0,g101eth0 -w " + annotation_node
	shell:
		"""
		if [ {config[cluster]} = "yes" ]; then
		# cluster content
		echo $TMPDIR
		fi
		
		database_dir="{params.anno_sources}"
		
		#store fasta to RAMdisk
		cp -v "{params.fasta_file}" /dev/shm/
		cp -v "{params.fasta_file}.fai" /dev/shm/
		cp -v "${{database_dir}}/gnomad.genomes.v3.1.1.sites.all_red.vcf.gz" /dev/shm/
		cp -v "${{database_dir}}/gnomad.genomes.v3.1.1.sites.all_red.vcf.gz.tbi" /dev/shm/
		
		fasta_name="$(basename -- '{params.fasta_file}')"
		
		
		vep -v \
		--cache \
		--offline \
		--fork 4 \
		--buffer 20000 \
		--hgvs \
		--everything \
		--force_overwrite \
		--offline \
		--pick_allele_gene \
		--tab \
		-a GRCh38 \
		--dir_cache "{params.anno_sources}" \
		--fasta /dev/shm/$fasta_name \
		--custom "${{database_dir}}/clinvar20220229.vcf.gz",ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
		--plugin dbNSFP,"${{database_dir}}/dbNSFP/dbNSFP4.0a_hg38.gz",Ensembl_transcriptid,Uniprot_acc,VEP_canonical,LRT_pred,SIFT_pred,MutationTaster_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,REVEL_score,CADD_raw,CADD_phred \
		--custom "/dev/shm/gnomad.genomes.v3.1.1.sites.all_red.vcf.gz",gnomAD_ge,vcf,exact,0,AF,AC,nhomalt \
		-i {input.vcf} \
		-o STDOUT | \
		sed "s/#Uploaded/Uploaded/g" | \
            	bgzip > {output}
		
		rm /dev/shm/${{fasta_name}}*
		rm /dev/shm/gnomad.genomes*
		
		#--custom "${{database_dir}}/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz",gnomAD_ex,vcf,exact,0,AF,AN,AC,AF_afr,AF_amr,AF_asj,AF_eas,AF_fin,AF_nfe,AF_oth,AF_sas,nhomalt \
		
		"""

rule get_freqs:
	input:
		vcf="data/rare_chrom_split/all_samples_{contig}.vcf.gz"
	output:
		"data/rare_chrom_annotated/variants_{contig}.frqx.gz"
	params:
		partition='long',
		plink_out=lambda wildcards, output: output[0][:-8],
	conda:
		"spark",
		#"envs/vep_101.yml"
	shell:
		"""
		plink \
		--vcf {input} \
		--freqx gz \
		--out {params.plink_out}
		"""

rule make_aaf_files:
	input:
		pathAnnot_tsv="data/rare_chrom_annotated/variants_{contig}.tsv.gz",
		pathAAF_tsv="data/rare_chrom_annotated/variants_{contig}.frqx.gz"
	output:
		regenie_anno_file="data/rare_chrom_categorized/{contig}.anno",
		regenie_aaf_file="data/rare_chrom_categorized/{contig}.aaf",
		regenie_set_list="data/rare_chrom_categorized/{contig}.set",
		variables_categorized="data/rare_chrom_categorized/{contig}.tsv.gz"
	script:
		"../scripts/make_aaf.R"




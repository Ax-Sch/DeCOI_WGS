
samples=""

rule get_bedfile:
	output:
		"data/bedfilesForPrefilter/ExonsClinVarSpliceAI.bed"
	params: 
		partition='batch',
		out_folder="data/bedfilesForPrefilter"
	resources: cpus=1, mem_mb=9000, time_job=720
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
	resources: cpus=1, mem_mb=4000, time_job=720, additional=""
	shell:
		"""
		bcftools view {input} -r chr{wildcards.contig} -Oz -o {output}	
		tabix -pvcf {output}
		"""

rule vep:
	input:
		vcf="data/rare_chrom_split/all_samples_{contig}.vcf.gz"
	output:
		"data/rare_chrom_annotated/variants_{contig}.vcf.gz"
	params:
		partition='long',
		anno_sources=config["database_dir"],
		fasta_file=config["fasta"],
	conda:
		"envs/vep_101.yml"
	resources: cpus=4, mem_mb=40000, time_job=10040, additional=" --gres localtmp:160G -N1-1" #--exclude=b101eth0,b102eth0,b103eth0,b104eth0,b105eth0,b106eth0,b107eth0,b201eth0,b203eth0,b204eth0,b205eth0,b206eth0,b207eth0,b208eth0,b209eth0,b210eth0,b211eth0,b212eth0,b213eth0,b214eth0,b301eth0,b302eth0,b303eth0,b305eth0,b306eth0,b307eth0,b308eth0,b309eth0,b310eth0,b311eth0,b312eth0,b313eth0,b314eth0,g101eth0 -w " + annotation_node
	shell:
		"""
		if [ {config[cluster]} = "yes" ]; then
		# cluster content
		echo $TMPDIR
		fi
		
		database_dir="{params.anno_sources}"
		
		#store fasta to RAMdisk
		cp -v "{params.fasta_file}"* /dev/shm/
		fasta_name="$(basename -- {params.fasta_file})
		
		#wget -N ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz  
		#wget -N ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
		
		cd $work_dir
		
		vep -v \
		--cache \
		--offline \
		--force_overwrite \
		--fork 4 \
		--buffer 20000 \
		--hgvs \
		--no_stats \
		--format vcf \
		--numbers \
		--vcf \
		--gene_phenotype \
		--pick_allele_gene \
		--af \
		--af_esp \
		--check_existing \
		-a GRCh38 \
		--dir_cache $TMPDIR \
		--fasta /dev/shm/$fasta_name \
		--custom "${{database_dir}}/clinvar20220229.vcf.gz",ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
		--plugin dbNSFP,"${{database_dir}}/dbNSFP/dbNSFP4.0a_hg38.gz",Ensembl_proteinid,MutationTaster_pred,MutationTaster_score,PROVEAN_pred,PROVEAN_score,PrimateAI_pred,PrimateAI_score,Polyphen2_HVAR_pred,SIFT_pred,UK10K_AF,REVEL_score,REVEL_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,BayesDel_noAF_score,BayesDel_noAF_pred,gnomAD_genomes_AC,gnomAD_genomes_nhomalt \
		-i {input.vcf} \
		-o -o STDOUT | \
            	cut -f1-8 | \
            	bgzip > {output}
		
		#--custom $database_dir'/gnomad.genomes.v3.1.1.sites.all_red.vcf.gz',gnomAD_ge,vcf,exact,0,AF,AC,nhomalt \
		#--custom $database_dir'/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz',gnomAD_ex,vcf,exact,0,AF,AN,AC,AF_afr,AF_amr,AF_asj,AF_eas,AF_fin,AF_nfe,AF_oth,AF_sas,nhomalt \
		
		rm /dev/shm/$fasta_name
		
		"""

rule join:
	input:
		in_snv="output/all_samples_genome_annotated_snv_{contig}.vcf",
		in_cnv="output/all_samples_genome_annotated_cnv_{contig}.vcf"
	output:
		"output/annotated_{contig}.vcf.gz"
	params:
		partition='batch'
	resources: cpus=1, mem_mb=4000, time_job=720, additional=""
	shell:
		"""
		cat {input.in_snv} | bgzip > {input.in_snv}.bgz 
		cat {input.in_cnv} | bgzip > {input.in_cnv}.bgz
		tabix -pvcf {input.in_snv}.bgz 
		tabix -pvcf {input.in_cnv}.bgz 
		
		bcftools concat -Oz -a {input.in_snv}.bgz {input.in_cnv}.bgz > {output}
		"""










rule filtering:
	input:
		"output/annotated_{contig}.vcf.gz"
	output:
		"output/annotated_split_vep_{contig}.tsv"
	params:
		partition='batch',
		vep_filtered="output/annotated_filtered_{contig}.vcf",
		vep_renamed="output/annotated_filtered_renamed_{contig}.vcf"
	conda:
		"envs/vep_101.yml"
	resources: cpus=1, mem_mb=9000, time_job=720
	shell:
		"""
		function return_AF_string_w_genomes
		{{
		echo "(AF_1000G < "$1" or not AF_1000G) and (gnomAD_ex_AF_afr < "$1" or not gnomAD_ex_AF_afr) and (gnomAD_ex_AF_amr < "$1" or not gnomAD_ex_AF_amr) and (gnomAD_ex_AF_eas < "$1" or not gnomAD_ex_AF_eas) and (gnomAD_ex_AF_fin < "$1" or not gnomAD_ex_AF_fin) and (gnomAD_ex_AF_nfe < "$1" or not gnomAD_ex_AF_nfe) and (gnomAD_ex_AF_oth < "$1" or not gnomAD_ex_AF_oth) and (gnomAD_ex_AF_sas < "$1" or not gnomAD_ex_AF_sas) and (gnomAD_ge_AF_AFR < "$1" or not gnomAD_ge_AF_AFR) and (gnomAD_ge_AF_AMR < "$1" or not gnomAD_ge_AF_AMR) and (gnomAD_ge_AF_EAS < "$1" or not gnomAD_ge_AF_EAS) and (gnomAD_ge_AF_FIN < "$1" or not gnomAD_ge_AF_FIN) and (gnomAD_ge_AF_NFE < "$1" or not gnomAD_ge_AF_NFE) and (gnomAD_ge_AF_SAS < "$1" or not gnomAD_ge_AF_SAS) and (AA_AF < "$1" or not AA_AF) and (EA_AF < "$1" or not EA_AF)"
		}}
		filter_vep_str="$(return_AF_string_w_genomes 0.02)"
		
		bcftools_exec="/gpfs01/homedirs/aschmidt/annotation_sources/bcftools/bcftools"
		export BCFTOOLS_PLUGINS='/gpfs01/homedirs/aschmidt/annotation_sources/bcftools/plugins'
		
		echo "filter_vep -i  {input} --format vcf -o {params.vep_filtered} --force_overwrite --filter \\""$filter_vep_str"\\"" | bash
		cat {params.vep_filtered} | sed 's:;AN=:;AN_cohort=:g' | sed 's:ID=AN,:ID=AN_cohort,:g' | sed 's:;AC=:;AC_cohort=:g' | sed 's:ID=AC,:ID=AC_cohort,:g' | sed 's:||:|0|:g'| sed 's:||:|0|:g' > {params.vep_renamed}
		
		split_string1='%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%Allele\\t%QUAL\\t%AC_cohort\\t%AN_cohort\\t%IMPACT\\t%Consequence\\t%HGVSc\\t%HGVSp\\t%SYMBOL\\t%Existing_variation\\t%QD\\t%FS\\t%MQ\\t%MQRankSum\\t%ReadPosRankSum\\t'
		split_string2='%AA_AF\\t%EA_AF\\t%CLIN_SIG\\t%SOMATIC\\t%PHENO\\t%BayesDel_noAF_pred\\t%BayesDel_noAF_score\\t%Ensembl_proteinid\\t%MutationTaster_pred\\t%MutationTaster_score\\t%PROVEAN_pred\\t%PROVEAN_score\\t%Polyphen2_HVAR_pred\\t%PrimateAI_pred\\t%PrimateAI_score\\t%REVEL_rankscore\\t%REVEL_score\\t%SIFT_pred\\t%UK10K_AF\\t%gnomAD_genomes_AC\\t%gnomAD_genomes_nhomalt\\t%phastCons100way_vertebrate\\t%phastCons100way_vertebrate_rankscore\\t%phyloP100way_vertebrate\\t%phyloP100way_vertebrate_rankscore\\t%CADD_PHRED\\t%CADD_RAW\\t%SV_overlap_AF\\t%SV_overlap_PC\\t%SV_overlap_name\\t%ClinVar\\t%ClinVar_CLNSIG\\t%ClinVar_CLNREVSTAT\\t%ClinVar_CLNDN\\t%gnomAD_ge\\t%gnomAD_ge_AF\\t%gnomAD_ge_AC\\t%gnomAD_ge_nhomalt\\t%gnomAD_ex\\t%gnomAD_ex_AF\\t%gnomAD_ex_AN\\t%gnomAD_ex_AC\\t%gnomAD_ex_AF_afr\\t%gnomAD_ex_AF_amr\\t%gnomAD_ex_AF_asj\\t%gnomAD_ex_AF_eas\\t%gnomAD_ex_AF_fin\\t%gnomAD_ex_AF_nfe\\t%gnomAD_ex_AF_oth\\t%gnomAD_ex_AF_sas\\t%gnomAD_ex_nhomalt\\t%SpliceAI\\t%SpliceAI_ALLELE\\t%SpliceAI_SYMBOL\\t%SpliceAI_DS_AG\\t%SpliceAI_DS_AL\\t%SpliceAI_DS_DG\\t%SpliceAI_DS_DL\\t%SVLEN\\t%SVTYPE\\t%END\\t'
		split_string3='[%GT\\t%TGT\\t%AD\\t%GQ\\t]' # 

		header1=$(echo $split_string1$split_string2 | sed 's:\\\\t:        :g')
		header2=$(bcftools query -u -H -f $split_string3'\\n' {params.vep_renamed} | head -n1 | cut -f2- -d " ") ||  echo ""
		echo "sample        ${{header1}}        $header2" | awk -v OFS="\\t" '$1=$1' | cut -f1-50 > {output}
		for sample in `bcftools query -l {params.vep_renamed}`; do
			bcftools view -c1 -s $sample {params.vep_renamed}  | \
			$bcftools_exec +split-vep -f $split_string1$split_string2$split_string3'\\n' -d -A tab | \
			awk -v var="$sample"  -F $'\t' 'BEGIN {{OFS = FS}} {{print var, $0 }}' >> {output}
		done
		rm -f {params.vep_filtered} {params.vep_renamed}
		"""


rule extract_tsv:
	input:
		"output/{sample}/{sample}.exome.vcf.gz.annotated.zip"
	output:
		"output/{sample}/{sample}_prejoin.tsv"
	params:
		partition='batch',
		folder="output/{sample}/",
		tsv_file="output/{sample}/output/{sample}/{sample}.exome.vcf.gz.VASE_FILTERED_TSV.anno_tab"
	conda:
		"envs/vep_101.yml"
	log:
		log="logs/{sample}/{sample}.extract_tsv.log"
	resources: cpus=1, mem_mb=9000, time_job=720
	shell:
		"""
		unzip -o {input} -d {params.folder}
		tail -n+3 {params.tsv_file} | awk '{{print "{wildcards.sample}\t",$0}}' > {output}
		"""


rule combine_tsv:
	input:
		expand("output/{sample}/{sample}_prejoin.tsv",sample=samples)
	output:
		"all_joinded.tsv"
	params:
		partition='batch',
	resources: cpus=1, mem_mb=9000, time_job=720
	shell:
		"""
		cat {input} > {output}
                """


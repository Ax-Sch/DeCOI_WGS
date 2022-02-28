#!/bin/bash
#SBATCH -J subsetVCF   
#SBATCH --time=07-00:00:00 
#SBATCH --partition=long    #long
#SBATCH --threads-per-core=1
#SBATCH --mem=4G
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --verbose
#SBATCH --begin=now
#SBATCH --chdir=/ceph01/scratch/aschmidt/WGS_GWAS_pipeline 

module load bcftools

mkdir -p data/subset_vcf

#bcftools view -S resources/samples.txt --min-ac=1 /ceph01/projects/covid_wgs/glnexus/output/covid_wgs.bcf -Oz > data/subset_vcf/subseted.vcf.gz

#tabix -p vcf data/subset_vcf/subseted.vcf.gz

#md5sum data/subset_vcf/subseted.vcf.gz > data/subset_vcf/md5sum_vcf.txt


bcftools view -S resources/samples.txt --force-samples --min-ac=1 --force-samples "data/HAIL_GWAS_vcf/HAIL_GWAS_vcf.vcf.bgz" -Oz > data/subset_vcf/subseted_Vcf_for_GWAS_after_first_QC.vcf.gz

tabix -p vcf data/subset_vcf/subseted_Vcf_for_GWAS_after_first_QC.vcf.gz

md5sum data/subset_vcf/subseted_Vcf_for_GWAS_after_first_QC.vcf.gz -Oz > data/subset_vcf/subseted_Vcf_for_GWAS_after_first_QC_md5sum.txt

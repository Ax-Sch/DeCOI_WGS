#!/bin/bash
#SBATCH -J  ala          # Job name
#SBATCH --time=07-00:00:00                  # Walltime
#SBATCH --partition=long    #long
#SBATCH --threads-per-core=1
#SBATCH --mem=1G
#SBATCH --ntasks=1                       # 1 tasks
             # number of cores per task                       # number of nodes
#SBATCH --cpus-per-task=1
#SBATCH --verbose
#SBATCH --begin=now
#SBATCH --chdir=/ceph01/scratch/aschmidt/GWAS_WGS

module load anaconda/4.6.11-py37

conda activate WGS_GWAS
#snakemake --unlock -s /ceph01/homedirs/aschmidt/WGS_GWAS_pipeline/Snakefile

snakemake -s /ceph01/homedirs/aschmidt/WGS_GWAS_pipeline/Snakefile --profile slurm --rerun-incomplete # --cleanup-metadata hail.full.normID.noChrM.mt

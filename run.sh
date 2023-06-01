#!/bin/bash
#SBATCH -J  COVID_WGS          # Job name
#SBATCH --time=02-00:00:00                  # Walltime
#SBATCH --partition=medium    #long
#SBATCH --threads-per-core=1
#SBATCH --mem=5G
#SBATCH --ntasks=1                       # 1 tasks
             # number of cores per task                       # number of nodes
#SBATCH --cpus-per-task=1
#SBATCH --verbose
#SBATCH --begin=now
#SBATCH --chdir=/home/aschmidt/COVID_WGS
#SBATCH --account=ag_ihg_ludwig

#conda activate snakemake6
snakemake --profile slurm --rerun-incomplete --use-conda --latency-wait 60  # --cleanup-metadata hail.full.normID.noChrM.mt
#snakemake --profile slurm --rerun-incomplete --use-conda --latency-wait 30 dev_test_rvtest # --cleanup-metadata hail.full.normID.noChrM.mt

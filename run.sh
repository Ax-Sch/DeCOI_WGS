#!/bin/bash
#SBATCH -J  COVID_WGS         			# Job name, ADJUST as you like
#SBATCH --time=02-00:00:00    			# Probably ADJUST THIS to the constraints of your slurm setup
#SBATCH --partition=medium    			#long, probably ADJUST THIS to the constraints of your slurm setup
#SBATCH --threads-per-core=1
#SBATCH --mem=5G
#SBATCH --ntasks=1				# 1 tasks
#SBATCH --cpus-per-task=1
#SBATCH --verbose
#SBATCH --begin=now
#SBATCH --chdir=/home/aschmidt/COVID_WGS    	#ADJUST THIS to the directory this file is located in
#SBATCH --account=ag_ihg_ludwig             	#ADJUST THIS or remove it if not needed

# depending on your cluster you might need the line "conda activate snakemake7" here.
snakemake --profile slurm --rerun-incomplete --use-conda --latency-wait 60  

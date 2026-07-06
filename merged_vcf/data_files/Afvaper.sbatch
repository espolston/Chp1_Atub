#!/bin/bash
#SBATCH --job-name=afvaper
#SBATCH --output=afvaper.out
#SBATCH --error=afvaper.err
#SBATCH --time=2-00:00:00
#SBATCH --partition=tier3q
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=0

module load gcc/12.1.0
module load R/4.5.1

Rscript /scratch/espolston/Afvaper_allchr.R

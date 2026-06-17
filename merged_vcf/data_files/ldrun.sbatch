#!/bin/bash
#SBATCH --job-name=LD
#SBATCH --output=LD.out
#SBATCH --error=LD.err
#SBATCH --time=00:30:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=20GB
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=espolston@rcc.uchicago.edu

module load plink
plink --bfile /scratch/midway2/espolston/binary_dcgm --fam /scratch/midway2/espolston/dcgm.fam --thin 0.1  --r2 gz --ld-window-kb 200 --ld-window-r2 0 --allow-extra-chr --allow-no-sex --make-founders --out dcgm_ld


#!/bin/bash
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate shapeit_env

for sc in {1..16}; do 
#make regions files from 10kb win file
awk -v scaf=${sc} 'BEGIN{OFS="\t"} $1 == "Scaffold_" scaf {$5 = $2 + 5000; print $1,$5}' dcgm_sig.bed > dcgm_sig_${sc}.txt
bcftools view --regions-file dcgm_sig_${sc}.txt /scratch/midway3/espolston/new_phased/phased_dcgm_Scaffold_${sc}_newhead.vcf.gz > /scratch/midway3/espolston/afvaper/phased_dcgm_Scaffold_${sc}_newhead_cmhsite.vcf.gz
tabix -p vcf /scratch/midway3/espolston/new_phased/phased_dcgm_Scaffold_${sc}_newhead_cmhsite.vcf.gz; done

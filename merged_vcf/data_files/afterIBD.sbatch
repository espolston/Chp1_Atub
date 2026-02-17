#!/bin/bash
#SBATCH --job-name=afterIBD
#SBATCH --output=afterIB.out
#SBATCH --error=afterIB.err
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=espolston@rcc.uchicago.edu 

module load plink
module load htslib
module load python/anaconda-2023.09
source activate /project/kreiner/espolston/

#remove related indivs from vcfs
bcftools view --samples ^4875,4828,4933,5171,5203,5244,5247,5381,556,P16_Ag_17_T,P17_Ag_4_T,P17_Nat_19_T,P4_Nat_16_T,P17_Nat_7_T,5066,5202,5209,4830,P4_Ag_26_T,4979,P7_Nat_1_T,5084,5220,5258,P4_Ag_1_T,P17_Ag_7_T,P4_Nat_3_T  -Oz -o /scratch/midway3/espolston/vcf/drought_commongarden_merged_final_ID_IBD.vcf.gz /scratch/midway3/espolston/vcf/drought_commongarden_merged_finalprocessed_ID.vcf.gz

bcftools view --samples ^4875,4828,4933,5171,5203,5244,5247,5381,556,P16_Ag_17_T,P17_Ag_4_T,P17_Nat_19_T,P4_Nat_16_T,P17_Nat_7_T,5066,5202,5209,4830,P4_Ag_26_T,4979,P7_Nat_1_T,5084,5220,5258,P4_Ag_1_T,P17_Ag_7_T,P4_Nat_3_T  -Oz -o /scratch/midway3/espolston/vcf/droughtcgmerged_allsites_final_ID_IBD.vcf.gz /scratch/midway3/espolston/vcf/droughtcgmerged_allsites_finalprocessed_ID.vcf.gz

#index and copy to project
tabix -p vcf /scratch/midway3/espolston/vcf/drought_commongarden_merged_final_ID_IBD.vcf.gz
tabix -p vcf /scratch/midway3/espolston/vcf/droughtcgmerged_allsites_final_ID_IBD.vcf.gz 

cp /scratch/midway3/espolston/vcf/drought_commongarden_merged_final_ID_IBD.vcf.gz /project/kreiner/pairedenv_commongarden/normalized_SNPsonly_vcfs/
cp /scratch/midway3/espolston/vcf/drought_commongarden_merged_final_ID_IBD.vcf.gz.tbi /project/kreiner/pairedenv_commongarden/normalized_SNPsonly_vcfs/
cp /scratch/midway3/espolston/vcf/droughtcgmerged_allsites_final_ID_IBD.vcf.gz /project/kreiner/pairedenv_commongarden/allsites/
cp /scratch/midway3/espolston/vcf/droughtcgmerged_allsites_final_ID_IBD.vcf.gz.tbi /project/kreiner/pairedenv_commongarden/allsites/

#make binary files (overwrite ones that were made before)
plink --vcf /project/kreiner/pairedenv_commongarden/normalized_SNPsonly_vcfs/drought_commongarden_merged_final_ID_IBD.vcf.gz --allow-extra-chr --make-bed --double-id --out /scratch/midway3/espolston/vcf/binary_dcgm

perl /scratch/midway3/espolston/vcf/fam_edit_merged.pl /scratch/midway3/espolston/vcf/binary_dcgm.fam /scratch/midway3/espolston/vcf/Chp1_DroughtandCommonGardenInfo_MergedMetaData.tsv /scratch/midway3/espolston/vcf/dcgm.fam

plink --bfile /scratch/midway3/espolston/vcf/binary_dcgm --bim /scratch/midway3/espolston/vcf/binary_dcgm.bim --fam /scratch/midway3/espolston/vcf/dcgm.fam --mh --within /scratch/midway3/espolston/vcf/dcgm.fam --allow-extra-chr --allow-no-sex --nonfounders --out /scratch/midway3/espolston/vcf/dcgm_cmh

#run pca again after IBD rem
plink --bfile /scratch/midway3/espolston/vcf/binary_dcgm --bim /scratch/midway3/espolston/vcf/binary_dcgm.bim --fam /scratch/midway3/espolston/vcf/dcgm.fam --pca --allow-no-sex --allow-extra-chr --make-founders --out /scratch/midway3/espolston/vcf/dcgm_pca_after_IBD

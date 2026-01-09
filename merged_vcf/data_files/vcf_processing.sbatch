#!/bin/bash
#SBATCH --job-name=processing
#SBATCH --output=processing.out
#SBATCH --error=processing.err
#SBATCH --time=12:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=50GB
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=espolston@rcc.uchicago.edu

module load htslib
module load plink
module load python/anaconda-2023.09
source activate /project/kreiner/espolston/

mkdir -p /scratch/midway3/espolston/vcf/

#VARIANT
#.05 missingness and MAF .01 filters (for variant file), remove 2 indivs (P16_Nat_1_T and P12_Nat_14_T), rename all drought P12 to P10 and check for duplicates
#add ID to vcf as well

#had to reheader the vcf file 
#bcftools view -h /project/kreiner/pairedenv_commongarden/normalized_SNPsonly_vcfs/drought_commongarden_merged_normalized_filtsnps_nomiss30.vcf.gz > /scratch/midway3/espolston/header.txt
#vim /scratch/midway3/espolston/header.txt
#i > edit FORMAT=<ID=GL,Number=. to Number=G > esc > :wq
#bcftools reheader -h /scratch/midway3/espolston/header.txt -o /project/kreiner/pairedenv_commongarden/normalized_SNPsonly_vcfs/drought_commongarden_merged_normalized_filtsnps_nomiss30_rehead.vcf.gz /project/kreiner/pairedenv_commongarden/normalized_SNPsonly_vcfs/drought_commongarden_merged_normalized_filtsnps_nomiss30.vcf.gz
#tabix -p vcf /project/kreiner/pairedenv_commongarden/normalized_SNPsonly_vcfs/drought_commongarden_merged_normalized_filtsnps_nomiss30_rehead.vcf.gz

bcftools view --samples ^P16_Nat_1_T,P12_Nat_14_T \
/project/kreiner/pairedenv_commongarden/normalized_SNPsonly_vcfs/drought_commongarden_merged_normalized_filtsnps_nomiss30_rehead.vcf.gz | \
bcftools filter -i 'F_MISSING<=0.05' \
  -Oz -o /scratch/midway3/espolston/vcf/drought_commongarden_merged_normalized_filtsnps_miss05.vcf.gz

bcftools filter -i 'AF >= 0.01 & AF <= 0.99' /scratch/midway3/espolston/vcf/drought_commongarden_merged_normalized_filtsnps_miss05.vcf.gz -Oz -o /scratch/midway3/espolston/vcf/drought_commongarden_merged_normalized_filtsnps_miss05_MAF.vcf.gz

#reheader vcf for all P12 to be P10
#made rehead_mergedvcf.txt from copying header order and then editing P12 to P10
bcftools reheader --samples rehead_mergedvcf.txt -o /scratch/midway3/espolston/vcf/drought_commongarden_merged_rename.vcf.gz -Oz /scratch/midway3/espolston/vcf/drought_commongarden_merged_normalized_filtsnps_miss05_MAF.vcf.gz

awk '{OFS="\t"} /^#/ {print; next} {split($1, array, "_"); $3 = array[2] ":" $2; print}' /scratch/midway3/espolston/vcf/drought_commongarden_merged_rename.vcf.gz > /scratch/midway3/espolston/vcf/drought_commongarden_merged_finalprocessed_ID.vcf.gz

bgzip -c /scratch/midway3/espolston/vcf/drought_commongarden_merged_finalprocessed_ID.vcf > /scratch/midway3/espolston/vcf/drought_commongarden_merged_finalprocessed_ID.vcf.gz

tabix -p vcf /scratch/midway3/espolston/vcf/drought_commongarden_merged_finalprocessed_ID.vcf.gz

#test integrity of file
bgzip -t drought_commongarden_merged_finalprocessed_ID.vcf.gz

plink --vcf /scratch/midway3/espolston/vcf/drought_commongarden_merged_finalprocessed_ID.vcf.gz --allow-extra-chr --make-bed --double-id --out /scratch/midway3/espolston/vcf/binary_dcgm

#check for duplicates
plink --vcf /scratch/midway3/espolston/vcf/drought_commongarden_merged_finalprocessed_ID.vcf.gz --allow-extra-chr --double-id --list-duplicate-vars --out duplicate_vars

#checked duplicate_vars.dupvar and none found

cp /scratch/midway3/espolston/vcf/drought_commongarden_merged_finalprocessed_ID.vcf.gz /project/kreiner/pairedenv_commongarden/normalized_SNPsonly_vcfs/

cp /scratch/midway3/espolston/vcf/drought_commongarden_merged_finalprocessed_ID.vcf.gz.tbi /project/kreiner/pairedenv_commongarden/normalized_SNPsonly_vcfs/

#run mds and pcs
#first edit fam file
perl /scratch/midway3/espolston/vcf/fam_edit_merged.pl /scratch/midway3/espolston/vcf/binary_dcgm.fam /scratch/midway3/espolston/vcf/Chp1_DroughtandCommonGardenInfo_MergedMetaData.tsv /scratch/midway3/espolston/vcf/dcgm.fam

plink --bfile /scratch/midway3/espolston/vcf/binary_dcgm --bim /scratch/midway3/espolston/vcf/binary_dcgm.bim --fam /scratch/midway3/espolston/vcf/dcgm.fam --pca --allow-no-sex --allow-extra-chr --make-founders --out /scratch/midway3/espolston/vcf/dcgm_pca

plink --bfile /scratch/midway3/espolston/vcf/binary_dcgm --bim /scratch/midway3/espolston/vcf/binary_dcgm.bim --fam /scratch/midway3/espolston/vcf/dcgm.fam --mds-plot 20 --within /scratch/midway3/espolston/vcf/dcgm.fam --cluster --allow-no-sex --allow-extra-chr --make-founders --out /scratch/midway3/espolston/vcf/dcgm_mds

plink --bfile /scratch/midway3/espolston/vcf/binary_dcgm --bim /scratch/midway3/espolston/vcf/binary_dcgm.bim --fam /scratch/midway2/espolston/vcf/dcgm.fam --mh --within /scratch/midway2/espolston/dcgm.fam --allow-extra-chr --allow-no-sex --nonfounders --out /scratch/midway3/espolston/vcf/dcgm_cmh

#INVARIANT
#.05 missingness, remove 2 indivs (P16_Nat_1_T and P10_Nat_14_T), rename all drought P12 to P10, add ID

bcftools view --samples ^P16_Nat_1_T,P12_Nat_14_T /project/kreiner/pairedenv_commongarden/allsites/droughtcgmerged_allsites_filtered.vcf.gz | bcftools filter -i 'F_MISSING<=0.05' -Oz -o /scratch/midway3/espolston/vcf/droughtcgmerged_allsites_filtered_.05miss.vcf.gz

bcftools reheader --samples rehead_mergedvcf.txt -Oz -o /scratch/midway3/espolston/vcf/droughtcgmerged_allsites_rename.vcf.gz /scratch/midway3/espolston/vcf/droughtcgmerged_allsites_filtered_.05miss.vcf.gz

bgzip -d /scratch/midway3/espolston/vcf/droughtcgmerged_allsites_rename.vcf.gz

awk '{OFS="\t"} /^#/ {print; next} {split($1, array, "_"); $3 = array[2] ":" $2; print}' /scratch/midway3/espolston/vcf/droughtcgmerged_allsites_rename.vcf > /scratch/midway3/espolston/vcf/droughtcgmerged_allsites_finalprocessed_ID.vcf

bgzip -c /scratch/midway3/espolston/vcf/droughtcgmerged_allsites_finalprocessed_ID.vcf > /scratch/midway3/espolston/vcf/droughtcgmerged_allsites_finalprocessed_ID.vcf.gz

tabix -p vcf /scratch/midway3/espolston/vcf/droughtcgmerged_allsites_finalprocessed_ID.vcf.gz

#test integrity of file
bgzip -t droughtcgmerged_allsites_finalprocessed_ID.vcf.gz 

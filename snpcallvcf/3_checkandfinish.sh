#!/bin/bash

#Check all regions have completed filtering
REGIONS_FILE="/scratch/midway3/espolston/region_list/scaf_1.txt"
DIR="/scratch/midway3/espolston/scaf1/"
DONE_FILE="/scratch/midway3/espolston/scaf_1_check.txt"

> "${DONE_FILE}"

while read line; do
  prefix=$(echo "$line" | sed 's/\.vcf.gz$//')
  echo "${DIR}/${prefix}_snps.vcf.gz.tbi" >> "${DONE_FILE}"
  echo "${DIR}/${prefix}_invariant.vcf.gz.tbi" >> "${DONE_FILE}"
done < "$REGIONS_FILE"

all_files_exist=true

while read file; do
  if [ ! -f "$file" ]; then
    echo "Error: File not found: $file"
    all_files_exist=false
  fi
done < "$DONE_FILE"


if "$all_files_exist"; then
  echo "All files in the list exist."
else
  echo "Some files in the list are missing."
fi


#Merge the regions into 1 variant and 1 invariant vcf per scaffold
#double check the options- allow overlaps etc
bcftools concat FILE LIST -O z -o Scaffold_1_invariant_oldname.vcf.gz

#reheader these files -- if from freebayes
#bcftools reheader --samples sample_rename.txt -O z -o Scaffold_1_snps.vcf.gz Scaffold_1_snps_oldname.vcf.gz
#bcftools reheader --samples sample_rename.txt -O z -o Scaffold_1_invariant.vcf.gz Scaffold_1_invariant_oldname.vcf.gz

#Edit vcf ID column for clumping down the line
#STILL NEED TO EDIT THIS LINE
awk '{OFS="\t"} /^#/ {print; next} {split($1, array, "_"); $3 = array[2] ":" $2; print}' /scratch/midway3/espolston/fixed_commongarden_allfiltsnps_193_hap2.vcf > /scratch/midway2/espolston/fixed_commongarden_allfiltsnps_193_hap2_ID.vcf

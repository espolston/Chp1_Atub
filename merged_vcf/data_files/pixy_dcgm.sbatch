#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=5
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --job-name=pixy
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=espolston@rcc.uchicago.edu 

#pixy: https://pixy.readthedocs.io/en/latest/index.html

#make population file
#get sample names from vcf
#module load python/anaconda-2023.09
#source activate /project/kreiner/espolston/
#bcftools query -l /project/kreiner/pairedenv_commongarden/allsites/droughtcgmerged_allsites_final_ID_IBD_rehead.vcf.gz  > /scratch/midway3/espolston/sample_ids.txt

#add ag/nat as second column
#PIXY NEEDS COMPARISONS AS STRINGS NOT 1/2!!!!
#perl pixy_pop.pl /scratch/midway3/espolston/sample_ids.txt /scratch/midway3/espolston/vcf/dcgm.fam /scratch/midway3/espolston/dcgm_namepopfile.txt /scratch/midway3/espolston/dcgm_pair

#make sites files
#bed file of dcgm clump sig sites - done in chp1_mergedvcf.Rmd
#/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/dcgm_sig.bsd

module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate pixy

module load parallel

#did run this one
#pixy on all pairs for clumped sites
#parallel -j 16 'pixy --stats pi fst dxy \
#  --vcf /project/kreiner/pairedenv_commongarden/allsites/droughtcgmerged_allsites_final_ID_IBD.vcf.gz \
#  --populations /scratch/midway3/espolston/dcgm_namepopfile.txt \
#  --bed_file /scratch/midway3/espolston/dcgm_sig.bed \
#  --chromosomes 'Scaffold_{1}' \
#  --output_folder /scratch/midway3/espolston/agnat_{1} \
#  --n_core 5' ::: {1..16}


#pixy on all sites
pixy --stats pi fst \
  --vcf /project/kreiner/pairedenv_commongarden/allsites/droughtcgmerged_allsites_final_ID_IBD.vcf.gz \
  --populations /scratch/midway3/espolston/dcgm_namepopfile.txt \
  --window_size 10000 \
  --output_folder /scratch/midway2/espolston/pixy_formanplots10 \
  --n_core 5

#AG VS NAT for pairs for clumped sites
parallel -j 17 'pixy --stats pi fst dxy \
  --vcf /project/kreiner/pairedenv_commongarden/allsites/droughtcgmerged_allsites_final_ID_IBD.vcf.gz \
  --populations /scratch/midway3/espolston/dcgm_pair_{1}.txt \
  --bed_file /scratch/midway3/espolston/dcgm_sig.bed \
  --chromosomes Scaffold_{2} \
  --output_folder /scratch/midway3/espolston/pair_{1}_scaf_{2} \
  --n_core 5' ::: {1..17} ::: {1..16}

#AG VS NAT for pairs for 10 kb windows
parallel -j 17 'pixy --stats fst \
  --vcf /project/kreiner/pairedenv_commongarden/allsites/droughtcgmerged_allsites_final_ID_IBD.vcf.gz \
  --populations /scratch/midway3/espolston/dcgm_pair_{1}.txt \
  --window_size 10000 \
  --output_folder /scratch/midway3/espolston/pair_{1}_windows \
  --n_core 5' ::: {1..17}

#Ag vs Nat for 50 mb windows
parallel -j 17 'pixy --stats fst \
  --vcf /project/kreiner/pairedenv_commongarden/allsites/droughtcgmerged_allsites_final_ID_IBD.vcf.gz \
  --populations /scratch/midway3/espolston/dcgm_pair_{1}.txt \
  --window_size 50000000 \
  --output_folder /scratch/midway3/espolston/pair_{1}_50Mb_windows \
  --n_core 5' ::: {1..17}

#pair ag and nat compared to all other pair ag nat for 50 mb windows
for p in {1..17}; do awk -v pa="$p" '{$3 = "pair_" pa "_" $2; print $1,$3}' dcgm_pair_$p.txt; done > dcgm_all.txt

#had to run this with pixy2
for CHROM in Scaffold_1 Scaffold_2 Scaffold_3 Scaffold_4 Scaffold_5 \
             Scaffold_6 Scaffold_7 Scaffold_8 Scaffold_9 Scaffold_10 \
             Scaffold_11 Scaffold_12 Scaffold_13 Scaffold_14 Scaffold_15 Scaffold_16; do

pixy --stats fst \
  --vcf /project/kreiner/pairedenv_commongarden/allsites/droughtcgmerged_allsites_final_ID_IBD.vcf.gz \
  --populations /scratch/midway3/espolston/dcgm_all.txt \
  --chromosomes $CHROM \
  --window_size 50000000 \
  --fst_components \
  --chunk_size 5000 \  #changes how much held in memory
  --output_folder /scratch/midway3/espolston/allpair_50Mb_${CHROM} \
  --n_cores 8
done

#0 AND 4 FOLD DEGENERACY SITES
#/project/kreiner/data/genome/degeneracy has 0 and 4 fold sites annotated with https://github.com/harvardinformatics/degenotate
pixy --stats pi \
  --vcf /project/kreiner/pairedenv_commongarden/allsites/droughtcgmerged_allsites_final_ID_IBD.vcf.gz \
  --populations /scratch/midway3/espolston/dcgm_namepopfile.txt \
  --bed_file /project/kreiner/data/genome/degeneracy/0fold_Atub193_hap2.bed \
  --output_folder /scratch/midway3/espolston/zerofold_pixy/ \
  --n_core 5
  
pixy --stats pi \
  --vcf /project/kreiner/pairedenv_commongarden/allsites/droughtcgmerged_allsites_final_ID_IBD.vcf.gz \
  --populations /scratch/midway3/espolston/dcgm_namepopfile.txt \
  --bed_file /project/kreiner/data/genome/degeneracy/4fold_Atub193_hap2.bed \
  --output_folder /scratch/midway3/espolston/fourfold_pixy/ \
  --n_core 5

pixy --stats watterson_theta \
  --vcf /project/kreiner/pairedenv_commongarden/allsites/droughtcgmerged_allsites_final_ID_IBD.vcf.gz \
  --populations /scratch/midway3/espolston/dcgm_namepopfile.txt \
  --bed_file /project/kreiner/data/genome/degeneracy/4fold_Atub193_hap2.bed \
  --output_folder /scratch/midway3/espolston/fourfold_pixy/ \
  --n_core 5

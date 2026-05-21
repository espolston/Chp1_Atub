#!/bin/bash
#SBATCH --job-name=selscan
#SBATCH --output=selscan.out
#SBATCH --error=selscan.err
#SBATCH --time=2-00:00:00
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5

#Download and install selscan 
#module load gcc/12.1.0
#module load python/3.10.5
#module load miniconda3/24.9.2
#conda create -n selscan_env selscan -c bioconda -c conda-forge
#conda activate selscan_env
#conda install bioconda::selscan #for some reason selscan didnt install first

module load gcc/12.1.0
module load miniconda3/24.9.2
source /apps/software/gcc-12.1.0/miniconda3/24.9.2/etc/profile.d/conda.sh
conda activate shapeit_env

#Make plink map file and integrate the LDbased map (done on midway)
#plink --bfile binary_dcgm --bim binary_dcgm.bim --fam dcgm.fam --recode --tab -allow-no-sex --allow-extra-chr --out binary_dcgm.map
#awk '{OFS="\t"} NR > 1 {$4 = split($2,a, "_"); $5 = a[2]":"$1; print $1,$2,$3,$5}' LDbased_cMestimates_bychrom.txt > cMestimate_wID.txt
#awk 'BEGIN{OFS="\t"} NR==FNR { a[$4] = $3; next} {if ($2 in a) {$3 = $3 + a[$2]} print}' cMestimate_wID.txt binary_dcgm.map.map > dcgm.map

#Make pop files from pixy dcgm_namepopfile.txt
#awk '{OFS="\t"} $2 == "AG" {print $1}' dcgm_namepopfile.txt > ag_indivs.txt
#awk '{OFS="\t"} $2 == "NAT" {print $1}' dcgm_namepopfile.txt > nat_indivs.txt

#map file and file with pops separated (ag_indivs.txt, nat_indivs.txt) copied over to randi

#XPEHH https://github.com/szpiech/selscan/
#split each scaffold vcf into ag and nat pop vcf
#convert drought samples into AT# that is in the phased vcf header
#awk '{OFS="\t"} {print $2,$1}' rehead_forwhatshap.txt > rehead_afterwhatshap.txt
#awk '{OFS="\t"} NR == FNR {a[$1]; next} $2 in a' ag_indivs.txt rehead_afterwhatshap.txt > ag_AT_header.txt
#awk '{OFS="\t"} NR == FNR {a[$1]; next} $2 in a' nat_indivs.txt rehead_afterwhatshap.txt > nat_AT_header.txt
#awk '{OFS="\t"} {print $1}' ag_AT_header.txt > ag_AT_indivs.txt
#awk '{OFS="\t"} {print $1}' nat_AT_header.txt > nat_AT_indivs.txt

parallel -j 5 'bcftools view --samples-file /scratch/espolston/sel_scan/ag_AT_indivs.txt -Oz -o /scratch/espolston/sel_scan/phased_ag_oldhead_Scaffold_{}.vcf.gz /scratch/espolston/phased/phased_dcgm_Scaffold_{}.vcf.gz' ::: {1..16}
parallel -j 5 'bcftools view --samples-file /scratch/espolston/sel_scan/nat_AT_indivs.txt -Oz -o /scratch/espolston/sel_scan/phased_nat_oldhead_Scaffold_{}.vcf.gz /scratch/espolston/phased/phased_dcgm_Scaffold_{}.vcf.gz' ::: {1..16}
parallel -j 5 'bcftools reheader --samples /scratch/espolston/sel_scan/ag_AT_header.txt -o /scratch/espolston/sel_scan/phased_ag_Scaffold_{}.vcf.gz /scratch/espolston/sel_scan/phased_ag_oldhead_Scaffold_{}.vcf.gz' ::: {1..16}
parallel -j 5 'bcftools reheader --samples /scratch/espolston/sel_scan/nat_AT_header.txt -o /scratch/espolston/sel_scan/phased_nat_Scaffold_{}.vcf.gz /scratch/espolston/sel_scan/phased_nat_oldhead_Scaffold_{}.vcf.gz' ::: {1..16}
parallel -j 5 'tabix -p vcf /scratch/espolston/sel_scan/phased_nat_Scaffold_{}.vcf.gz' ::: {1..16}
parallel -j 5 'tabix -p vcf /scratch/espolston/sel_scan/phased_ag_Scaffold_{}.vcf.gz' ::: {1..16} 

rm /scratch/espolston/sel_scan/phased_nat_oldhead_Scaffold_*.vcf.gz
rm /scratch/espolston/sel_scan/phased_ag_oldhead_Scaffold_*.vcf.gz

#make chrom specific map file of the nonzero cM sites
#awk '{OFS="\t"} $3 != 0 {print}' /scratch/espolston/sel_scan/dcgm.map > /scratch/espolston/sel_scan/dcgm_filtered.map
#for scaf in {1..16}; do awk -v s="$scaf" 'BEGIN{OFS="\t"} $1 == "Scaffold_" s {print}' /scratch/espolston/sel_scan/dcgm_filtered.map > /scratch/espolston/sel_scan/dcgm_Scaffold_${scaf}.map; done

#make chrom specific map file
for scaf in {1..16}; do awk -v s="$scaf" 'BEGIN{OFS="\t"} $1 == "Scaffold_" s {print}' /scratch/espolston/sel_scan/dcgm_cM.map > /scratch/espolston/sel_scan/dcgm_Scaffold_${scaf}.map; done

#sort the map file
#for scaf in {1..16}; do awk '{print $0}' /scratch/espolston/sel_scan/dcgm_Scaffold_${scaf}.map | sort -k4n > /scratch/espolston/sel_scan/dcgm_Scaffold_${scaf}_sorted.map; done
#for scaf in {1..16}; do rm /scratch/espolston/sel_scan/dcgm_Scaffold_${scaf}.map; done


conda activate selscan_env

#XPEHH
#vcf = vcf pop1
#vcfref = vcf pop2
#selscan --xpehh --vcf <pop1 vcf> --vcf-ref <pop2 vcf> --map <chrom specific map file from plink> --out <outfile starter> --threads <num threads>
#ehh cutoff, the haplotype homozygosity decays from 1 to 0 as moves away from focal loci, once the cutoff value is reached the ehh decay curve is truncated and integrated
parallel -j 5 'selscan --xpehh --vcf /scratch/espolston/sel_scan/phased_ag_Scaffold_{}.vcf.gz --vcf-ref /scratch/espolston/sel_scan/phased_nat_Scaffold_{}.vcf.gz --map /scratch/espolston/sel_scan/dcgm_Scaffold_{}.map --cutoff .05 --out /scratch/espolston/sel_scan/Scaffold_{} --threads 5' ::: {1..16}

#run norm after to normalize and analyze selection outputs
#read pg 19 https://github.com/szpiech/selscan/blob/main/Selscan_Manual_v2_1_1.pdf for more options on filtering of this
#norm --xpehh --files <outfiles produced, can run all at once> --log <outputname>
norm --xpehh --files sel_scan/Scaffold_*.xpehh.out  --log sel_scan/xpehh_norm_1_16 --bp-win
norm --xpehh --files sel_scan/Scaffold_*.xpehh.out  --log sel_scan/xpehh_norm_1_16


#running per pair

#use pair indiv lists from pixy
for p in {1..17}; do awk '{OFS="\t"} $2 == "AG" {print $1}' /scratch/espolston/sel_scan_pair/dcgm_pair_${p}.txt > /scratch/espolston/sel_scan_pair/dcgm_pair_${p}_ag.txt; done
for p in {1..17}; do awk '{OFS="\t"} $2 == "NAT" {print $1}' /scratch/espolston/sel_scan_pair/dcgm_pair_${p}.txt > /scratch/espolston/sel_scan_pair/dcgm_pair_${p}_nat.txt; done

#create pair ag/nat vcfs
parallel -j 5 'bcftools view --samples-file /scratch/espolston/sel_scan_pair/dcgm_pair_{2}_nat.txt -Oz -o /scratch/espolston/sel_scan_pair/phased_dcgm_Scaffold_{1}_nat_pair_{2}.vcf.gz /scratch/espolston/sel_scan/phased_nat_Scaffold_{1}.vcf.gz' ::: {1..16} ::: {1..17}
parallel -j 5 'bcftools view --samples-file /scratch/espolston/sel_scan_pair/dcgm_pair_{2}_ag.txt -Oz -o /scratch/espolston/sel_scan_pair/phased_dcgm_Scaffold_{1}_ag_pair_{2}.vcf.gz /scratch/espolston/sel_scan/phased_ag_Scaffold_{1}.vcf.gz' ::: {1..16} ::: {1..17}
parallel -j 5 'tabix -p vcf /scratch/espolston/sel_scan_pair/phased_dcgm_Scaffold_{1}_nat_pair_{2}.vcf.gz'  ::: {1..16} ::: {1..17}
parallel -j 5 'tabix -p vcf /scratch/espolston/sel_scan_pair/phased_dcgm_Scaffold_{1}_ag_pair_{2}.vcf.gz'  ::: {1..16} ::: {1..17} 

conda activate selscan_env

#XPEHH
#vcf = vcf pop1
#vcfref = vcf pop2
parallel -j 5 'selscan --xpehh --vcf /scratch/espolston/sel_scan_pair/phased_dcgm_Scaffold_{1}_ag_pair_{2}.vcf.gz --vcf-ref /scratch/espolston/sel_scan_pair/phased_dcgm_Scaffold_{1}_nat_pair_{2}.vcf.gz --map /scratch/espolston/sel_scan/dcgm_Scaffold_{1}.map --cutoff .05 --out /scratch/espolston/sel_scan_pair/Scaffold_{1}_pair_{2} --threads 5' ::: {1..16}  ::: {1..17}

#run norm after to normalize and analyze selection outputs
norm --xpehh --files sel_scan_pair/Scaffold_*.xpehh.out  --log sel_scan_pair/xpehh_norm_1_16

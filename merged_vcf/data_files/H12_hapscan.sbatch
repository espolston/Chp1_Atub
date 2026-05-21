#!/bin/bash
#SBATCH --job-name=H12_hapscan
#SBATCH --output=H12_hapscan.out
#SBATCH --error=H12_hapscan.err
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=espolston@rcc.uchicago.edu 
#SBATCH --account=pi-kreiner

module load python/anaconda-2022.05
conda activate /project/ecev35300/steinrue/conda_env/week6_selection_env

#for some reason didnt have phased vcf index files
for sc in {1..16}; do
tabix /scratch/midway3/espolston/phased/phased_dcgm_Scaffold_${sc}.vcf.gz; done

#create input file: pos, ind1_nucleotide, ind2_nucleotide ... (with any non ATCG = N)

#run H12
#https://github.com/ngarud/SelectionHapStats/blob/master/README_selectionHapStats.pdf 

#Make file of POS and then comma separated list of all indivs nucleotides at that position (phased) with anyother letters (Not ATCG) converted to an N
for sc in {1..16}; do
bcftools query --format '%POS[,%TGT]\n' /scratch/midway3/espolston/phased/phased_dcgm_Scaffold_${sc}.vcf.gz > /scratch/midway3/espolston/garud_h12/h12_input_Scaffold_${sc}.txt; done

for sc in {1..16}; do awk -F',' '{
	printf $1;	
    for (n=2; n<=NF; ++n) {
        split($n, part, "|");
        printf ","part[1]","part[2];
    }
    printf "\n";
}' /scratch/midway3/espolston/garud_h12/h12_input_Scaffold_${sc}.txt > /scratch/midway3/espolston/garud_h12/h12_final_Scaffold_${sc}.txt; done

#Run H12 and H2H1
#python H12_H2H1.py <file of all nucleotides at each site> <sample size> -o <outfilename.txt> -w <window site in snps, set at 400 in ex> -j <distance between centers of analysis windows in snps, default 50> -d <tolerable hamming distance between haplotypes to be single haplotype, default 0>
for sc in {1..16}; do python /scratch/midway3/espolston/h12/H12_H2H1.py  /scratch/midway3/espolston/garud_h12/h12_final_Scaffold_${sc}.txt 440 -o /scratch/midway3/espolston/garud_h12/dcgm_hap_stats_Scaffold_${sc}.txt -w 400; done

#detect the h12 peaks
#python H12peakFinder.py <output from h12_h1h2.py> -o <output name> can also use -t to define h12 cutoff
for sc in {1..16}; do python /scratch/midway3/espolston/h12/H12peakFinder.py dcgm_hap_stats_Scaffold_${sc}.txt -o /scratch/midway3/espolston/garud_h12/dcgm_hap_peaks_Scaffold_${sc}.txt; done

#run visualizations - there are more as well https://github.com/ngarud/SelectionHapStats/blob/master/README_selectionHapStats.pdf 
#ran these locally
for sc in {1..16}; do Rscript /Users/libbypolston/Desktop/h12/H12_viz.R /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12_H2H1/dcgm_hap_stats_Scaffold_${sc}.txt /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12_H2H1/dcgm_hap_peaks_Scaffold_${sc}.txt
/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/H12_H2H1/dcgm_Scaffold_${sc}_h12scan.pdf 10; done

for sc in {1..16}; do Rscript /Users/libbypolston/Desktop/h12/hapSpectrum_viz.R /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12_H2H1/dcgm_hap_peaks_Scaffold_${sc}.txt /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/H12_H2H1/dcgm_Scaffold_${sc}_spectrum.pdf 10 440; done

#run h12 PER PAIR!
#first made phased vcfs have the updated sample names
module load python/anaconda-2022.05
source activate /project/kreiner/espolston/

parallel -j 5 'bcftools reheader --samples /scratch/midway3/espolston/rehead_afterwhatshap.txt -o /scratch/midway3/espolston/new_phased/phased_dcgm_Scaffold_{}_newhead.vcf.gz /scratch/midway3/espolston/phased/phased_dcgm_Scaffold_{}.vcf.gz' ::: {1..16}
parallel -j 5 'tabix -p vcf /scratch/midway3/espolston/new_phased/phased_dcgm_Scaffold_{}_newhead.vcf.gz' ::: {1..16}

#split ag and nat vcfs - with dcgm_namepopfile.txt
awk '{OFS="\t"} $2 == "AG" {print $1}' dcgm_namepopfile.txt > ag_indivs.txt
awk '{OFS="\t"} $2 == "NAT" {print $1}' dcgm_namepopfile.txt > nat_indivs.txt

parallel -j 5 'bcftools view --samples-file ag_indivs.txt -Oz -o /scratch/midway3/espolston/garud_h12_agnat/phased_dcgm_Scaffold_{}_ag.vcf.gz /scratch/midway3/espolston/new_phased/phased_dcgm_Scaffold_{}_newhead.vcf.gz' ::: {1..16}
parallel -j 5 'tabix -p vcf /scratch/midway3/espolston/garud_h12_agnat/phased_dcgm_Scaffold_{}_ag.vcf.gz' ::: {1..16}
parallel -j 5 'bcftools view --samples-file nat_indivs.txt -Oz -o /scratch/midway3/espolston/garud_h12_agnat/phased_dcgm_Scaffold_{}_nat.vcf.gz /scratch/midway3/espolston/new_phased/phased_dcgm_Scaffold_{}_newhead.vcf.gz' ::: {1..16}
parallel -j 5 'tabix -p vcf /scratch/midway3/espolston/garud_h12_agnat/phased_dcgm_Scaffold_{}_nat.vcf.gz' ::: {1..16}

for sc in {1..16}; do
bcftools query --format '%POS[,%TGT]\n' /scratch/midway3/espolston/garud_h12_agnat/phased_dcgm_Scaffold_${sc}_ag.vcf.gz > /scratch/midway3/espolston/garud_h12_agnat/h12_input_Scaffold_${sc}_ag.txt; done

for sc in {1..16}; do
bcftools query --format '%POS[,%TGT]\n' /scratch/midway3/espolston/garud_h12_agnat/phased_dcgm_Scaffold_${sc}_nat.vcf.gz > /scratch/midway3/espolston/garud_h12_agnat/h12_input_Scaffold_${sc}_nat.txt; done

for sc in {1..16}; do awk -F',' '{
	printf $1;	
    for (n=2; n<=NF; ++n) {
        split($n, part, "|");
        printf ","part[1]","part[2];
    }
    printf "\n";
}' /scratch/midway3/espolston/garud_h12_agnat/h12_input_Scaffold_${sc}_ag.txt > /scratch/midway3/espolston/garud_h12_agnat/h12_final_Scaffold_${sc}_ag.txt; done

for sc in {1..16}; do awk -F',' '{
	printf $1;	
    for (n=2; n<=NF; ++n) {
        split($n, part, "|");
        printf ","part[1]","part[2];
    }
    printf "\n";
}' /scratch/midway3/espolston/garud_h12_agnat/h12_input_Scaffold_${sc}_nat.txt > /scratch/midway3/espolston/garud_h12_agnat/h12_final_Scaffold_${sc}_nat.txt; done

#get # indivs per env type
#awk '$2=="AG" {count++} END {print count}' dcgm_namepopfile.txt #219
#awk '$2=="NAT" {count++} END {print count}' dcgm_namepopfile.txt #221

module load python/anaconda-2022.05
conda activate /project/ecev35300/steinrue/conda_env/week6_selection_env

#Run H12 and H2H1
for sc in {1..16}; do python /scratch/midway3/espolston/h12/H12_H2H1.py  /scratch/midway3/espolston/garud_h12_agnat/h12_final_Scaffold_${sc}_ag.txt 219 -o /scratch/midway3/espolston/garud_h12_agnat/dcgm_hap_stats_Scaffold_${sc}_ag.txt -w 400; done

#detect the h12 peaks
for sc in {1..16}; do python /scratch/midway3/espolston/h12/H12peakFinder.py /scratch/midway3/espolston/garud_h12_agnat/dcgm_hap_stats_Scaffold_${sc}_ag.txt -o /scratch/midway3/espolston/garud_h12_agnat/dcgm_hap_peaks_Scaffold_${sc}_ag.txt; done

#Run H12 and H2H1
for sc in {1..16}; do python /scratch/midway3/espolston/h12/H12_H2H1.py  /scratch/midway3/espolston/garud_h12_agnat/h12_final_Scaffold_${sc}_nat.txt 220 -o /scratch/midway3/espolston/garud_h12_agnat/dcgm_hap_stats_Scaffold_${sc}_nat.txt -w 400; done

#detect the h12 peaks
for sc in {1..16}; do python /scratch/midway3/espolston/h12/H12peakFinder.py /scratch/midway3/espolston/garud_h12_agnat/dcgm_hap_stats_Scaffold_${sc}_nat.txt -o /scratch/midway3/espolston/garud_h12_agnat/dcgm_hap_peaks_Scaffold_${sc}_nat.txt; done

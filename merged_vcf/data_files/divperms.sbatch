#!/bin/bash
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=10
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --job-name=pixy
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=espolston@rcc.uchicago.edu 

module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner/espolston/

#make bed file for 10 kb windows across genome
#get fasta_generate_regions.py from https://github.com/freebayes/freebayes/blob/master/scripts/fasta_generate_regions.py 
#cp /cds3/kreiner/longread_data_assemblies/final_assemblies/Atub_193_hap2.fasta 
#make regions
#python3 fasta_generate_regions.py Atub_193_hap2.fasta 10000 > dcgm_10kbsites_start.bed
#make regions in correct format
#awk '{OFS="\t"; split($1, part, ":"); $2=part[1]; $3=part[2]; split($3, region, "-"); $4=region[1]; $5=region[2]; print $2,$4,$5}' dcgm_10kbsites_start.bed > dcgm_10kbsites.bed

#get AF from VCF
#bcftools query -f '%CHROM\t%POS\t%AF\n' /project/kreiner/pairedenv_commongarden/normalized_SNPsonly_vcfs/drought_commongarden_merged_final_ID_IBD.vcf.gz > allsnp_AF.txt

awk 'BEGIN{OFS="\t"} {$4 = "Scaffold_"$1; print $4,$2,$3}' dcgm_sig_xpehh.bed > dcgm_sig_xpehh_1.bed
mv dcgm_sig_xpehh_1.bed dcgm_sig_xpehh.bed

#for chr spec files
#dcm_sig.bed and dcgm_sig_xpehh.bed both made in chp1_analyses.Rmd
for sc in {1..16}; do 
	awk -v scaf=${sc} 'BEGIN{OFS="\t"} $1 == "Scaffold_" scaf {$4 = $2+1; print $1,$2,$4,$3}' allsnp_AF.txt > allsnp_AF_scaf_${sc}.txt
	awk -v scaf=${sc} 'BEGIN{OFS="\t"} $1 == "Scaffold_" scaf {print $1,$2,$3}' dcgm_10kbsites.bed > dcgm_10kbsites_${sc}.bed
	awk -v scaf=${sc} 'BEGIN{OFS="\t"} $1 == "Scaffold_" scaf {print $1,$2,$3}' dcgm_sig.bed > dcgm_sig_${sc}.bed
	awk -v scaf=${sc} 'BEGIN{OFS="\t"} $1 == "Scaffold_" scaf {print $1,$2,$3}' dcgm_sig_xpehh.bed > dcgm_sig_xpehh_${sc}.bed
	
	#Get mean AF in all 10 kb windows
	bedtools map -a dcgm_10kbsites_${sc}.bed -b allsnp_AF_scaf_${sc}.txt -c 4 -o mean > dcgm_10kbsites_${sc}_AF
	
	#Get mean AF in 10kb windows around CMH clumped sites
#dcgm_sig.bed = 10 kb window around cmh site
	bedtools map -a dcgm_sig_${sc}.bed -b allsnp_AF_scaf_${sc}.txt -c 4 -o mean > dcgm_sig_${sc}_AF
	
	bedtools map -a dcgm_sig_xpehh_${sc}.bed -b allsnp_AF_scaf_${sc}.txt -c 4 -o mean > dcgm_sig_${sc}_AF
	
	rm allsnp_AF_scaf_${sc}.txt dcgm_10kbsites_${sc}.bed dcgm_sig_${sc}.bed dcgm_sig_xpehh_${sc}.bed; done

#merge all of these files into one
cat dcgm_sig_* > dcgm_sig_all_AF.txt
cat dcgm_10kbsites_* > dcgm_10kbsites_all_AF.txt
cat dcgm_sig_xpehh_* > dcgm_sig_xpehh_all_AF.txt

#sites with no snps were "." not NA
sed 's/\t\.\t/\tNA\t/g; s/\t\.$/\tNA/' dcgm_10kbsites_all_AF.txt > dcgm_10kbsites_all_AF_final.txt
sed 's/\t\.\t/\tNA\t/g; s/\t\.$/\tNA/' dcgm_sig_all_AF.txt > dcgm_sig_all_AF_final.txt
sed 's/\t\.\t/\tNA\t/g; s/\t\.$/\tNA/' dcgm_sig_xpehh_all_AF.txt > dcgm_sig_xpehh_all_AF_final.txt


for sc in {1..16}; do rm dcgm_sig_${sc}_AF dcgm_sig_xpehh_${sc}_AF dcgm_10kbsites_${sc}_AF; done

#conda conda create -n pixy -c bioconda -c conda-forge::pixy
#conda install bioconda::bedtools

conda activate pixy
module load parallel

#get mean pi in 10kb windows that we have our AF freq for
#perl pixy_pop.pl sample_ids.txt dcgm.fam dcgm_namepopfile.txt dcgm_pair
parallel -j 10 'pixy --stats pi \
  --vcf /project/kreiner/pairedenv_commongarden/allsites/droughtcgmerged_allsites_final_ID_IBD.vcf.gz \
  --populations /scratch/midway2/espolston/pixyforperm/dcgm_pair_{1}.txt \
  --bed_file /scratch/midway2/espolston/pixyforperm/dcgm_10kbsites.bed \
  --chromosomes Scaffold_{2} \
  --output_folder /scratch/midway2/espolston/pixyforperm/pair_{1}_scaf_{2}_10kbwin \
  --chunk_size 5000 \
  --n_core 5' ::: {1..17} ::: {1..16}

#get mean pi in 10kb windows for cmh sites
parallel -j 10 'pixy --stats pi \
  --vcf /project/kreiner/pairedenv_commongarden/allsites/droughtcgmerged_allsites_final_ID_IBD.vcf.gz \
  --populations /scratch/midway2/espolston/pixyforperm/dcgm_pair_{1}.txt \
  --bed_file /scratch/midway2/espolston/pixyforperm/dcgm_sig.bed \
  --chromosomes Scaffold_{2} \
  --output_folder /scratch/midway2/espolston/pixyforperm/pair_{1}_scaf_{2}_cmh \
  --chunk_size 5000 \
  --n_core 5' ::: {1..17} ::: {1..16}
  
parallel -j 10 'pixy --stats pi \
  --vcf /project/kreiner/pairedenv_commongarden/allsites/droughtcgmerged_allsites_final_ID_IBD.vcf.gz \
  --populations /scratch/midway2/espolston/pixyforperm/dcgm_pair_{1}.txt \
  --bed_file /scratch/midway2/espolston/pixyforperm/dcgm_sig_xpehh.bed \
  --chromosomes Scaffold_{2} \
  --output_folder /scratch/midway2/espolston/pixyforperm/pair_{1}_scaf_{2}_xpehh \
  --chunk_size 5000 \
  --n_core 5' ::: {1..17} ::: {1..16}

for pnum in {1..17}; do
    for snum in {1..16}; do
        dir="/scratch/midway2/espolston/pixyforperm/pair_${pnum}_scaf_${snum}_10kbwin"
        for file in "$dir"/*.txt; do
            [ -f "$file" ] || continue  # skip if no .txt files found
            base=$(basename "$file")
            mv "$file" "/scratch/midway2/espolston/pixyforperm/finished_10kbwin/${base%.txt}_pair_${pnum}_scaf_${snum}_10kbwin.txt"
            rm -rf $dir
        done
    done
done

for pnum in {1..17}; do
    for snum in {1..16}; do
        dir="/scratch/midway2/espolston/pixyforperm/pair_${pnum}_scaf_${snum}_cmh"
        for file in "$dir"/*.txt; do
            [ -f "$file" ] || continue  # skip if no .txt files found
            base=$(basename "$file")
            mv "$file" "/scratch/midway2/espolston/pixyforperm/finished_cmh/${base%.txt}_pair_${pnum}_scaf_${snum}_cmh.txt"
            rm -rf $dir
        done
    done
done

for pnum in {1..17}; do
    for snum in {1..16}; do
        dir="/scratch/midway2/espolston/pixyforperm/pair_${pnum}_scaf_${snum}_xpehh"
        for file in "$dir"/*.txt; do
            [ -f "$file" ] || continue  # skip if no .txt files found
            base=$(basename "$file")
            mv "$file" "/scratch/midway2/espolston/pixyforperm/finished_xpehh/${base%.txt}_pair_${pnum}_scaf_${snum}_xpehh.txt"
            rm -rf $dir
        done
    done
done

#Choose 865 of these windows that have the same AF (within + or - .01) and calculate slope, mean dif in ag v nat pi and errors of each
#did in pairs in case I want to break out more

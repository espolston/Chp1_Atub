#!/bin/bash
#SBATCH --job-name=shapeit
#SBATCH --output=shapeitprep.out
#SBATCH --error=shapeitprep.err
#SBATCH --time=2-00:00:00
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=50GB
#SBATCH --cpus-per-task=1

#for randi
module load gcc/12.1.0
module load miniconda3/24.9.2
source /apps/software/gcc-12.1.0/miniconda3/24.9.2/etc/profile.d/conda.sh
conda activate shapeit_env

#filter recomb map to just that scaffold before shapeit
awk '{OFS="\t"} $2 == "Scaffold_1" {print}' /scratch/espolston/LDbased_cMestimates_bychrom.txt > /scratch/espolston/LDbased_cMestimates_bychrom_Scaffold_1.txt
awk '{OFS="\t"} $2 == "Scaffold_2" {print}' /scratch/espolston/LDbased_cMestimates_bychrom.txt > /scratch/espolston/LDbased_cMestimates_bychrom_Scaffold_2.txt
awk '{OFS="\t"} $2 == "Scaffold_3" {print}' /scratch/espolston/LDbased_cMestimates_bychrom.txt > /scratch/espolston/LDbased_cMestimates_bychrom_Scaffold_3.txt
awk '{OFS="\t"} $2 == "Scaffold_4" {print}' /scratch/espolston/LDbased_cMestimates_bychrom.txt > /scratch/espolston/LDbased_cMestimates_bychrom_Scaffold_4.txt
awk '{OFS="\t"} $2 == "Scaffold_5" {print}' /scratch/espolston/LDbased_cMestimates_bychrom.txt > /scratch/espolston/LDbased_cMestimates_bychrom_Scaffold_5.txt
awk '{OFS="\t"} $2 == "Scaffold_6" {print}' /scratch/espolston/LDbased_cMestimates_bychrom.txt > /scratch/espolston/LDbased_cMestimates_bychrom_Scaffold_6.txt
awk '{OFS="\t"} $2 == "Scaffold_7" {print}' /scratch/espolston/LDbased_cMestimates_bychrom.txt > /scratch/espolston/LDbased_cMestimates_bychrom_Scaffold_7.txt
awk '{OFS="\t"} $2 == "Scaffold_8" {print}' /scratch/espolston/LDbased_cMestimates_bychrom.txt > /scratch/espolston/LDbased_cMestimates_bychrom_Scaffold_8.txt
awk '{OFS="\t"} $2 == "Scaffold_9" {print}' /scratch/espolston/LDbased_cMestimates_bychrom.txt > /scratch/espolston/LDbased_cMestimates_bychrom_Scaffold_9.txt
awk '{OFS="\t"} $2 == "Scaffold_10" {print}' /scratch/espolston/LDbased_cMestimates_bychrom.txt > /scratch/espolston/LDbased_cMestimates_bychrom_Scaffold_10.txt
awk '{OFS="\t"} $2 == "Scaffold_11" {print}' /scratch/espolston/LDbased_cMestimates_bychrom.txt > /scratch/espolston/LDbased_cMestimates_bychrom_Scaffold_11.txt
awk '{OFS="\t"} $2 == "Scaffold_12" {print}' /scratch/espolston/LDbased_cMestimates_bychrom.txt > /scratch/espolston/LDbased_cMestimates_bychrom_Scaffold_12.txt
awk '{OFS="\t"} $2 == "Scaffold_13" {print}' /scratch/espolston/LDbased_cMestimates_bychrom.txt > /scratch/espolston/LDbased_cMestimates_bychrom_Scaffold_13.txt
awk '{OFS="\t"} $2 == "Scaffold_14" {print}' /scratch/espolston/LDbased_cMestimates_bychrom.txt > /scratch/espolston/LDbased_cMestimates_bychrom_Scaffold_14.txt
awk '{OFS="\t"} $2 == "Scaffold_15" {print}' /scratch/espolston/LDbased_cMestimates_bychrom.txt > /scratch/espolston/LDbased_cMestimates_bychrom_Scaffold_15.txt
awk '{OFS="\t"} $2 == "Scaffold_16" {print}' /scratch/espolston/LDbased_cMestimates_bychrom.txt > /scratch/espolston/LDbased_cMestimates_bychrom_Scaffold_16.txt


#need to make single scaffold vcf
#cd /scratch/espolston/whatshap_Scaffold_1
#for file in *.gz; do tabix -p vcf "$file"; done 

#cd /scratch/espolston/whatshap_Scaffold_2
#for file in *.gz; do tabix -p vcf "$file"; done

#cd /scratch/espolston/whatshap_Scaffold_3
#for file in *.gz; do tabix -p vcf "$file"; done

#cd /scratch/espolston/whatshap_Scaffold_4
#for file in *.gz; do tabix -p vcf "$file"; done

#cd /scratch/espolston/whatshap_Scaffold_5
#for file in *.gz; do tabix -p vcf "$file"; done

#cd /scratch/espolston/whatshap_Scaffold_6
#for file in *.gz; do tabix -p vcf "$file"; done

#cd /scratch/espolston/whatshap_Scaffold_7
#for file in *.gz; do tabix -p vcf "$file"; done

#cd /scratch/espolston/whatshap_Scaffold_8
#for file in *.gz; do tabix -p vcf "$file"; done

#cd /scratch/espolston/whatshap_Scaffold_9
#for file in *.gz; do tabix -p vcf "$file"; done

#cd /scratch/espolston/whatshap_Scaffold_10
#for file in *.gz; do tabix -p vcf "$file"; done

#cd /scratch/espolston/whatshap_Scaffold_11
#for file in *.gz; do tabix -p vcf "$file"; done

#cd /scratch/espolston/whatshap_Scaffold_12
#for file in *.gz; do tabix -p vcf "$file"; done

#cd /scratch/espolston/whatshap_Scaffold_13
#for file in *.gz; do tabix -p vcf "$file"; done

#cd /scratch/espolston/whatshap_Scaffold_14
#for file in *.gz; do tabix -p vcf "$file"; done

#cd /scratch/espolston/whatshap_Scaffold_15
#for file in *.gz; do tabix -p vcf "$file"; done

#cd /scratch/espolston/whatshap_Scaffold_16
#for file in *.gz; do tabix -p vcf "$file"; done

#need to remove GL bc whatshap messes up GL field, then merge
parallel -j 5 'cd /scratch/espolston/whatshap_Scaffold_{}/; for file in *.gz; do bcftools annotate -x FORMAT/GL "$file" -Oz -o noGL_"$file"; done' ::: {1..16}

parallel -j 5 'cd /scratch/espolston/whatshap_Scaffold_{}/; for file in noGL_*.gz; do tabix -p vcf "$file"; done' ::: {1..16}

parallel -j 5 'cd /scratch/espolston/whatshap_Scaffold_{}/; bcftools merge noGL_whatshapout_*_Scaffold_{}.vcf.gz -Oz -o /scratch/espolston/whatshapout_Scaffold_{}.vcf.gz' ::: {1..16}

parallel -j 5 'tabix -p vcf /scratch/espolston/whatshapout_Scaffold_{}.vcf.gz' ::: {1..16}

!/bin/bash
#SBATCH --job-name=IBD
#SBATCH --output=IBD.out
#SBATCH --error=IBD.err
#SBATCH --time=5:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=espolston@rcc.uchicago.edu 

module load plink

#LD thin before IBD calcs
#https://www.cog-genomics.org/plink/1.9/ld#indep
#https://zzz.bwh.harvard.edu/plink/summary.shtml
#window size in SNPs (e.g. 50), the number of SNPs to shift the window at each step (e.g. 5), the VIF threshold. The VIF is 1/(1-R^2) where R^2 is the multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously

plink --bfile /scratch/midway3/espolston/vcf/binary_dcgm --fam /scratch/midway3/espolston/vcf/dcgm.fam --allow-no-sex --allow-extra-chr --make-founders --indep 50 5 2 

plink --bfile /scratch/midway3/espolston/vcf/binary_dcgm --fam /scratch/midway3/espolston/vcf/dcgm.fam --extract plink.prune.in --allow-no-sex --allow-extra-chr  --make-bed --out /scratch/midway3/espolston/vcf/dcgm_thinned

#IBD calcs
#https://www.cog-genomics.org/plink/1.9/ibd

plink --bfile /scratch/midway3/espolston/vcf/dcgm_thinned --fam /scratch/midway3/espolston/vcf/dcgm.fam --allow-no-sex --allow-extra-chr --make-founders --genome gz --min .05 --out /scratch/midway3/espolston/vcf/dcgm_IBD

cp /cds3/kreiner/dcgm_filteredvcf/plink_pixy/binary_dcgm.* /scratch/midway3/espolston/
cp /cds3/kreiner/dcgm_filteredvcf/plink_pixy/dcgm.fam /scratch/midway3/espolston/
scp espolston@midway3.rcc.uchicago.edu:/scratch/midway3/espolston/dcgm.fam espolston@randi.cri.uchicago.edu:/scratch/espolston/
scp espolston@midway3.rcc.uchicago.edu:/scratch/midway3/espolston/binary_dcgm.\* espolston@randi.cri.uchicago.edu:/scratch/espolston/

module load gcc/11.3.0
module load gcc/12.1.0
module load intel/2022.2
module load llvm/14.0.5
module load plink/2.0

#estimating relatedness from IBD for afvaper results
plink2 --bfile /scratch/espolston/related/binary_dcgm --bim /scratch/espolston/related/binary_dcgm.bim --fam /scratch/espolston/related/dcgm.fam --make-king square0 --allow-extra-chr --out /scratch/espolston/related/dcgm_relatednessIBD

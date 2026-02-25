#!/bin/bash
#SBATCH --job-name=shapeit
#SBATCH --output=shapeit.out
#SBATCH --error=shapeit.err
#SBATCH --time=5-00:00:00
#SBATCH --partition=tier1q
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4

#recombination map from ~/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Code/rhotocM_forlibby.R and allchroms_193reference_ldmapresults.txt2

#module load python/anaconda-2023.09
#conda config --add channels defaults
#conda config --add channels bioconda
#conda config --add channels conda-forge
#conda config --set channel_priority strict
#conda create -n shapeit_env shapeit4
#conda activate shapeit_env
#conda install -c conda-forge bcftools

#Installing on randi
#module load gcc/12.1.0
#module load python/3.10.5
#module load miniconda3/24.9.2
#conda create -n shapeit_env shapeit4 -c bioconda -c conda-forge
#conda activate shapeit_env 
#conda install -c conda-forge bcftools

#recomb map made in R with rhotocM_forlibby.R and uses allchroms_193reference_ldmapresults.txt2

#for midway
#module load python/anaconda-2023.09
#source activate shapeit_env
#module load parallel
#module load apptainer

#for randi
module load gcc/12.1.0
module load miniconda3/24.9.2
source /apps/software/gcc-12.1.0/miniconda3/24.9.2/etc/profile.d/conda.sh
conda activate shapeit_env

# ============================================
# CONFIGURATION - Edit these paths
# ============================================
SAMPLES_FILE="/scratch/espolston/scafs.txt"
MAP_FILE="/scratch/espolston/LDbased_cMestimates_bychrom"
WHATSHAP_DIR_START="/scratch/espolston/"
OUT_DIR_START="/scratch/espolston/phased"
# ============================================

#prep of making scaffold pos files and merging whatshap output in shapeitprep.sbatch

run_shapeit() {
	region=$1
	
	mkdir -p "${OUT_DIR_START}"
	mkdir -p "${OUT_DIR_START}/completed/"
	
	done_file="${OUT_DIR_START}/completed/${samp}_${region}.done"
	
	#Run shapeit on that chromosome
	#from this: https://odelaneau.github.io/shapeit4/#installation
	if shapeit4 --input "${WHATSHAP_DIR_START}/whatshapout_${region}.vcf.gz" --map "${MAP_FILE}_${region}.txt" --region "${region}" --use-PS 0.0001 --thread "${SLURM_CPUS_PER_TASK}" --output "${OUT_DIR_START}/phased_dcgm_${region}.vcf"; then
        # Mark as complete only if successful
        echo "Completed at $(date)" > "${done_file}"
        echo "Region: ${region}" >> "${done_file}"
    fi
    
    bgzip -c "${OUT_DIR_START}_${region}/phased_dcgm_${region}.vcf" > "${OUT_DIR_START}/phased_dcgm_${region}.vcf.gz"
    
    tabix -p vcf "${OUT_DIR_START}/phased_dcgm_${region}.vcf.gz"

    rm "${OUT_DIR_START}/phased_dcgm_${region}.vcf"
}

export -f run_shapeit
export SAMPLES_FILE WHATSHAP_DIR_START OUT_DIR_START MAP_FILE

# one srun per job; -c matches --cpus-per-task
SRUN="srun -N1 -n1 -c ${SLURM_CPUS_PER_TASK}"

# parallel caps total concurrency to your allocation size
# $SLURM_NTASKS equals nodes * ntasks-per-node
parallel --delay 0.2 \
--colsep '\t' \
--j "${SLURM_NTASKS}" \
--joblog "runtask-${SLURM_JOBID}.log" \
--resume \
--env SAMPLES_FILE,MAP_FILE,OUT_DIR_START,WHATSHAP_DIR_START,run_shapeit \
"${SRUN} bash -lc 'source /apps/software/gcc-12.1.0/miniconda3/24.9.2/etc/profile.d/conda.sh
conda activate shapeit_env
run_shapeit {1}'" \
:::: "${SAMPLES_FILE}"

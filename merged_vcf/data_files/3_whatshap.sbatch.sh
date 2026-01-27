#!/bin/bash
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=12
#SBATCH --time=0:10:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --job-name indiv_gatk_parallel
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=espolston@rcc.uchicago.edu 

module load python/anaconda-2023.09
source activate whatshap-env
module load parallel
module load apptainer

# ============================================
# CONFIGURATION - Edit these paths
# ============================================
REFERENCE="/scratch/midway3/espolston/Atub_193_hap2.fasta"
SAMPLES_FILE="/scratch/midway3/espolston/failed_regions_whatshap.txt"
MAIN_DIR="/scratch/midway3/espolston/"
OUT_DIR_START="/scratch/midway2/espolston/whatshap"
BAM_DIR="/scratch/midway2/espolston/bams_commongardenanddrought/"
# ============================================

mkdir -p "${OUT_DIR_START}"
run_whatshap() {
    samp=$1
	bam=$2
	region=$3
	
	mkdir -p "${OUT_DIR_START}_${region}"
	
	#check if vcf was already made, otherwise make sample/scaffold specific vcf to speed up whatshap
	if [ -f "${MAIN_DIR}sample_${samp}_${region}.vcf.gz" ]; then
    	echo "files already created"
	else
    	bcftools view --samples "${samp}" --regions "${region}" /project/kreiner/pairedenv_commongarden/normalized_SNPsonly_vcfs/drought_commongarden_whatshap.vcf.gz -Oz -o "${MAIN_DIR}sample_${samp}_${region}.vcf.gz"
    	
    	tabix -p vcf "${MAIN_DIR}sample_${samp}_${region}.vcf.gz"
	fi

	done_file="${OUT_DIR_START}/${samp}_${region}.done"
	
	#run whatshap
	cd ${BAM_DIR}
	if whatshap phase --reference="${REFERENCE}" -o "${OUT_DIR_START}_${region}/whatshapout_${samp}_${region}.vcf" --sample="${samp}" --chromosome="${region}" "${MAIN_DIR}sample_${samp}_${region}.vcf.gz" "${BAM_DIR}${bam}"; then
        # Mark as complete only if successful
        echo "Completed at $(date)" > "${done_file}"
        echo "Region: ${region}" >> "${done_file}"
    fi
    
	rm "${MAIN_DIR}sample_${samp}_${region}.vcf.gz"
	rm "${MAIN_DIR}sample_${samp}_${region}.vcf.gz.tbi"
}

export -f run_whatshap
export REFERENCE SAMPLES_FILE MAIN_DIR OUT_DIR_START BAM_DIR

# one srun per job; -c matches --cpus-per-task
SRUN="srun -N1 -n1 -c ${SLURM_CPUS_PER_TASK}"

# parallel caps total concurrency to your allocation size
# $SLURM_NTASKS equals nodes * ntasks-per-node
parallel --delay 0.2 \
--colsep '\t' \
--j "${SLURM_NTASKS}" \
--joblog "runtask-${SLURM_JOBID}.log" \
--resume \
--env REFERENCE,SAMPLES_FILE,MAIN_DIR,OUT_DIR_START,BAM_DIR,run_whatshap \
"${SRUN} bash -lc 'run_whatshap {1} {2} {3}'" \
:::: "${SAMPLES_FILE}"
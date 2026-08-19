#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --job-name gatk_vcf_parallel
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=espolston@rcc.uchicago.edu 

module load python/anaconda-2023.09
source activate /project/kreiner/espolston/
module load gcc
module load parallel
module load samtools/
module load apptainer

# ============================================
# CONFIGURATION - Edit these paths
# ============================================
REFERENCE="/scratch/midway3/espolston/Atub_193_hap2.fasta"
FAILED_REGIONS="/scratch/midway3/espolston/failed_regions_vcf.txt"
SAMPLES_FILE="/scratch/midway3/espolston/gvcf_names.txt"
INPUT_DIR="/scratch/midway2/espolston/gvcfs"
OUTPUT_DIR="/scratch/midway3/espolston"
# ============================================

run_gatk_vcf() {
    region=$1
    region_safe=$(echo "$region" | tr ':' '_' | tr '-' '_')
    
    #make map file for gatk
	#form of: sample1\tPATH/sample1.g.vcf.gz
	mkdir -p "${OUTPUT_DIR}/db"

#Make region specific databases
apptainer exec /project/kreiner/espolston/gatk_4.6.2.0.sif gatk GenomicsDBImport \
       --genomicsdb-workspace-path "${OUTPUT_DIR}/db/${region_safe}" \
       --batch-size 50 \
       -L "${region}" \
       --sample-name-map "${OUTPUT_DIR}/sample_map.txt" \
       --reader-threads 3

	output_vcf="${OUTPUT_DIR}/vcf/Atub_drought_commongarden_merged_${region_safe}.vcf.gz"
    done_file="${OUTPUT_DIR}/vcf/completed/Atub_drought_commongarden_merged_${region_safe}.done"
	
	#Make single vcf file for region       
	if apptainer exec /project/kreiner/espolston/gatk_4.6.2.0.sif gatk GenotypeGVCFs \
   		-R "${REFERENCE}" \
   		-V "gendb://${OUTPUT_DIR}/db/${region_safe}" \
   		-O "${output_vcf}"; then
        # Mark as complete only if successful
        echo "Completed at $(date)" > "${done_file}"
        echo "Region: ${region}" >> "${done_file}"
    fi
}
export -f run_gatk_vcf
export REFERENCE BAM_DIR INPUT_DIR OUTPUT_DIR REGIONS_FILE SAMPLES_FILE

# one srun per job; -c matches --cpus-per-task
SRUN="srun -N1 -n1 -c ${SLURM_CPUS_PER_TASK}"

# parallel caps total concurrency to your allocation size
# $SLURM_NTASKS equals nodes * ntasks-per-node
parallel --delay 0.2 \
-j "${SLURM_NTASKS}" \
--joblog "runtask-${SLURM_JOBID}.log" \
--resume \
--env REFERENCE,BAM_DIR,INPUT_DIR,OUTPUT_DIR,run_gatk_vcf,SAMPLES_FILE \
"${SRUN} bash -lc 'run_gatk_vcf {}'" \
:::: "${REGIONS_FILE}"

#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=8
#SBATCH --time=0:10:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --job-name indiv_gatk_parallel
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=espolston@rcc.uchicago.edu 

module load python/anaconda-2023.09
source activate /project/kreiner/espolston/
module load gcc
module load parallel
module load samtools
module load apptainer

# ============================================
# CONFIGURATION - Edit these paths
# ============================================
REFERENCE="/scratch/midway3/espolston/Atub_193_hap2.fasta"
BAM_DIR="/scratch/midway2/espolston/bams_commongardenanddrought"
REGIONS_FILE="/scratch/midway3/espolston/bamsdone_samplenames.txt"
OUTPUT_DIR="/scratch/midway2/espolston/gvcfs"
# ============================================

# Wrapper function to mark completion

#Region file bams as 1, sample name as 2 (bc of drought bam rename) 
run_gatk() {
    indiv=$1
    outputname=$2
    
    output_vcf="${OUTPUT_DIR}/${outputname}.g.vcf.gz"
    done_file="${OUTPUT_DIR}/completed/${outputname}.done"
	
    export TMPDIR=/scratch/midway2/espolston/tmp 
    mkdir -p ${TMPDIR}    

    # Run gatk and rename the drought samples
    if apptainer exec --bind /scratch:/scratch /project/kreiner/espolston/gatk_4.6.2.0.sif gatk --java-options "-Xmx60g -Djava.io.tmpdir=${TMPDIR}" HaplotypeCaller \
		-R "${REFERENCE}" \
    	-I "${BAM_DIR}/${indiv}" \
    	-ERC GVCF \
   	-O "${output_vcf}" \
	--native-pair-hmm-threads 8; then
        # Mark as complete only if successful
        echo "Completed at $(date)" > "${done_file}"
    fi
}

export -f run_gatk
export REFERENCE BAM_DIR OUTPUT_DIR REGIONS_FILE

# one srun per job; -c matches --cpus-per-task
SRUN="srun -N1 -n1 -c ${SLURM_CPUS_PER_TASK}"

# parallel caps total concurrency to your allocation size
# $SLURM_NTASKS equals nodes * ntasks-per-node
parallel --delay 0.2 \
--colsep '\t' \
-j "${SLURM_NTASKS}" \
--joblog "runtask-${SLURM_JOBID}.log" \
--resume \
--env REFERENCE,BAM_DIR,OUTPUT_DIR,REGIONS_FILE,run_gatk \
"${SRUN} bash -lc 'run_gatk {1} {2}'" \
:::: "${REGIONS_FILE}"

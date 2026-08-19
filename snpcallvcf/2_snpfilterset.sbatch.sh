#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --job-name filter_parallel

module load python/anaconda-2023.09
source activate /project/kreiner/espolston/
module load gcc
module load parallel

# ============================================
# CONFIGURATION - Edit these paths
# ============================================
REGIONS_FILE="/scratch/midway3/espolston/region_list/scaf_1.txt"
INPUT_DIR="/scratch/midway2/espolston/scaf1_unfilt/"
OUTPUT_DIR="/scratch/midway3/espolston/scaf1/"
# ============================================

# Create directories
mkdir -p "${OUTPUT_DIR}"

#copy all files from scaf_1_region10.txt and scaf_1_region100.txt
while read line; do
  cp /cds3/kreiner/"$line" $INPUT_DIR
done < scaf1_unfilt.txt

#writing script to run SNPFiltering.sh (editted from Julia) on each scaffolds regions for variant and invariant sites

#run snp filtering wrapper
run_filter() {
    input_vcf=$1
    prefix=`echo $input_vcf | sed 's/.vcf//g'`

    # Run invariant filter
    bcftools view -O z --include 'N_ALT = 0' $INPUT_DIR/${prefix}.vcf | \
    bcftools filter -i 'F_MISSING <= 0.25' -Oz >  "$OUTPUT_DIR/${prefix}_filt_invariant.vcf.gz" && sleep 5 && tabix $OUTPUT_DIR/${prefix}_filt_invariant.vcf.gz
    
    #run variant filter
    bcftools filter -i 'F_MISSING <= 0.25' $INPUT_DIR/${prefix}.vcf | \
                bcftools filter -i 'QUAL >= 30' | \
                bcftools filter -i 'AB >= 0.25 & AB <= 0.75 | AB <= 0.01' | \
                bcftools filter -i 'SAF > 0 & SAR > 0' | \
                bcftools filter -i 'MQM >=30 & MQMR >= 30' | \
                bcftools filter -i '((PAIRED > 0.05) & (PAIREDR > 0.05) & (PAIREDR / PAIRED < 1.75 ) & (PAIREDR / PAIRED > 0.25)) | ((PAIRED < 0.05) & (PAIREDR < 0.05))' |
                bcftools filter -o "$OUTPUT_DIR/${prefix}_filtered_snps.vcf.gz" -O z -i '((AF > 0) & (AF < 1))' && sleep 5 && tabix "OUTPUT_DIR/${prefix}_filtered_snps.vcf.gz"
    fi
}

export -f run_filter
export INPUT_DIR OUTPUT_DIR REGIONS_FILE

# one srun per job; -c matches --cpus-per-task
SRUN="srun --exclusive -N1 -n1 -c ${SLURM_CPUS_PER_TASK}"

# parallel caps total concurrency to your allocation size
# $SLURM_NTASKS equals nodes * ntasks-per-node
parallel --delay 0.2 \
-j "${SLURM_NTASKS}" \
--joblog "runtask-${SLURM_JOBID}.log" \
--resume \
--env INPUT_DIR, OUTPUT_DIR, REGIONS_FILE, run_filter \
"${SRUN} bash -lc 'run_filter {1}'" \
:::: "${REGIONS_FILE}"

#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --time=0:10:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --job-name=whatshap
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=espolston@rcc.uchicago.edu 

#Whatshap for shapeit4
#whatshap v2.8
#https://whatshap.readthedocs.io/en/latest/guide.html and https://app.paperpile.com/view/?id=a31a47c0-d4f2-4e23-abd6-4a792987fdfd
#module load python/anaconda-2023.09
#source create -n whatshap-env whatshap
#source activate whatshap-env
#conda install -c conda-forge bcftools

#for making env on randi:
#module load gcc/12.1.0
#module load python/3.10.5
#module load miniconda3/24.9.2
#conda create -n whatshap-env whatshap -c bioconda -c conda-forge
#conda init
#conda activate whatshap-env 
#conda install -c conda-forge bcftools

#samples in vcf need to make sample name in bam file (for cg:# for d:AT#_#)
#rehead_forwhatshap.txt: col 1=current_samplename col2=new_samplename
#module load htslib
#module load python/anaconda-2023.09
#source activate /project/kreiner/espolston/
##bcftools reheader --samples /scratch/midway3/espolston/vcf/rehead_forwhatshap.txt -o /project/kreiner/pairedenv_commongarden/normalized_SNPsonly_vcfs/drought_commongarden_whatshap.vcf.gz /project/kreiner/pairedenv_commongarden/normalized_SNPsonly_vcfs/drought_commongarden_merged_final_ID_IBD.vcf.gz

#tabix -p vcf /project/kreiner/pairedenv_commongarden/normalized_SNPsonly_vcfs/drought_commongarden_whatshap.vcf.gz

#from this: https://github.com/whatshap/whatshap/issues/478 add --sample=NAME and split vcf by sample and scaf
#only have to do scaf 1-16
#bcftools view --samples AT1_10 --regions Scaffold_1 /project/kreiner/pairedenv_commongarden/normalized_SNPsonly_vcfs/drought_commongarden_whatshap.vcf.gz -Oz -o /scratch/midway3/espolston/sample.vcf.gz

#tabix -p vcf /scratch/midway3/espolston/sample.vcf.gz

#whatshap phase --reference=/scratch/midway3/espolston/Atub_193_hap2.fasta -o /scratch/midway3/espolston/dcgm_whatshapout_test.vcf --sample=AT1_10 --chromosome=Scaffold_1 /scratch/midway3/espolston/sample.vcf.gz /scratch/midway2/espolston/bams_commongardenanddrought/AT1_10_193_2.dd.bam


#HERE-- make this a loop with make vcf for that sample, run only that sample, delete original 1 sample vcf
#then have a checkin like gvcf steps
#Make file of all samples and all scaffold combos from dcgm_bams.txt
#awk '{OFS="\t"; split($1, part, "_193_2.dd.bam"); $2=part[1]; print $2,$1}' /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm_bams.txt > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm_whatshap_samples.txt
#scafs.txt is just Scaffold_1 \n Scaffold_2 ... for 1-16
#awk '{OFS="\t"} NR == FNR {
#        lines[NR] = $0;
#        count = NR;
#        next
#    }
#    {
#        for (i = 1; i <= count; i++) {
#            print lines[i], $0
#        }
#    }
#' /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm_whatshap_samples.txt /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/scafs.txt > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm_whatshap_runs.txt

#for midway
module load python/anaconda-2023.09
source activate whatshap-env
module load parallel
module load apptainer

#for randi
conda activate whatshap-env


# ============================================
# CONFIGURATION - Edit these paths
# ============================================
REFERENCE="/scratch/midway3/espolston/Atub_193_hap2.fasta"
SAMPLES_FILE="/scratch/midway3/espolston/dcgm_whatshap_runs.txt"
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
    
    if bgzip -c "${OUT_DIR_START}_${region}/whatshapout_${samp}_${region}.vcf" > "${OUT_DIR_START}_${region}/whatshapout_${samp}_${region}.vcf.gz"; then
    	rm "${OUT_DIR_START}_${region}/whatshapout_${samp}_${region}.vcf"
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
	

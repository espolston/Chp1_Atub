#!/bin/bash
# ============================================
# CONFIGURATION - Edit these paths
# ============================================
#Regions file for 1st run
SAMPLES_FILE="/scratch/midway3/espolston/dcgm_whatshap_runs.txt"
OUT_DIR_START="/scratch/midway2/espolston/whatshap"
FAILED_REGIONS="/scratch/midway3/espolston/failed_regions_whatshap.txt"
# ============================================

> "$FAILED_REGIONS"  # Clear the file

#separates out each line
while IFS=$'\t' read -r samp bam region; do
    #check for done file
    if [ ! -f "${OUT_DIR_START}/${samp}_${region}.done" ]; then
	echo "${samp}"$'\t'"${bam}"$'\t'"${region}" >> "$FAILED_REGIONS"  # Write original region format
    fi
done < "$REGIONS_FILE"

num_failed=$(wc -l < "$FAILED_REGIONS")
num_total=$(wc -l < "$SAMPLES_FILE")
num_completed=$((num_total - num_failed))

echo "Summary:"
echo "  Total whatshap runs: $num_total"
echo "  Completed: $num_completed"
echo "  Failed/Not Run: $num_failed"
echo
echo "Failed regions written to: $FAILED_REGIONS"
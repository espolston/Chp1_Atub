#!/bin/bash
# ============================================
# CONFIGURATION - Edit these paths
# ============================================
#Regions file for 1st run
REGIONS_FILE="/scratch/midway3/espolston/bamsdone_samplenames.txt"
#regions file for 2nd run bc need to move some scaffolds off of scratch
#REGIONS_FILE="/scratch/midway3/espolston/failed_regions_gvcf.txt"
OUTPUT_DIR="/scratch/midway2/espolston/"
FAILED_REGIONS="/scratch/midway3/espolston/failed_regions_gvcf.txt"
# ============================================

> "$FAILED_REGIONS"  # Clear the file

#separates out each line
while IFS=$'\t' read -r file outputname; do
    #check for done file
    if [ ! -f "${OUTPUT_DIR}/gvcfs/completed/${outputname}.done" ]; then
	echo "${file}"$'\t'"${outputname}" >> "$FAILED_REGIONS"  # Write original region format
    fi
done < "$REGIONS_FILE"

num_failed=$(wc -l < "$FAILED_REGIONS")
num_total=$(wc -l < "$REGIONS_FILE")
num_completed=$((num_total - num_failed))

echo "Summary:"
echo "  Total gvcfs: $num_total"
echo "  Completed: $num_completed"
echo "  Failed/Not Run: $num_failed"
echo
echo "Failed regions written to: $FAILED_REGIONS"
#!/bin/bash
# ============================================
# CONFIGURATION - Edit these paths
# ============================================
#regions file for first run
REGIONS_FILE="/scratch/midway3/espolston/100kbregions"
#regions file for all other runs (bc not all files in scratch)
#REGIONS_FILE="/scratch/midway3/espolston/failed_regions_vcf.txt"
OUTPUT_DIR="/scratch/midway3/espolston"
FAILED_REGIONS="/scratch/midway3/espolston/failed_regions_vcf.txt"
# ============================================

> "$FAILED_REGIONS"  # Clear the file

#separates out each line
while IFS= read -r region; do

    # Sanitize region name same way as main script
    region_safe=$(echo "$region" | tr ':' '_' | tr '-' '_')
    
    #check for done file
    if [ ! -f "${OUTPUT_DIR}/vcf/completed/Atub_drought_commongarden_merged_${region_safe}.done" ]; then
        echo "$region" >> "$FAILED_REGIONS"  # Write original region format
    fi
done < "$REGIONS_FILE"

num_failed=$(wc -l < "$FAILED_REGIONS")
num_total=$(wc -l < "$REGIONS_FILE")
num_completed=$((num_total - num_failed))

echo "Summary:"
echo "  Total regions: $num_total"
echo "  Completed: $num_completed"
echo "  Failed/Not Run: $num_failed"
echo
echo "Failed regions written to: $FAILED_REGIONS"

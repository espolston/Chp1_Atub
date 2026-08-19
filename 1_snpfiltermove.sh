#!/bin/bash
#Run this to copy all scaffold 1 files into scratch, then change headings in 
mkdir /scratch/midway3/espolston/scaf1_unfilt/

#copy all files from scaf_1_region10.txt and scaf_1_region100.txt
while read line; do
  cp "/cds3/kreiner/commongarden_drought_mergedsnp/region_vcfs/$line" /scratch/midway2/espolston/scaf1_unfilt/
done < scaf_1_region100.txt

while read line; do
  cp "/cds3/kreiner/commongarden_drought_mergedsnp/10kb_region_vcfs/$line" /scratch/midway2/espolston/scaf1_unfilt/
done < scaf_1_region10.txt
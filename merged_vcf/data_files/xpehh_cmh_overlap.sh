#!/bin/bash
#input files
#windows <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm_10kbsites.bed")
#xpehh <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/sel_scan/Aggregated_xpehh.txt", header = T)
#cmh_loci <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/P.05FDR_dcgm")

#make cmh and xpehh bed files
#filter to only the sites that are above the critical XPEHH threshold
#10 is norm_xpehh, 13 is the new pos 2, so end up w chrom pos1 pos2 xpehh
#awk 'BEGIN{OFS="\t"} NR > 1 && $10 > 5 {$13 = $2+1; print $3,$2,$13,$10}' /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/sel_scan/Aggregated_xpehh.txt > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/sel_scan/xpehh_up.bed
#awk 'BEGIN{OFS="\t"} NR > 1 && $10 < -5 {$13 = $2+1; print $3,$2,$13,$10}' /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/sel_scan/Aggregated_xpehh.txt > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/sel_scan/xpehh_down.bed
#cat /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/sel_scan/xpehh_up.bed /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/sel_scan/xpehh_down.bed > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/xpehhhits_dcgm.bed

#filter on 14 for bonferroni p < .05, 15 is the new pos 2, so end up w chrom pos1 pos2 bon_pval
#awk 'BEGIN{OFS="\t"} NR > 1 && $14 < .05 {$15 = $3+1; print $1,$3,$15,$14}' /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/P.05FDR_dcgm > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/cmhhits_dcgm.bed

#now see how many of each are in the 10kb windows
for sc in {1..16}; do 
#filter by scaf and then remove tab at end of line
	awk -v scaf=${sc} 'BEGIN{OFS="\t"; target="Scaffold_" scaf} $1 == target {sub(/\t+$/, ""); print $1,$2,$3,$4}' /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/cmhhits_dcgm.bed > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/cmh_${sc}.bed
	awk -v scaf=${sc} 'BEGIN{OFS="\t"; target= scaf} $1 == target {$5 = "Scaffold_" scaf; sub(/\t+$/, ""); print $5,$2,$3,$4}' /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/xpehhhits_dcgm.bed > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/xpehh_${sc}.bed
	awk -v scaf=${sc} 'BEGIN{OFS="\t"; target="Scaffold_" scaf} $1 == target {sub(/\t+$/, ""); print $1,$2,$3}' /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm_10kbsites.bed > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/windows_${sc}.bed
	
	#sort by pos
	sort -k1,1 -k2,2n /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/cmh_${sc}.bed > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/cmh_${sc}_sort.bed
	sort -k1,1 -k2,2n /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/xpehh_${sc}.bed > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/xpehh_${sc}_sort.bed
	sort -k1,1 -k2,2n /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/windows_${sc}.bed > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/windows_${sc}_sort.bed
	
	#Get how many hits of cmh and xpehh are within a 10 kb window
	bedtools map -a /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/windows_${sc}_sort.bed -b /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/xpehh_${sc}_sort.bed -c 2 -o count > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_xpehh_${sc}.bed

	bedtools map -a /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/windows_${sc}_sort.bed -b /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/cmh_${sc}_sort.bed -c 2 -o count > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_cmh_${sc}.bed
	
	rm /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/cmh_${sc}* /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/xpehh_${sc}* /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/windows_${sc}*; done

#merge all of these files into one
cat /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_cmh_* > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_cmh_all.txt
cat /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_xpehh_* > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_xpehh_all.txt


#check which of the cmh windows are the clumped sites
#clump_loci <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clumphits_dcgm.txt")
awk 'BEGIN{OFS="\t"} NR > 1 {$6 = $2+1; print $1,$2,$6,$3}' /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clumphits_dcgm.txt > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/clumpcmhhits_dcgm.bed

for sc in {1..16}; do 
#filter by scaf and then remove tab at end of line
	awk -v scaf=${sc} 'BEGIN{OFS="\t"; target=scaf} $1 == target {sub(/\t+$/, ""); $5 = "Scaffold_" scaf; print $5,$2,$3,$4}' /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/clumpcmhhits_dcgm.bed > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/cmh_${sc}.bed
	awk -v scaf=${sc} 'BEGIN{OFS="\t"; target="Scaffold_" scaf} $1 == target {sub(/\t+$/, ""); print $1,$2,$3}' /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm_10kbsites.bed > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/windows_${sc}.bed
	
	#sort by pos
	sort -k1,1 -k2,2n /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/cmh_${sc}.bed > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/cmh_${sc}_sort.bed
	sort -k1,1 -k2,2n /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/windows_${sc}.bed > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/windows_${sc}_sort.bed
	
	#Get how many hits of cmh and xpehh are within a 10 kb window
	bedtools map -a /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/windows_${sc}_sort.bed -b /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/cmh_${sc}_sort.bed -c 2 -o count > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_clumpcmh_${sc}.bed
	
	rm /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/cmh_${sc}* /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/windows_${sc}*; done

cat /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_clumpcmh_* > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_clumpcmh_all.txt


#Run clumping on xpehh on cluster
#raw p cutoffs from chp1_mergedvcf.Rmd (after nonclumped manhattan plots)
#--clump-p1 from xpehh cutoff of + or - 5 (took max pvalue from zscore) --clump-p2 from xpehh cutoff of + or - 2
module load plink
plink --bfile /scratch/midway3/espolston/binary_dcgm --fam /scratch/midway3/espolston/dcgm.fam --clump /scratch/midway3/espolston/xpehh_plink --clump-p1 .0000005732734 --clump-kb 1000 --clump-r2 0.04550026 --clump-field pval --allow-no-sex --allow-extra-chr --clump-snp-field locus_id --make-founders --out /scratch/midway3/espolston/clumped_dcgm_xpehh

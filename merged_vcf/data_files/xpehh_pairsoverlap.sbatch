#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --time=2-00:00:00
#SBATCH --partition=tier1q
#SBATCH --job-name=pairoverlap

set -euo pipefail

#i did have claude write this script

#=======================================================================
# Purpose: For each paired-population XPEHH scan (pair 1-17), count how
# many sites with |norm_xpehh| > 5 fall into each 10kb window (chrom 1-16),
# using bedtools map, mirroring the logic in xpehh_cmh_overlap.sh but
# looping over the 17 Aggregated_pair_*_xpehh.txt files instead of a
# single aggregated/pooled file.
#=======================================================================

# ---- EDIT THESE PATHS AS NEEDED ----
INDIR=/scratch/espolston/sel_scan_pair          # where Aggregated_pair_*_xpehh.txt live
WINDOWS_BED=/scratch/espolston/dcgm_10kbsites.bed
OUTDIR=/scratch/espolston/sel_scan_pair/window_counts
TMPDIR=${OUTDIR}/tmp
mkdir -p "${OUTDIR}" "${TMPDIR}"

# column indices in Aggregated_pair_*_xpehh.txt (space-delimited, header present):
# 1 pair  2 locus_id  3 POS  4 CHROM  5 gpos  6 popA_freq  7 ihhA
# 8 popB_freq  9 ihhB  10 raw_xpehh  11 norm_xpehh  12 crit  13 P
POS_COL=3
CHROM_COL=4
NORMXPEHH_COL=11

module load gcc/12.1.0
module load bedtools/2.30.0

for p in {1..17}; do
    infile="${INDIR}/Aggregated_pair_${p}_xpehh.txt"
    if [[ ! -f "${infile}" ]]; then
        echo "WARNING: ${infile} not found, skipping pair ${p}" >&2
        continue
    fi

    # build a bed of sites with |norm_xpehh| > 5 for this pair:
    # chrom(prefixed), start, end, norm_xpehh
    sig_bed="${TMPDIR}/pair${p}_sigxpehh.bed"
    awk -v pos_col="${POS_COL}" -v chrom_col="${CHROM_COL}" -v xp_col="${NORMXPEHH_COL}" \
        'BEGIN{OFS="\t"} NR > 1 {
            xp = $xp_col
            if (xp > 5 || xp < -5) {
                print "Scaffold_" $chrom_col, $pos_col, $pos_col+1, xp
            }
        }' "${infile}" > "${sig_bed}"

    for sc in {1..16}; do
        target="Scaffold_${sc}"

        sig_sc="${TMPDIR}/pair${p}_sigxpehh_${sc}.bed"
        win_sc="${TMPDIR}/windows_${sc}.bed"

        awk -v target="${target}" 'BEGIN{OFS="\t"} $1 == target {print $1,$2,$3,$4}' \
            "${sig_bed}" | sort -k1,1 -k2,2n > "${sig_sc}"

        # windows file only needs to be built once per chrom, reuse if present
        if [[ ! -s "${win_sc}" ]]; then
            awk -v target="${target}" 'BEGIN{OFS="\t"} $1 == target {print $1,$2,$3}' \
                "${WINDOWS_BED}" | sort -k1,1 -k2,2n > "${win_sc}"
        fi

        bedtools map -a "${win_sc}" -b "${sig_sc}" -c 2 -o count \
            > "${OUTDIR}/window_pair${p}_xpehh_${sc}.bed"
    done
done

# merge everything into one file, tagging each row with its pair number
: > "${OUTDIR}/window_xpehh_allpairs.txt"
for p in {1..17}; do
    for sc in {1..16}; do
        f="${OUTDIR}/window_pair${p}_xpehh_${sc}.bed"
        [[ -f "${f}" ]] || continue
        awk -v pair="${p}" 'BEGIN{OFS="\t"} {print pair, $0}' "${f}" \
            >> "${OUTDIR}/window_xpehh_allpairs.txt"
    done
done

# header for the merged file: pair, chrom, win_start, win_end, n_sig_sites
sed -i '1i pair\tchrom\twin_start\twin_end\tn_sig_sites' "${OUTDIR}/window_xpehh_allpairs.txt"

# cleanup intermediates
rm -rf "${TMPDIR}"

echo "Done. Per-pair/per-chrom window counts in ${OUTDIR}/window_pair<P>_xpehh_<SC>.bed"
echo "Combined table: ${OUTDIR}/window_xpehh_allpairs.txt"

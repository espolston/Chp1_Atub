#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

#Rscript h12_recombcorrec.R recomb_file input_file_start end
#ex: Rscript h12_recombcorrec.R recomb_ratebins.txt /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12_H2H1/dcgm_hap_stats_Scaffold_ .txt (or ag.txt...)

if(length(args) == 0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args) == 3){
	recomb_file <- args[1]
	infile_start <- args[2]
	end <- args[3]
} else{
	stop("Check arguments")
}

library(data.table)
library(tidyverse)

#first- need to filter H12 based on recombination rate bins - this almost the same as Xpehh recomb rate check
#map <- fread("dcgm_cM.map")
#colnames(map) <- c("Scaffold", "SNP", "cM", "Locus")
#map$chr <- as.numeric(gsub("Scaffold_([0-9]+).*", "\\1", map$Scaffold))

#covert cM to r (cM in rhotocM_forlibby_new.R is just r*100)
#map$r <- map$cM/100

#make recomb rate bins
#bins <- quantile(map$r)

#assign recomb rates to these bins and set window positions
#rec_bins <- map %>%
#  group_by(chr) %>%
#  mutate("win_pos_1" = case_when(
#    Locus == min(Locus) ~ 0, 
#    Locus > min(Locus) ~ lag(Locus))) %>%
#  mutate("win_pos_2" = Locus) %>%
#  select(chr, SNP, win_pos_1, win_pos_2, r) %>%
#  mutate("bin" = case_when(
#    r <= bins[2] ~ "1",
#    r > bins[2] & r <= bins[3] ~ "2",
#    r > bins[3] & r <= bins[4] ~ "3",
#    r > bins[4] & r <= bins[5] ~ "4"
#  )) %>%
#  ungroup()

#write.table(rec_bins, "recomb_ratebins.txt", row.names = F, quote = F)

#load in ag and nat combined run
for(s in 1:16){
	H12Scan_all = read.table(paste(infile_start, s, end, sep = ""))
  colnames(H12Scan_all) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2/H1", "H123")
  H12Scan_all$chr <- rep(s, length(H12Scan_all$win_center))

  #get rec rate for each of our H12 scores
  rec_bins <- read.table(recomb_file, header = T)
  colnames(rec_bins)
  
  rec_bins_wr <- H12Scan_all %>%
    left_join(
    rec_bins %>% select(chr, win_pos_1, win_pos_2, r, bin),
    join_by(chr == chr, win_center >= win_pos_1, win_center <= win_pos_2)) %>%
    rename(rec_rate = r) %>%
    distinct(win_center, .keep_all = TRUE) %>%
    filter(bin > 1) %>%
    mutate(win_pos_1 = win_pos_1.x, win_pos_2 = win_pos_2.x) %>%
    select(win_center, win_pos_1, win_pos_2, num_unique_hap, hap_freq_spectrum, num_hap_in_each_freqbin, H1, H2, H12, 'H2/H1', H123)
  
  write.table(rec_bins_wr, paste(infile_start, s, "_recombfiltered", end, sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
}

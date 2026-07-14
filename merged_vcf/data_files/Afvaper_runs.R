#AFvapeR to look at amount of parallelism - write here and then run on server in R/4.6.0 
#from: https://github.com/JimWhiting91/afvaper
#to install on randi:
  #mkdir -p ~/Rtmp
  #export TMPDIR=~/Rtmp
  #R
#inR: 
  # install.packages(c("data.table", "memuse", "vcfR"))
  # install.packages("remotes",repos = "http://cran.us.r-project.org")
  # remotes::install_github("JimWhiting91/afvaper")
  # install.packages(c("ragg", "systemfonts", "textshaping"))

#to get installed on cluster (midway)
#on command line: git clone https://github.com/JimWhiting91/afvaper.git 
#then in R: remotes::install_local("/scratch/midway3/espolston/afvaper")
library(afvaper,verbose = F)
library(vcfR,verbose = F)
library(data.table)
library(dplyr)

#making popmap
#awk 'NR > 1 {print $1,$3,$5}' /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/Chp1_DroughtandCommonGardenInfo_Merged-MetaData.tsv > /Users/libbypolston/Desktop/popmap.txt

#import fam file for sample name, pair, contrast
popmap <- read.table("/scratch/espolston/dcgm.fam")

#maps sample name in vcf to pair_environment
popmap <- popmap %>%
  mutate("indiv" = case_when(
    grepl("_", V1) == T ~ paste(V1, "_", V2, "_T", sep =""),
    grepl("_", V1) == F ~ V1), 
    "env" = case_when(
      V6 == 1 ~ "NAT",
      V6 == 2 ~ "AG"), "contrast" = paste(V3, "_", env, sep ="")) %>%
  select(indiv, contrast)

#just for this one bc it has AT - comment to popmap <- pop_convert when done
#maps sample name in vcf (when with AT) to pair_environment
#popmap <- popmap %>%
#  mutate("indiv" = case_when(
#    grepl("_", V1) == T ~ paste(V1, "_", V2, "_T", sep =""),
#    grepl("_", V1) == F ~ V1), 
#    "env" = case_when(
#      V6 == 1 ~ "NAT",
#      V6 == 2 ~ "AG"), "contrast" = paste(V3, "_", env, sep ="")) %>%
#  select(indiv, contrast)
#convert_AT <- read.table("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/bam_samplenames.txt")
#pop_convert <- merge(popmap, convert_AT, by.x = "indiv", by.y = "V2") %>% select(V1,contrast)
#popmap <- pop_convert

#individuals in column 1 and population/habitat/ecotype in column 2 -> i think this would need to be pair1_ag etc
vector_list <- list(c("17_AG","17_NAT"),
                    c("16_AG","16_NAT"),
                    c("15_AG","15_NAT"),
                    c("14_AG","14_NAT"),
                    c("13_AG","13_NAT"),
                    c("12_AG","12_NAT"),
                    c("11_AG","11_NAT"),
                    c("10_AG","10_NAT"),
                    c("9_AG","9_NAT"),
                    c("8_AG","8_NAT"),
                    c("7_AG","7_NAT"),
                    c("6_AG","6_NAT"),
                    c("5_AG","5_NAT"),
                    c("4_AG","4_NAT"),
                    c("3_AG","3_NAT"),
                    c("2_AG","2_NAT"),
                    c("1_AG","1_NAT"))
# Name vectors for contrasts
names(vector_list) <- c("pair_17", "pair_16","pair_15","pair_14","pair_13","pair_12","pair_11","pair_10","pair_9", "pair_8","pair_7","pair_6","pair_5","pair_4","pair_3","pair_2","pair_1")

#------ Calculate allele frequency change vectors -------

# Set our window size - "physical genome length at which linkage decays below a certain value, and matching that to the median physical size  of AF-vapeR windows."
# Did LD decay calc, got median window size and then calculated avg # snps in that window size genomewide
window_snps = 200

# How many permutations do we want in total?
total_perms <- 10000

#index with chromosome sizes
#got chrom lengths with: bcftools view -h /project/kreiner/pairedenv_commongarden/normalized_SNPsonly_vcfs/drought_commongarden_merged_final_ID_IBD.vcf.gz > header.txt
genome_fai <- data.frame(chr=c("Scaffold_1","Scaffold_2","Scaffold_3", "Scaffold_4","Scaffold_5","Scaffold_6","Scaffold_7","Scaffold_8","Scaffold_9","Scaffold_10","Scaffold_11","Scaffold_12","Scaffold_13","Scaffold_14","Scaffold_15","Scaffold_16"),
                         size=c(59452825,52828123,46342860,36395640,37853689,36916374,32560229,34047692,28498443,30665661,36917966,32568677,28832420,25823852,26125370,27492738)) 

# Fetch proportional size of all chromosomes
chr_props <- genome_fai$size/sum(genome_fai$size)
chr_perms <- data.frame(chr=genome_fai$chr,
                        perms=round(chr_props * total_perms))

# This gives us approximately 10000 null perms in total, distributed across the genome according to relative size of chromosomes...

all_chr_res <- lapply(1:16,function(i){
  vcf_in <- read.vcfR(paste("/scratch/espolston/new_phased/phased_dcgm_Scaffold_", i, "_newhead.vcf.gz", sep = ""), verbose = F)
  
  chr_AF_input <- calc_AF_vectors(vcf = vcf_in,
                                  window_size = window_snps,
                                  popmap = popmap,
                                  vectors = vector_list,
                                  n_cores = 4)
  
  chr_null_input <- calc_AF_vectors(vcf = vcf_in,
                                    window_size = window_snps,
                                    popmap = popmap,
                                    vectors = vector_list,
                                    n_cores = 4,
                                    null_perms = chr_perms$perms[i],
                                    data_type = "vcf")
  cat(paste("scaffold ", i, "done\n"), file = "/scratch/espolston/afvaper_output.txt", append = TRUE)
  
  return(list(chr_AF_input,chr_null_input))
})

saveRDS(all_chr_res, "/scratch/espolston/afvapr/all_chr_res.txt")
#all_chr_res <- readRDS("/scratch/espolston/afvapr/all_chr_res.txt")
  cat("merged res file done\n", file = "/scratch/espolston/afvaper_output.txt", append = TRUE)
  
  AF_input <- merge_eigen_res(lapply(all_chr_res,'[[',1))
  null_input <- merge_eigen_res(lapply(all_chr_res,'[[',2))
  
  #------ Perform eigen analysis over allele frequency matrices  -------
  
  ###eigen analysis over AF matrices
  # Perform eigen analysis
  #The eigen_res output is a list containing, for each matrix, the eigenvalue distribution, the eigenvector loadings, and the projected A matrix that shows per-SNP scores for each eigenvector
  eigen_res <- lapply(AF_input,eigen_analyse_vectors)
  # View chromosomal regions:
  #head(names(eigen_res))
  
  saveRDS(eigen_res, "/scratch/espolston/afvapr/Afvaper_eigenresiduals.txt")
  #eigen_res <- readRDS("/scratch/espolston/afvapr/Afvaper_eigenresiduals.txt")
  cat("eig file done\n", file = "/scratch/espolston/afvaper_output.txt", append = TRUE)
  
  # View eigenvalue distribution of first matrix
  #eigen_res[[1]]$eigenvals
  # View eigenvector loadings of first matrix
  #eigen_res[[1]]$eigenvecs
  # View head of SNP scores
  #head(eigen_res[[1]]$A_matrix)
  
  #------ Calculate null cutoffs -------
  # Get cutoffs for 95%, 99% and 99.9%
  null_cutoffs <- find_null_cutoff(null_input,cutoffs = c(0.90,0.95,0.99))
  null_cutoffs
  
  #------ Calculate empirical pvalues -------
  # Calculate p-vals
  pvals <- eigen_pvals(eigen_res,null_input)
  
  # Showpvals
  #head(pvals)
  saveRDS(pvals, "/scratch/espolston/afvapr/Afvaper_pvals.txt")
  cat("p file file done\n", file = "/scratch/espolston/afvaper_output.txt", append = TRUE)
  
  #------ Plot eigenvalues along chromosome -------
  library(ggplot2)
  # Plot the raw eigenvalues, and visualise the cutoff of 90%
  all_plots <- eigenval_plot(eigen_res,cutoffs = null_cutoffs[,"90%"])
  
  # Show the plots for eigenvalue 1
  eigenval_all <- all_plots[[1]] + ggtitle("Raw eigenvalues for Eig 1 w 90% cutoff")
  ggsave("/scratch/espolston/afvapr/afvapr_eigenvalfig.png", eigenval_all, device = "png", dpi = 300, height = 6, width = 10, units = "in")
  #problems with randi ggsave version
  #ggsave("/scratch/espolston/afvapr/afvapr_eigenvalfig.png",
  #       eigenval_all, dpi = 300, height = 6, width = 10, units = "in",
  #       device = function(filename, width, height, ...) {
  #         grDevices::png(filename, width = width, height = height, res = 300, type = "cairo", ...)
  #       })
  
  # Plot empirical p-values, -log10(p) of 2 ~ p=0.01, 3 ~ p=0.001 etc.
  all_plots_p <- eigenval_plot(eigen_res,null_vectors = null_input,plot.pvalues = T)
  
  # Show the plots for eigenvalue 1
  allchr_plot <- all_plots_p[[1]] + ggtitle("Empirical pvalues for Eig 1") 
  ggsave(allchr_plot, file= "/scratch/espolston/afvapr/afvapr_pvalfig.png", device = "png", dpi = 300, height = 6, width = 10, units = "in")
  
  #problems with randi ggsave
  #ggsave(allchr_plot, file= "/scratch/espolston/afvapr/afvapr_pvalfig.png", dpi = 300, height = 6, width = 10, units = "in", device = function(filename, width, height, ...) {
  #         grDevices::png(filename, width = width, height = height, res = 300, type = "cairo", ...)
  #       })
  
  #------ Pull Significant windows -------
  
  # Find significant windows above 99.9% null permutation
  significant_windows <- signif_eigen_windows(eigen_res,null_cutoffs[,"90%"])
  
  # Display 'outliers'
  #significant_windows
  
  saveRDS(significant_windows, file ="/scratch/espolston/afvapr/afvapr_sigwindows.txt")
  #significant_windows <- readRDS("/scratch/espolston/afvapr/afvapr_sigwindows.txt")
  cat("sigwin file done\n", file = "/scratch/espolston/afvaper_output.txt", append = TRUE)
  
  
  significant_windows95 <- signif_eigen_windows(eigen_res,null_cutoffs[,"95%"])
  saveRDS(significant_windows95, file ="/scratch/espolston/afvapr/afvapr_sigwindows_95.txt")
  significant_windows99 <- signif_eigen_windows(eigen_res,null_cutoffs[,"99%"])
  saveRDS(significant_windows99, file ="/scratch/espolston/afvapr/afvapr_sigwindows_99.txt")
  
  #------ Summarize Outliers -------
  
  # Summarise parallel evolution in windows that are significant on eigenvector 1
  eig1_parallel <- summarise_window_parallelism(window_id = significant_windows[[1]],
                                                eigen_res = eigen_res,
                                                loading_cutoff = 0.3,
                                                eigenvector = 1)
  # Show results
  head(eig1_parallel)
  
  write.table(eig1_parallel, file = "/scratch/espolston/afvapr/afvapr_eig1_parallel.txt", quote = F, row.names = F)
  cat("eig file done\n", file = "/scratch/espolston/afvaper_output.txt", append = TRUE)
  
  
  #For outliers on eigenvectors 2+ we have an additional eigenvalue_sum column that describes the sum of eigenvalues 1 + 2, as well as the individual eigenvalue 1 and eigenvalue 2 scores for every window. These tell us that most of these windows are exhibiting a signature closer to full-parallelism (large eigenvalue 1) rather than multi-parallelism (more balanced eigenvalue 1 + 2), which is expected for this simulation (these regions are around the focal 10 Mb fully parallel sweep).
  eig2_parallel <- summarise_window_parallelism(window_id = significant_windows[[2]],
                                                eigen_res = eigen_res,
                                                loading_cutoff = 0.3,
                                                eigenvector = 2)
  # Show results
  head(eig2_parallel)
  
  write.table(eig2_parallel, file = "/scratch/espolston/afvapr/afvapr_eig2_parallel.txt", quote = F, row.names = F)
  cat("eig file done\n", file = "/scratch/espolston/afvaper_output.txt", append = TRUE)
  
  #------ Explore candidate regions -------
  # Fetch an A matrix
  A_mat <- eigen_res[[1]]$A_matrix
  head(A_mat)
  
  to_plot <- data.frame(snp=rownames(A_mat),
                        eig1_score=A_mat[,1])
  
  to_plot <- to_plot %>% tidyr::separate("snp",into=c("sc","scaf","pos"),sep="_")
  to_plot$pos <- as.integer(to_plot$pos)
  
  to_plot2 <- ggplot(to_plot,aes(x=pos,y=abs(eig1_score)))+
    geom_point()+
    labs(y="Eig1 Score",x="Pos (bp)") + ggtitle("Candidate regions eig 1")
  
  ggsave(to_plot2, file = "/scratch/espolston/afvapr/afvapr_candidatefig.png", device = "png", dpi = 300, height = 6, width = 10, units = "in")
  
  #problems with randi ggsave
  #ggsave(to_plot2, file = "/scratch/espolston/afvapr/afvapr_candidatefig.png", dpi = 300, height = 6, width = 10, units = "in", device = function(filename, width, height, ...) {
  #  grDevices::png(filename, width = width, height = height, res = 300, type = "cairo", ...)
  #})
  
  cat("done\n", file = "/scratch/espolston/afvaper_output.txt", append = TRUE)

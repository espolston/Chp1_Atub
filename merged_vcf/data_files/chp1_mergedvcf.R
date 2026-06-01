---
title: "Chp1_mergedvcf"
author: "Libby Polston"
date: "2025-12-19"
output: pdf_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(data.table)
library(tidyverse)
library(viridis) 
mako(6)

#colors: 
colors <- c("#0B0405FF", "#382A54FF", "#395D9CFF", "#3497A9FF", "#60CEACFF", "#DEF5E5FF")
#ag = "#60CEACFF"
#nat = "#382A54FF"
#light gray for highlighting things: "#C8C8C8" 
```

```{r}
#sex: 1=M, 2=F, 0=unknown
#enviro: 1=natural/control, 2=ag/case
dcgm_metadata <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/Chp1_DroughtandCommonGardenInfo_MergedMetaData.tsv")
dcg_metadata <- dcg_metadata %>%
  mutate("Env"= case_when(
    Environment == 1 ~ "Nat",
    Environment == 2 ~ "Ag"
  )) %>%
  mutate("Sex_code" = case_when(
    Sex == 1 ~ "M",
    Sex == 2 ~ "F",
    Sex == 0 ~ "Unknown"))
```

#MDS and PCA for pre relatedness removal (this fam, bim etc files overwritten for 440 indiv but mds and pca files not rerun)
```{r}
mds_dat <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm_mds.mds")
fam <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm.fam")
colnames(fam) <- c("Name", "Num", "Pair", "a1", "Sex", "Enviro")
colnames(mds_dat) <- c("Name", "Num", "SOL", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12", "D13", "D14", "D15", "D16", "D17", "D18", "D19", "D20")

mds_dat <- mds_dat %>%
  mutate("Pair" = fam$Pair, "en_num" = fam$Enviro) %>%
  rowwise() %>%
  mutate("Enviro"= case_when(
    en_num == 1 ~ "Nat",
    en_num == 2 ~ "Ag"
  )) %>%
  mutate("Dataset" = case_when(
    str_detect(Name, "_") == T ~ "Drought", 
    str_detect(Name, "_") == F ~ "Common Garden"
  )) %>%
    mutate("lab" = case_when(
    Dataset == "Drought" ~ as.character(paste(Name, Num, "T", sep = "_")),
    Dataset == "Common Garden" ~ as.character(paste(Name, sep = ""))))
mds_dat$Pair <- as.factor(mds_dat$Pair)
mds_dat$Enviro <- as.factor(mds_dat$Enviro)
mds_dat$Dataset <- as.factor(mds_dat$Dataset)

cgd <- read.table("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/Chp1_DroughtandCommonGardenInfo_Merged-MetaData.tsv", sep="\t",na.strings = c("","NA"),header=T)
cols <- c(1,7,8)
cgd <- cgd[,cols]
mds_winfo <- merge(mds_dat, cgd, by.x = "lab", by.y = "Sample_Name", all= T)

mds_plot <- ggplot(mds_dat, aes(x = D1, y = D2)) + geom_point(aes(color = Dataset, shape = Enviro, alpha = .4)) + ggtitle("MDS with .01 MAF .05 missingness")
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/MDS_dcgm.png", plot = mds_plot, device = "png", height = 6, width = 12, units = "in") 

#checking proportion of variance from PCA
pca_pro <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm_pca.eigenval")
pca_pro$per_var <- pca_pro$V1/sum(pca_pro$V1)

#Looking at plots now
pca_dat <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm_pca.eigenvec")
colnames(pca_dat) <- c("Name", "Num", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")

pca_dat <- pca_dat %>%
  mutate("Pair" = fam$Pair, "en_num" = fam$Enviro) %>%
  rowwise() %>%
  mutate("Enviro"= case_when(
    en_num == 1 ~ "Nat",
    en_num == 2 ~ "Ag"
  )) %>%
  mutate("Dataset" = case_when(
    str_detect(Name, "_") == T ~ "Drought", 
    str_detect(Name, "_") == F ~ "Common Garden"
  )) %>%
    mutate("lab" = case_when(
    Dataset == "Drought" ~ as.character(paste(Name, Num, "T", sep = "_")),
    Dataset == "Common Garden" ~ as.character(paste(Name, sep = ""))))
pca_dat$Pair <- as.factor(pca_dat$Pair)
pca_dat$Enviro <- as.factor(pca_dat$Enviro)
pca_dat$Dataset <- as.factor(pca_dat$Dataset)

pca_winfo <- merge(pca_dat, cgd, by.x = "lab", by.y = "Sample_Name", all= T)

#PCA for dataset
pca_plot <- ggplot(pca_winfo, aes(x = PC1, y = PC2)) + geom_point(aes(color = Dataset, shape = Enviro)) + ggtitle("PCA with .01 MAF .05 missingness")
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/PCA_dcgm.png", plot = pca_plot, device = "png", height = 6, width = 12, units = "in") 

model <- lm(PC1 ~ Dataset + Enviro + Lat + Long, pca_winfo)
summary(model)

model2 <- lm(PC2 ~ Dataset + Enviro + Lat + Long, pca_winfo)
summary(model2)

```

#PCA with related indivs removed
```{r}
library(data.table)
library(tidyverse)

#PCA w related indivs removed
fam <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm.fam")
colnames(fam) <- c("Name", "Num", "Pair", "a1", "Sex", "Enviro")

#checking proportion of variance from PCA
pca_pro <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm_pca_after_IBD.eigenval")
pca_pro$per_var <- pca_pro$V1/sum(pca_pro$V1)

#Looking at plots now
pca_dat <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm_pca_after_IBD.eigenvec")
colnames(pca_dat) <- c("Name", "Num", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")

#adding in metadata
pca_dat <- pca_dat %>%
  mutate("Pair" = fam$Pair, "en_num" = fam$Enviro) %>%
  rowwise() %>%
  mutate("Enviro"= case_when(
    en_num == 1 ~ "Nat",
    en_num == 2 ~ "Ag"
  )) %>%
  mutate("Dataset" = case_when(
    str_detect(Name, "_") == T ~ "Drought", 
    str_detect(Name, "_") == F ~ "Common Garden"
  )) %>%
  mutate("lab" = case_when(
    Dataset == "Drought" ~ as.character(paste(Name, Num, "T", sep = "_")),
    Dataset == "Common Garden" ~ as.character(paste(Name, sep = ""))))
pca_dat$Pair <- as.factor(pca_dat$Pair)
pca_dat$Enviro <- as.factor(pca_dat$Enviro)
pca_dat$Dataset <- as.factor(pca_dat$Dataset)

cgd <- read.table("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/Chp1_DroughtandCommonGardenInfo_Merged-MetaData.tsv", sep="\t",na.strings = c("","NA"),header=T)
cols <- c(1,7,8)
cgd <- cgd[,cols]

pca_winfo <- merge(pca_dat, cgd, by.x = "lab", by.y = "Sample_Name", all= T)

#PCA for dataset
library(ggrepel)
pca_plot <- ggplot(pca_winfo, aes(x = PC1, y = PC2)) + geom_point(aes(color = Pair, shape = Enviro)) + ggtitle("PCA with .01 MAF .05 missingness") 
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/PCA_dcgm.png", plot = pca_plot, device = "png", height = 6, width = 12, units = "in") 

model <- lm(PC1 ~ Dataset + Enviro + Lat + Long, pca_winfo)
summary(model)

model2 <- lm(PC2 ~ Dataset + Enviro + Lat + Long, pca_winfo)
summary(model2)

#checking which indivs are in this clump
pca_plot + geom_label_repel(data = subset(pca_winfo, PC1 < -.05 & PC2 < -.1), aes(label = lab), nudge_x = .1, box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps = 20)
pca_plot

#P1_Ag_1, 5015, P1_Ag_7, P1_Ag_10, 5032, P1_Ag_4, P1_Ag_19, 4917, P1_Ag_5, P1_Ag_12, P1_Ag_2, 5333, 5086, P1_Ag_17, 4955, P1_Ag_3, P1_Ag_16, P1_Ag_18, P1_Nat_6
```

#Assessing IBD
```{r}
#More info: https://www.cog-genomics.org/plink/1.9/ibd
dcgm_IBD <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm_IBD.genome")

dcgm_IBD <- dcgm_IBD %>%
   mutate("Dataset_1" = case_when(
    str_detect(FID1, "_") == T ~ "Drought", 
    str_detect(FID1, "_") == F ~ "Common Garden"
  )) %>%
    mutate("lab_1" = case_when(
    Dataset_1 == "Drought" ~ as.character(paste(FID1, IID1, "T", sep = "_")),
    Dataset_1 == "Common Garden" ~ as.character(paste(FID1, sep = "")))) %>%
   mutate("Dataset_2" = case_when(
    str_detect(FID2, "_") == T ~ "Drought", 
    str_detect(FID2, "_") == F ~ "Common Garden"
  )) %>%
    mutate("lab_2" = case_when(
    Dataset_2 == "Drought" ~ as.character(paste(FID2, IID2, "T", sep = "_")),
    Dataset_2 == "Common Garden" ~ as.character(paste(FID2, sep = "")))) 

cgd <- read.table("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/Chp1_DroughtandCommonGardenInfo_Merged-MetaData.tsv", sep="\t",na.strings = c("","NA"),header=T)

IBD_winfo <- dcgm_IBD %>%
  rowwise() %>%
  mutate("Pair_1" = cgd$Pair[which(cgd$Sample_Name == lab_1)]) %>%
  mutate("Pair_2" = cgd$Pair[which(cgd$Sample_Name == lab_2)]) %>%
  mutate("Pair_match" = case_when(
    Pair_1 == Pair_2 ~ T,
    Pair_1 != Pair_2 ~ F
  )) %>%
  mutate("Lat_1" = cgd$Lat[which(cgd$Sample_Name == lab_1)]) %>%
  mutate("Lat_2" = cgd$Lat[which(cgd$Sample_Name == lab_2)]) %>%
  mutate("Envn_1" = cgd$Environment[which(cgd$Sample_Name == lab_1)]) %>%
  mutate("Envn_2" = cgd$Environment[which(cgd$Sample_Name == lab_2)]) %>%
  mutate("Env_1"= case_when(
    Envn_1 == 1 ~ "Nat",
    Envn_1 == 2 ~ "Ag"
  )) %>%
  mutate("Env_2"= case_when(
    Envn_2 == 1 ~ "Nat",
    Envn_2 == 2 ~ "Ag"
  )) %>%
  select(FID1, IID1, FID2, IID2, lab_1, lab_2, Pair_1, Pair_2, PI_HAT, Pair_match, Env_1, Env_2, Lat_1, Lat_2)

#Box plot
ggplot(data = dcgm_IBD) + geom_boxplot(aes(x= PI_HAT))
ggplot(data = IBD_winfo) + geom_boxplot(aes(x= PI_HAT, color = Pair_match))
mean(dcgm_IBD$PI_HAT)

#Histogram
ggplot(data = dcgm_IBD) + geom_histogram(aes(x = PI_HAT))

#Heatmap
heatmap_dcgm <- IBD_winfo %>%
  mutate("x" = lab_1, "y" = lab_2, "value" = PI_HAT) %>%
  arrange(Lat_1) %>%
  select(x,y,value)

pair_dcgm <-  IBD_winfo %>%
  mutate("x" = paste("P",Pair_1, lab_1, sep = "_"), "y" = paste("P",Pair_2, lab_2, sep = "_"), "value" = PI_HAT) %>%
  arrange(Lat_1) %>%
  select(x,y,value)

all <- ggplot(heatmap_dcgm, aes(x = x, y = y, fill = value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(angle = 45, hjust = .5, vjust = 1)) + ggtitle("IBD heatmap (with .95 indiv)")
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/IBD_heatmap_allindivs_dcgm.png", plot = all, device = "png", height = 14, width = 14, units = "in")


heatmap_dcgm <- heatmap_dcgm[-203,] #remove the one that is .95 so that heatmap gives more info

all_minus <- ggplot(heatmap_dcgm, aes(x = x, y = y, fill = value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(angle = 45, hjust = .5, vjust = 1)) + ggtitle("IBD heatmap (without .95 indiv)")
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/IBD_heatmap_dcgm.png", plot = all_minus, device = "png", height = 14, width = 14, units = "in")

all_pair <- ggplot(pair_dcgm, aes(x = x, y = y, fill = value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(angle = 45, hjust = .5, vjust = 1)) + ggtitle("IBD heatmap (with .95 indiv)")
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/IBD_heatmap_allindivs_pair_dcgm.png", plot = all_pair, device = "png", height = 14, width = 14, units = "in")

#only pairs with Pi_hat > .2
dcgm_.2 <- filter(pair_dcgm, value > .2)
ggplot(dcgm_.2, aes(x = x, y = y, fill = value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(angle = 0, hjust = .5, vjust = 1))

#removing an indiv from pairs where the within pair ibd is high
pt2_ibd <- IBD_winfo %>%
  filter(PI_HAT > .2)
#looked at these and removed based on repeats, then if diff pairs rem one from larger sample size pop, then last 10 remove randomly
#remfromwhichcol <- sample(c(1,2), 10, replace = T)
#Indivs to remove: 
rems <- c('4875','4828','4933','5171','5203','5244','5247','5381','556','P16_Ag_17_T','P17_Ag_4_T','P17_Nat_19_T','P4_Nat_16_T','P17_Nat_7_T','5066','5202','5209','4830','P4_Ag_26_T','4979','P7_Nat_1_T','5084','5220','5258','P4_Ag_1_T','P17_Ag_7_T','P4_Nat_3_T')
#looking at pop stats to adjust pop sizes in chp1_metadata tab 2
cgd_rem <- cgd[cgd$Sample_Name %in% rems,]
cgd_rem <- cgd_rem %>%
  mutate("Env"= case_when(
    Environment == 1 ~ "Nat",
    Environment == 2 ~ "Ag"
  )) %>%
  select(Sample_Name, Pair, Env)
```

#CMH scan after related individuals removed
#p.fdr set at 10% (confirmed same threshold at Science 2022 paper)
```{r}
library(FDRestimation)
library(topr)

filtered <- function(cmh_file, name, FDR_cut = .1){
  FDR <- p.fdr(pvalues = cmh_file$P, threshold = 0.1, adjust.method = "BH", na.rm = F)
  ps <- unlist(FDR[[2]][,2])
  FDR_cmh <- add_column(cmh_file, "FDR_p" = ps)
  rm(FDR, ps)
  print("ps")
  
  FDR_bon <- p.fdr(pvalues = cmh_file$P, threshold = 0.1, adjust.method = "Bon", na.rm = F)
  psb <- unlist(FDR_bon[[2]][,2])
  FDR_cmh <- add_column(FDR_cmh, "Bon_p" = psb)
  remove(FDR_bon, psb, cmh_file)
  print("psb")
  
  write.table(FDR_cmh, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/FDR_", name, sep = "") , row.names = F, quote = F)

  FDR_cmh_filter <- FDR_cmh[which(FDR_cmh$FDR_p < .05),]
  write.table(FDR_cmh_filter, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/P.05FDR_", name, sep = ""), row.names = F, quote = F)
}

man_plot <- function(cmh_file, name){
  #p_fdr_cut <- cmh_file %>% #look for highest raw pval that corresponds to .05 FDR
  #  filter(FDR_p < .05) 
  #p_fd <- max(p_fdr_cut$P)
  
  p_bon_cut <- dcgm %>% #look for highest raw pval that corresponds to .05 Bon
    filter(Bon_p < .05) 
  p_bc <- max(p_bon_cut$P)
    
  cmh_file <- cmh_file %>%
    filter(P < .1) %>%
    separate(CHR, into = c("Sc", "Scaffold"), sep = "_", remove = T) %>%
    mutate("CHROM" = as.numeric(Scaffold), "POS" = BP) %>%
    select(CHROM, POS, P)
  
  if(name == "dcgm"){
    name_end <- "Drought Common Garden Merged"
  }
  
  #man_p <- manhattan(cmh_file, sign_thresh = c(p_fd, p_bc), alpha = .3, xaxis_label = "Scaffold", title = name_end)
  man_p <- manhattan(cmh_file, sign_thresh = p_bc, alpha = .3, xaxis_label = "Scaffold", title = name_end)
  ggsave(man_p, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ManhattanP_", name, ".png", sep = ""), device = "png", height = 6, width = 12, units = "in")
  
  chrs <- unique(cmh_file$CHROM)
  for(i in 1:16){
    #man_p <- manhattan(cmh_file, sign_thresh = c(p_fd, p_bc), alpha = .3, xaxis_label = "Scaffold", title = name_end, chr = chrs[i])
    man_p <- manhattan(cmh_file, sign_thresh = p_bc, alpha = .3, xaxis_label = "Scaffold", title = name_end, chr = chrs[i])
    ggsave(man_p, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ManhattanP_", name, "_", chrs[i], ".png", sep = ""), device = "png", height = 6, width = 12, units = "in")
  }
}

#Sorting files by p value
cmh_dcgm <- fread(file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm_cmh.cmh")
filtered(cmh_dcgm, "dcgm")

#Manhattan plots
dcgm <- fread(file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/FDR_dcgm")
man_plot(dcgm, "dcgm")

#Bon raw p cutoff (p_bc) = .000000002995
#FDR raw p cutoff (p_fd) = .002202

#Num snps
dcgm_b <- dcgm %>%
  filter(P <= .000000002995)
dcgm_f <- dcgm %>%
  filter(P <= .002202)
#w bon padj <= .05 3849 loci
#w fdr padj <= .05 734093 loci
```

#Get list of indivs (and their bams) for phasing
```{r}
fam_dcgm <- fread(file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm.fam")
colnames(fam_dcgm) <- c("Name", "Num", "Pair", "a1", "Sex", "Enviro")
fam_dcgm <- fam_dcgm %>%
  mutate("en_num" = fam_dcgm$Enviro) %>%
  rowwise() %>%
  mutate("Enviro"= case_when(
    en_num == 1 ~ "Nat",
    en_num == 2 ~ "Ag"
  )) %>%
  mutate("Dataset" = case_when(
    str_detect(Name, "_") == T ~ "Drought", 
    str_detect(Name, "_") == F ~ "Common Garden"
  )) %>%
    mutate("lab" = case_when(
    Dataset == "Drought" ~ as.character(paste(Name, Num, "T", sep = "_")),
    Dataset == "Common Garden" ~ as.character(paste(Name, sep = ""))))
bam_names <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/bamsdone_samplenames.txt", header = F)

dcgm_samples <- fam_dcgm$lab
bams_dcgm <- bam_names[bam_names$V2 %in% dcgm_samples,]

#makes list of bams for copying onto server and to add to command
write.table(bams_dcgm$V1, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm_bams.txt", col.names = F, row.names = F, quote = F)

#makes sample list for reheader for whatshap
bams_dcgm <- bams_dcgm %>%
  rowwise() %>%
  mutate("V3" = str_split(V1, "_193")[[1]][1]) %>%
  select(V2, V3)
write.table(bams_dcgm, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/rehead_forwhatshap.txt", col.names = F, row.names = F, quote = F)
```

#Clumping
```{r}
#https://speciationgenomics.github.io/ld_decay/
#https://www.biostars.org/p/84443/#88498
library(topr)
man_plot_clump <- function(cmh_file,clump_file, name){
  p_fdr_cut <- cmh_file %>% #look for highest raw pval that corresponds to .05 FDR
    filter(FDR_p < .05)
  p_fd <- max(p_fdr_cut$P)
  
  p_bon_cut <- cmh_file %>% #look for highest raw pval that corresponds to .05 Bon
    filter(Bon_p < .05) 
  p_bc <- max(p_bon_cut$P)
  
  cmh_file <- cmh_file %>%
    filter(P < .1) %>%
    separate(CHR, into = c("Sc", "Scaffold"), sep = "_", remove = T) %>%
    mutate("CHR" = as.numeric(Scaffold), "POS" = BP)
  
  cmh_file <- cmh_file %>%
    mutate("SNP" = paste(CHR, ":", POS, sep = ""))
  
  clump_file <- clump_file %>%
    separate(CHR, into = c("Sc", "Scaffold"), sep = "_", remove = T) %>%
    mutate("CHR" = as.numeric(Scaffold), "POS" = BP)
  
  clump <- cmh_file %>%
    filter(SNP %in% clump_file$SNP)
  
  if(name == "dcgm"){
    name_end <- "Common Garden and Drought with Clumped SNPs highlighted"
  }
  
  cmh_file <- dplyr::select(cmh_file, c(CHR, BP, P))
  clump_file <- dplyr::select(clump_file, c(CHR, BP, P))
  
  #man_p <- manhattan(list(cmh_file, clump_file), color = list(c("#3497A9FF", "#DEF5E5FF"), c("#0B0405FF")), legend_labels = c("CMH test", "Clumped SNPs"), sign_thresh = p_bc, alpha = .3, xaxis_label = "Scaffold", title = name_end, sign_thresh_label_size = 0)
  manhattan(
  cmh_file,
  color = c("#3497A9FF", "#DEF5E5FF"),
  sign_thresh = p_bc,
  alpha = .3,
  xaxis_label = "Scaffold",
  title = name_end,
  sign_thresh_label_size = 0
)

# Then: add clumped SNPs in black
manhattan(
  clump_file,
  color = "#0B0405FF",
  add = TRUE
)
  ggsave(man_p, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/colornew_ManhattanP_clumped_", name, "_16.png", sep = ""), device = "png", height = 6, width = 16, units = "in")
  #man_p <- manhattan(list(cmh_file, clump_file), color = c("darkblue", "red"), legend_labels = c("CMH test", "Clumped SNPs"), sign_thresh = c(p_fd, p_bc), alpha = .3, xaxis_label = "Scaffold", title = name_end)
  #ggsave(man_p, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ManhattanP_clumped_", name, ".png", sep = ""), device = "png", height = 6, width = 12, units = "in")
  #ggsave(man_p, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ManhattanP_clumped_", name, "_16.png", sep = ""), device = "png", height = 6, width = 16, units = "in")
  #ggsave(man_p, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ManhattanP_clumped_", name, "_18.png", sep = ""), device = "png", height = 6, width = 18, units = "in")
  
  #for figure 2
  #library(viridis)
  #man_p_fig2 <- manhattan(list(cmh_file, clump_file), color = c("#21908CFF", "#440154FF"), legend_labels = c("CMH test", "Clumped SNPs"), sign_thresh = NA, alpha = .3, xaxis_label = "Scaffold", title = name_end)
  #final <- man_p_fig2 + geom_hline(aes(yintercept = -log10(p_bc), colour = "#FDE725FF"), linetype = 2)
  #ggsave(final, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Presentations/Figures/Fig2_ManhattanP_clumped_", name, "_16.png", sep = ""), device = "png", height = 6, width = 16, units = "in")
}

man_plot_clump_gg <- function(cmh_file, clump_file, name){
    p_bon_cut <- cmh_file %>% #look for highest raw pval that corresponds to .05 Bon
    filter(Bon_p < .05) 
    p_bc <- max(p_bon_cut$P)
  
    cmh_file <- cmh_file %>%
      filter(P < .1) %>%
      separate(CHR, into = c("Sc", "Scaffold"), sep = "_", remove = T) %>%
      mutate("CHR" = as.numeric(Scaffold), "POS" = BP)
  
    cmh_file <- cmh_file %>%
      mutate("SNP" = paste(CHR, ":", POS, sep = ""))
    
    cmh_file_cum <- cmh_file %>%
        group_by(Scaffold) %>%
        summarise(chr_len = max(POS)) %>%
        arrange(as.integer(Scaffold)) %>%
        mutate(tot = cumsum(chr_len) - chr_len) %>%
        dplyr::select(-chr_len) %>%
        left_join(cmh_file, ., by = c("Scaffold" = "Scaffold")) %>%
        arrange(Scaffold, POS) %>%
        mutate(BPcum = POS + tot) %>% 
        select(CHR, POS, BPcum, SNP, P, Bon_p)
    
    clump_file_cum <- clump_file %>%
        group_by(CHROM) %>%
        summarise(chr_len = max(BP)) %>%
        arrange(as.integer(CHROM)) %>%
        mutate(tot = cumsum(chr_len) - chr_len) %>%
        dplyr::select(-chr_len) %>%
        left_join(clump_file, ., by = c("CHROM" = "CHROM")) %>%
        arrange(CHROM, BP) %>%
        mutate(BPcum = BP + tot) %>% 
        select(CHROM, BP, BPcum) %>%
        mutate("SNP" = paste(CHROM, ":", BP, sep = ""))
    
    axisdf <- cmh_file_cum %>%
        group_by(CHR) %>%
        summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
    clump <- cmh_file_cum %>%
      filter(SNP %in% clump_file_cum$SNP)
  
    if(name == "dcgm"){
      name_end <- "Common Garden and Drought with Clumped SNPs highlighted"
    }
    
    
    
    man_p_ggplo <- ggplot(cmh_file_cum, aes(x = BPcum, y = -log(P, base = 10))) + geom_point(aes(color=as.factor(CHR)), alpha=0.3, size=1.3) + scale_color_manual(values = rep(c("#3497A9FF", "#A0DFB9FF"), 8)) + geom_point(data = clump, aes(x = BPcum, y = -log(P, base = 10)), color= "#0B0405FF", alpha=.3 , size=1.3) + scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center, expand = c(0.02, 0.02)) + scale_y_continuous(expand = c(0,0)) + ggtitle( paste("Manhattan Plot with clump sites", name, sep = "")) + labs(x = "Scaffold") + geom_hline(aes(yintercept = -log(p_bc, base = 10)), linetype = "dashed", color="#0B0405FF") + theme_bw() + theme(axis.text = element_text(size = 15), legend.position = "none")#+ theme(legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank())
    #man_p_ggplo
    
  ggsave(man_p_ggplo, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/colornew_ManhattanP_clumped_", name, "_16.png", sep = ""), device = "png", height = 6, width = 16, units = "in")
  
  #also plotting without the clump loci just for my proposal
  man_p_ggplo_noclump <- ggplot(cmh_file_cum, aes(x = BPcum, y = -log(P, base = 10))) + geom_point(aes(color=as.factor(CHR)), alpha=0.3, size=1.3) + scale_color_manual(values = rep(c("#3497A9FF", "#A0DFB9FF"), 8))  + scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center, expand = c(0.02, 0.02)) + scale_y_continuous(expand = c(0,0)) + ggtitle( paste("Manhattan Plot with clump sites", name, sep = "")) + labs(x = "Scaffold") + geom_hline(aes(yintercept = -log(p_bc, base = 10)), linetype = "dashed", color="#0B0405FF") + theme_bw() + theme(axis.text = element_text(size = 15), legend.position = "none")#+ theme(legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank())
    #man_p_ggplo
    
  ggsave(man_p_ggplo_noclump, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/colornew_ManhattanP_", name, "_16.png", sep = ""), device = "png", height = 6, width = 16, units = "in")
}

dcgm_clump <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/clumped_dcgm.clumped")

#clumping
#P = FDP adj p value
#NSIG = number that are not significant in that clump
dcgm_table <- dcgm_clump %>%
  separate(SNP, into = c("CHROM", "BP_n"), sep = ":", remove = T) %>%
  mutate("SIG" = TOTAL - NSIG) %>%
  filter(SIG >= 1) %>%
  select(CHROM, BP, P, TOTAL, SIG)
write.table(dcgm_table, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clumphits_dcgm.txt", col.names = T, row.names = F, sep = "\t", quote = F)

#Num rows = # of independent loci
clump_loci <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clumphits_dcgm.txt")
dcgm <- fread(file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/FDR_dcgm")

#Make Manhattan Plots
man_plot_clump(dcgm, dcgm_clump, "dcgm")
man_plot_clump_gg(dcgm, clump_loci, "dcgm")
```

#Getting clumping info for #match to fdr 10% of Sci 2022 paper
```{r}
dcgm <- fread(file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/FDR_dcgm")
p_fdr_cut <- filter(dcgm, FDR_p < .1)
p_fd <- max(p_fdr_cut$P)

library(topr)
man_plot_clump_FDR.1 <- function(cmh_file,clump_file, name){
  p_fdr_cut <- cmh_file %>% #look for highest raw pval that corresponds to .1 FDR
    filter(FDR_p < .1)
  p_fd <- max(p_fdr_cut$P)
  
  cmh_file <- cmh_file %>%
    filter(P < .15) %>%
    separate(CHR, into = c("Sc", "Scaffold"), sep = "_", remove = T) %>%
    mutate("CHR" = as.numeric(Scaffold), "POS" = BP)
  
  cmh_file <- cmh_file %>%
    mutate("SNP" = paste(CHR, ":", POS, sep = ""))
  
  clump_file <- clump_file %>%
    separate(CHR, into = c("Sc", "Scaffold"), sep = "_", remove = T) %>%
    mutate("CHR" = as.numeric(Scaffold), "POS" = BP)
  
  clump <- cmh_file %>%
    filter(SNP %in% clump_file$SNP)
  
  if(name == "dcgm_FDR.1"){
    name_end <- "Common Garden and Drought with Clumped SNPs highlighted FDR .1"
  }
  
  cmh_file <- select(cmh_file, c(CHR, BP, P))
  clump_file <- select(clump_file, c(CHR, BP, P))

  man_p <- manhattan(list(cmh_file, clump_file), color = c("darkblue", "red"), legend_labels = c("CMH test", "Clumped SNPs"), sign_thresh = p_fd, alpha = .3, xaxis_label = "Scaffold", title = name_end)
  ggsave(man_p, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ManhattanP_clumped_", name, ".png", sep = ""), device = "png", height = 6, width = 12, units = "in")
}

#run clump command in github
dcgm_clump_FDR.1 <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/clumped_FDR.1_dcgm.clumped")
man_plot_clump_FDR.1(dcgm, dcgm_clump_FDR.1, "dcgm_FDR.1")

#clumping
#P = FDP adj p value
#NSIG = number that are not significant in that clump
dcgm_table <- dcgm_clump_FDR.1 %>%
  separate(SNP, into = c("CHROM", "BP_n"), sep = ":", remove = T) %>%
  mutate("SIG" = TOTAL - NSIG) %>%
  filter(SIG >= 1) %>%
  select(CHROM, BP, P, TOTAL, SIG)
write.table(dcgm_table, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clumphits_FDR.1_dcgm.txt", col.names = T, row.names = F, sep = "\t", quote = F)

#Num rows = # of independent loci
clump_loci_FDR.1 <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clumphits_FDR.1_dcgm.txt")
```

# Checking for matches in gff file
```{r}
gff <- fread(file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/Atub_193_hap2.all.gff")
#scaffold = V1
#source = V3 (gene, mRNA etc)
#start = V4
#stop = V5

#make gene files numeric chrom, start, stop and tab deliminated
gff_gene <- gff %>%
  filter(V3 == "gene") %>%
  select(V1, V4, V5, V9)

#clumped hits and bon .05 cutoff
dcgm_clump <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clumphits_dcgm.txt")
dcgm_clump_bon.05_bed <- dcgm_clump %>%
  mutate("BP_end" = as.numeric(BP+1), "Scaffold" = paste("Scaffold_", CHROM, sep = "")) %>%
  select(Scaffold, BP, BP_end, P)


write.table(gff_gene, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/gffgene.bed", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(dcgm_clump_bon.05_bed, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/dcgm_clump_bon.05.bed", sep = "\t", row.names = F, col.names = F, quote = F)

#Bon padj < .05 - using
#bedtools intersect -a /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/gffgene.bed -b /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/dcgm_clump_bon.05.bed -wo > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/dcgm_clump_bon.05_bedtools.txt

#Grep for PPO
#grep -e Protoporphyrinogen -e PPO /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/dcgm_clump_bon.05_bedtools.txt > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/PPO_grep.txt

#Now have V1-V4 from the gff file and V5-V8 from cmh file, V8 = FDR p val of that site
dcgm_clump_int <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/dcgm_clump_bon.05_bedtools.txt")

#filter these files to pull out V1,2,3, filter down to matching V4s and have a count of the number of snps in that gene and a list of the FDRps and rank orders from those snps
#RANKS CAN BE REPEATED!!!!!!!
dcgm_clump_rank <- dcgm_clump_int %>%
  arrange(V8) %>%
  mutate(Bonp_rank = dense_rank(V8)) %>%
  group_by(V4) %>%
  summarize(Scaffold = unique(V1), Start = min(V2), Stop = min(V3), n_snp = length((V8)), snp_min = min(V8), snp_max = max(V8), rank_min = min(Bonp_rank), rank_max = max(Bonp_rank)) %>%
  rename(Info = V4)

write.table(dcgm_clump_rank, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/dcgm_clump_intersectBONp.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#viewing
dcgm_clump_rank <- fread(file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/dcgm_clump_intersectBONp.txt")

#see num loci < bon .1
dcgm <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/FDR_dcgm")

dcgm_.1 <- filter(dcgm, FDR_p <= .1)

#List of % genes with consistently highest rank adjp
matches_info <- tibble("File" = c("CommonGarden_vcfmerged"), "Num_bonp_loci" = c(length(dcgm_b.05$SNP)), "Num_indep_loci" = c(length(dcgm_clump$BP)), "Num_genes" = c(length(dcgm_clump_rank$Info)))

dcgm_ranksim <- dcgm_clump_rank %>%
  separate(Info, into = c("Info","Note"), sep = ";Note=") %>%
  separate(Note, into = c("Note", "Info2"), sep = ":") %>%
  group_by(Note) %>%
  summarise(Note = unique(Note), n_snp = sum(n_snp), rank_min = min(rank_min))

matches_info <- add_column(matches_info, "Num_unique_Note" = c(length(dcgm_ranksim$Note)))

#Adding snp id
#awk -F"\t" '{OFS="\t"; split($1, array, "_"); print $0, array[2] ":" $6}' /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/dcgm_clump_bon.05_bedtools.txt > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/dcgm_clump_bon.05_bedtools_ID.txt
dcgm_clump_int_ID <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/dcgm_clump_bon.05_bedtools_ID.txt")

dcgm_intID <- dcgm_clump_int_ID %>%
  mutate("SNP" = V10, "Dataset" = "clumped")
dcgm_gene <- merge(dcgm_.1, dcgm_intID, by = "SNP", all.x = T)
dcgm_gene <- dcgm_gene %>%
  mutate("Info" = V4) %>%
  select(CHR, SNP, BP, P, FDR_p, Bon_p, Info, Dataset)
write.table(dcgm_gene, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/dcgm_.1ID", sep = "\t", row.names = F, col.names = T, quote = F)
```

#manhattan plots with top hits
```{r}
#Copied man plot down and editted to add top hits so I can add in the ranks and similar to info
library(topr)

man_plotlabs <- function(cmh_file, name){
  p_bon_cut <- cmh_file %>% #look for highest raw pval that corresponds to .05 Bon
    filter(Bon_p < .05) 
  p_bc <- max(p_bon_cut$P)
  
  ann_p <- cmh_file$P[which(cmh_file$Bonp_rank == 419)]
  
  cmh_file <- cmh_file %>%
    mutate("Gene_Symbol" = case_when(
    Gene_Symbol == "Protein of unknown function;" ~ NA,
    Gene_Symbol != "Protein of unknown function;" ~ Gene_Symbol
  ))
  
  cmh_file <- cmh_file %>%
    separate(CHR, into = c("Sc", "Scaffold"), sep = "_", remove = T) %>%
    mutate("CHROM" = as.numeric(Scaffold), "POS" = BP) %>%
    select(CHROM, POS, P, Gene_Symbol)
  
  name_end <- "Common Garden Drought CMH"
  
  man_p <- manhattan(cmh_file, sign_thresh = p_bc, alpha = .3, xaxis_label = "Scaffold", title = name_end, annotate = ann_p)
  ggsave(man_p, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ManhattanP_", name, ".png", sep = ""), device = "png", height = 6, width = 12, units = "in")
}

dcgm_clump_ID_rank <- fread(file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/dcgm_.1ID")

dcgm_rank <- dcgm_clump_ID_rank %>%
  arrange(Bon_p) %>%
  mutate(Bonp_rank = dense_rank(Bon_p)) %>%
  separate(Info, into = c("Info","Note"), sep = ";Note=") %>%
  separate(Note, into = c("Gene_Symbol", "Info2"), sep = ":") %>%
  select(CHR, SNP, BP, P, Bon_p, Bonp_rank, Gene_Symbol, Dataset)

man_plotlabs(dcgm_rank, "dcgm_ID")
```

#Make sig sites bed file for pixy
```{r}
dcgm_.05 <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/P.05FDR_dcgm")
dcgm_sig <- filter(dcgm_.05, Bon_p < .05)
dcgm_clump <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clumphits_dcgm.txt") 
dcgm_clump <- dcgm_clump %>%
  mutate("SNP" = paste(CHROM, ":", BP, sep = ""))

#bed file of 10 kb window around each snp
dcgm_sigb <- dcgm_sig %>%
  filter(SNP %in% dcgm_clump$SNP) %>%
  mutate("START" = BP - 5000, "STOP" = BP + 5000) %>%
  rowwise() %>%
  mutate("Num_chrom" = as.integer(str_split(CHR, "_")[[1]][2])) %>%
  arrange(Num_chrom) %>%
  select(CHR, START, STOP)

write.table(dcgm_sigb, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/dcgm_sig.bed", sep = "\t", row.names = F, col.names = F, quote = F)
```

#get genome wide pi and fst from pixy
```{r}
pi_file <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy_allsites/pixy_pi_", "forman10kb", ".txt", sep = ""), sep="\t",header=T)

#genomewide_pi <- sum(pi_file$count_diffs, na.rm = T)/sum(pi_file$count_comparisons, na.rm = T)

gen_pi_ag <- sum(pi_file$count_diffs[which(pi_file$pop == "AG")], na.rm = T)/sum(pi_file$count_comparisons[which(pi_file$pop == "AG")], na.rm = T)
gen_pi_nat <- sum(pi_file$count_diffs[which(pi_file$pop == "NAT")], na.rm = T)/sum(pi_file$count_comparisons[which(pi_file$pop == "NAT")], na.rm = T)

fst_file <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy_allsites/pixy_fst_", "forman10kb", ".txt", sep = ""), sep="\t",header=T)

genomewide_fst <- sum(fst_file$avg_wc_fst, na.rm = T)/length(is.na(fst_file$avg_wc_fst) == F)
```

#Plotting genome wide pi and fst
```{r}
library(cowplot)
manhattans <- function(ext){
  for(i in 1:2){
    if(i == 1){
        file <- paste(ext, "pixy_fst_", sep = "")
        name <- "forman10kb"
        inp <- read.table(paste(file, name, ".txt", sep = ""), sep="\t",header=T)
        pix_file <- inp %>%
          separate(chromosome, into = c("Sc", "Scaffold"), sep = "_", remove = T) %>%
          mutate("POS" = window_pos_1, "FST" = avg_wc_fst, Scaffold = as.numeric(Scaffold)) 
        pix_file <- pix_file %>%
          group_by(Scaffold) %>%
          summarise(chr_len = max(POS)) %>%
          mutate(tot = cumsum(chr_len) - chr_len) %>%
          dplyr::select(-chr_len) %>%
          left_join(pix_file, ., by = c("Scaffold" = "Scaffold")) %>%
          arrange(Scaffold, POS) %>%
          mutate(BPcum = POS + tot) %>% 
          dplyr::select(Scaffold, POS, window_pos_1, window_pos_2, BPcum, FST)
        axisdf <- pix_file %>%
          group_by(Scaffold) %>%
          summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
        
        dcgfclump_add <- dcgclumpsites %>%
          rowwise() %>%
          mutate("POS" = pix_file$BPcum[which(pix_file$window_pos_1 < BP & pix_file$window_pos_2 > BP & pix_file$Scaffold == CHROM)], "FST" = pix_file$FST[which(pix_file$window_pos_1 < BP & pix_file$window_pos_2 > BP & pix_file$Scaffold == CHROM)])
        man_p_fst <- ggplot(pix_file, aes(x = BPcum, y = FST)) + geom_point(aes(color=as.factor(Scaffold)), alpha=0.3, size=1.3) + scale_color_manual(values = rep(c("#3497A9FF", "#A0DFB9FF"), 8)) + geom_point(data = dcgfclump_add, aes(x = POS, y = FST), color="#0B0405FF", alpha=.3 , size=1.3) + scale_x_continuous(label = axisdf$Scaffold, breaks= axisdf$center, expand = c(0.02, 0.02)) + scale_y_continuous(expand = c(0,0)) + ggtitle( paste("FST across the genome ", name, sep = "")) + labs(x = "Scaffold") + theme_bw() + theme(axis.text = element_text(size = 15), legend.position = "none") #+ theme(legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank())
        ggsave(man_p_fst, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Newcolor_ManhattanP_FST_", name, ".png", sep = ""), device = "png", height = 4, width = 18, units = "in")
              
        #man_p_fst <- ggplot(pix_file, aes(x = BPcum, y = FST)) + geom_point(aes(color=as.factor(Scaffold)), alpha=0.3, size=1.3) + scale_color_manual(values = rep(c("#00009E", "#A6A6E3"), 8)) + geom_point(data = dcgfclump_add, aes(x = POS, y = FST), color="red", alpha=.3 , size=1.3) + scale_x_continuous(label = axisdf$Scaffold, breaks= axisdf$center, expand = c(0.02, 0.02)) + scale_y_continuous(expand = c(0,0)) + ggtitle( paste("FST across the genome ", name, sep = "")) + labs(x = "Scaffold") + theme_bw() + theme(legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank())
        #ggsave(man_p_fst, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ManhattanP_FST_", name, ".png", sep = ""), device = "png", height = 4, width = 18, units = "in")
        }else{
        file <- paste(ext, "pixy_pi_", sep = "")
        name <- "forman10kb"
        inp <- read.table(paste(file, name, ".txt", sep = ""), sep="\t",header=T)
        pix_file <- inp %>%
          separate(chromosome, into = c("Sc", "Scaffold"), sep = "_", remove = T) %>%
          mutate(Scaffold = as.numeric(Scaffold), "POS" = window_pos_1) 
        pix_ag <- pix_file %>%
          filter(pop == "AG") %>%
          group_by(Scaffold) %>%
          summarise(chr_len = max(POS)) %>%
          mutate(tot = cumsum(chr_len) - chr_len) %>%
          dplyr::select(-chr_len) %>%
          left_join(pix_file, ., by = c("Scaffold" = "Scaffold")) %>%
          arrange(Scaffold, POS) %>%
          mutate(BPcum = POS + tot) %>%
          filter(pop == "AG") %>%
          dplyr::select(Scaffold, POS, window_pos_1, window_pos_2, BPcum, avg_pi)
        pix_nat <- pix_file %>%
          filter(pop == "NAT") %>%
          group_by(Scaffold) %>%
          summarise(chr_len = max(POS)) %>%
          mutate(tot = cumsum(chr_len) - chr_len) %>%
          dplyr::select(-chr_len) %>%
          left_join(pix_file, ., by = c("Scaffold" = "Scaffold")) %>%
          arrange(Scaffold, POS) %>%
          mutate(BPcum = POS + tot) %>% 
          filter(pop == "NAT") %>%
          dplyr::select(Scaffold, POS, window_pos_1, window_pos_2, BPcum, avg_pi)
        
        axisdf_ag <- pix_ag %>%
          group_by(Scaffold) %>%
          summarize(center=( max(BPcum) + min(BPcum) ) / 2 ) 
        axisdf_nat <- pix_nat %>%
          group_by(Scaffold) %>%
          summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
        
        dcgfclump_add <- dcgclumpsites %>%
          rowwise() %>%
          mutate("POS" = pix_ag$BPcum[which(pix_ag$window_pos_1 < BP & pix_ag$window_pos_2 > BP & pix_ag$Scaffold == CHROM)[1]], "Pi_ag" = pix_ag$avg_pi[which(pix_ag$window_pos_1 < BP & pix_ag$window_pos_2 > BP & pix_ag$Scaffold == CHROM)], "Pi_nat" = pix_nat$avg_pi[which(pix_nat$window_pos_1 < BP & pix_nat$window_pos_2 > BP & pix_nat$Scaffold == CHROM)])

        man_p_pi_ag <- ggplot(pix_ag, aes(x = BPcum, y = avg_pi)) + geom_point(aes(color=as.factor(Scaffold)), alpha=0.3, size=1.3) + scale_color_manual(values = rep(c("#3497A9FF", "#A0DFB9FF"), 8)) + ggtitle(paste("Pi across the genome ", name, sep = "")) + geom_point(data = dcgfclump_add, aes(x = POS, y = Pi_ag), color="#0B0405FF", alpha=.3 , size=1.3) + scale_x_continuous(label = axisdf_ag$Scaffold, breaks= axisdf_ag$center, expand = c(0.02, 0.02)) + scale_y_continuous(expand = c(0,0)) + labs(x = "Scaffold", y = "pi") + theme_bw() + theme(axis.text = element_text(size = 15), legend.position = "none") #+ theme(legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank()) + geom_smooth(method = "loess", span = .01) 
        man_p_pi_nat <- ggplot(pix_nat, aes(x = BPcum, y = avg_pi)) + geom_point(aes(color=as.factor(Scaffold)), alpha=0.3, size=1.3) + scale_color_manual(values = rep(c("#3497A9FF", "#A0DFB9FF"), 8)) + geom_point(data = dcgfclump_add, aes(x = POS, y = Pi_nat), color="#0B0405FF", alpha=.3 , size=1.3) + scale_x_continuous(label = axisdf_nat$Scaffold, breaks= axisdf_nat$center, expand = c(0.02, 0.02)) + scale_y_continuous(expand = c(0,0)) + ggtitle(paste("Pi across the genome ", name, sep = "")) + labs(x = "Scaffold", y = "pi") + theme_bw() + theme(axis.text = element_text(size = 15), legend.position = "none")
        ggsave(plot_grid(man_p_pi_ag, man_p_pi_nat, labels = c("AG", "NAT"), nrow = 2), file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Newcolor_ManhattanP_Pi_", name, ".png", sep = ""), device = "png", height = 4, width = 18, units = "in")
    
        #man_p_pi_ag <- ggplot(pix_ag, aes(x = BPcum, y = avg_pi)) + geom_point(aes(color=as.factor(Scaffold)), alpha=0.3, size=1.3) + scale_color_manual(values = rep(c("#00009E", "#A6A6E3"), 8)) + ggtitle(paste("Pi across the genome ", name, sep = "")) + geom_point(data = dcgfclump_add, aes(x = POS, y = Pi_ag), color="red", alpha=.3 , size=1.3) + scale_x_continuous(label = axisdf_ag$Scaffold, breaks= axisdf_ag$center, expand = c(0.02, 0.02)) + scale_y_continuous(expand = c(0,0)) + labs(x = "Scaffold", y = "pi") + theme_bw() + theme(legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank()) + geom_smooth(method = "loess", span = .01)
        #man_p_pi_nat <- ggplot(pix_nat, aes(x = BPcum, y = avg_pi)) + geom_point(aes(color=as.factor(Scaffold)), alpha=0.3, size=1.3) + scale_color_manual(values = rep(c("#00009E", "#A6A6E3"), 8)) + geom_point(data = dcgfclump_add, aes(x = POS, y = Pi_nat), color="red", alpha=.3 , size=1.3) + scale_x_continuous(label = axisdf_nat$Scaffold, breaks= axisdf_nat$center, expand = c(0.02, 0.02)) + scale_y_continuous(expand = c(0,0)) + ggtitle(paste("Pi across the genome ", name, sep = "")) + labs(x = "Scaffold", y = "pi") + theme_bw() + theme(legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank()) + geom_smooth(method = "loess", span = .01)
        #ggsave(plot_grid(man_p_pi_ag, man_p_pi_nat, labels = c("AG", "NAT"), nrow = 2), file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ManhattanP_Pi_", name, ".png", sep = ""), device = "png", height = 4, width = 18, units = "in")
    }
  }
}

ext <- "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy_allsites/"
dcgclumpsites <- read.table("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clumphits_dcgm.txt", header = T)

dcgclumpsites <- dcgclumpsites %>%
  mutate("Scaffold" = paste("Scaffold_", CHROM, sep = ""))

manhattans(ext)
```
#Histogram of genomewide fst and per pair genomewide fst
```{r}
file <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy_allsites/pixy_fst_forman10kb.txt", sep="\t",header=T)

#setting all fst < 0 to 0
file <- file %>%
  mutate("fst_nonneg" = case_when(
    avg_wc_fst < 0 ~ 0,
    avg_wc_fst >= 0 ~ avg_wc_fst
  ))

meanfst <- mean(file$fst_nonneg)

fst_hist <- ggplot(file) + geom_histogram(aes(x = fst_nonneg), colour = "#3497A9FF", fill = "#3497A9FF") + ggtitle("Genomewide FST in 10kb windows") + labs(x = "FST", y = "Count") + theme_bw() + geom_vline(xintercept = meanfst, linetype = 2)

ggsave(fst_hist, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/FST_hist.png", device = "png", height = 6, width = 8, units = "in")

#histograms per pair
fst_all <- tibble("pop1" = character(), "pop2" = character(), "chromosome" = character(), "window_pos_1" = numeric(), "window_pos_2" = numeric(), "avg_wc_fst" = numeric(), "pair" = numeric(), "fst_nonneg" = numeric())

for(p in 1:17){
    file <- fread(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/pixy_fst_pair_", p, "_windows.txt", sep = ""), sep="\t", header=T)
    pix_file <- file %>%
      mutate("pair" = p) %>%
      select(chromosome, window_pos_1, window_pos_2, pair, pop1, pop2, avg_wc_fst) %>%
      mutate("fst_nonneg" = case_when(
        avg_wc_fst < 0 ~ 0,
        avg_wc_fst >= 0 ~ avg_wc_fst
      ))
    fst_all <- fst_all %>%
      add_row("pop1" = pix_file$pop1, "pop2" = pix_file$pop2, "chromosome" = pix_file$chromosome, "window_pos_1" = pix_file$window_pos_1, "window_pos_2" = pix_file$window_pos_2, "avg_wc_fst" = pix_file$avg_wc_fst, "pair" = pix_file$pair, "fst_nonneg" = pix_file$fst_nonneg)
}
  
fst_plot <- fst_all %>%
  filter(!is.na(fst_nonneg)) %>%
  mutate(pair = factor(pair)) %>%
  ggplot(aes(x = fst_nonneg, color = pair, fill = pair)) +
  geom_density(alpha = 0.25, linewidth = 0.6) +
  labs(
    x = "FST in 10kb Windows",
    y = "density",
    color = "pair",
    fill  = "pair"
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave(fst_plot, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/FST_dist_bypair.png", device = "png", height = 6, width = 8, units = "in")

zoom <- fst_plot + xlim(0, .05)
ggsave(zoom, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/FST_dist_bypair_zoom.png", device = "png", height = 6, width = 8, units = "in")

pair_means <- fst_all %>%
  filter(!is.na(fst_nonneg)) %>%
  group_by(pair) %>%
  summarise(mean_fst = mean(fst_nonneg, na.rm = TRUE), .groups = "drop")

#repeat with negative fsts
fst_plot <- fst_all %>%
  filter(!is.na(avg_wc_fst)) %>%
  mutate(pair = factor(pair)) %>%
  ggplot(aes(x = avg_wc_fst, color = pair, fill = pair)) +
  geom_density(alpha = 0.25, linewidth = 0.6) +
  labs(
    x = "FST in 10kb Windows",
    y = "density",
    color = "pair",
    fill  = "pair"
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave(fst_plot, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/FST_dist_bypair_wneg.png", device = "png", height = 6, width = 8, units = "in")

zoom <- fst_plot + xlim(-.01, .05)
ggsave(zoom, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/FST_dist_bypair_zoom_wneg.png", device = "png", height = 6, width = 8, units = "in")

pair_means <- fst_all %>%
  filter(!is.na(avg_wc_fst)) %>%
  group_by(pair) %>%
  summarise(mean_fst = mean(avg_wc_fst, na.rm = TRUE), .groups = "drop")

#get overall mean across all windows and paris
mean(fst_all$avg_wc_fst, na.rm = T)
```
#Histogram of per pair genomewide fst in 50 Mb windows
```{r}
#histograms per pair
fst_all_50 <- tibble("pop1" = character(), "pop2" = character(), "chromosome" = character(), "window_pos_1" = numeric(), "window_pos_2" = numeric(), "avg_wc_fst" = numeric(), "pair" = numeric(), "fst_nonneg" = numeric())

for(p in 1:17){
    file <- fread(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/pixy_fst_pair_", p, "_50Mb_windows.txt", sep = ""), sep="\t", header=T)
    pix_file <- file %>%
      mutate("pair" = p) %>%
      select(chromosome, window_pos_1, window_pos_2, pair, pop1, pop2, avg_wc_fst) %>%
      mutate("fst_nonneg" = case_when(
        avg_wc_fst < 0 ~ 0,
        avg_wc_fst >= 0 ~ avg_wc_fst
      ))
    fst_all_50 <- fst_all_50 %>%
      add_row("pop1" = pix_file$pop1, "pop2" = pix_file$pop2, "chromosome" = pix_file$chromosome, "window_pos_1" = pix_file$window_pos_1, "window_pos_2" = pix_file$window_pos_2, "avg_wc_fst" = pix_file$avg_wc_fst, "pair" = pix_file$pair, "fst_nonneg" = pix_file$fst_nonneg)
}

#repeat with negative fsts -- HOW SHOULD I DO THIS WITH 1 CALC PER CHROM??
fst_plot <- fst_all_50 %>%
  filter(!is.na(avg_wc_fst)) %>%
  mutate(pair = factor(pair)) %>%
  ggplot(aes(x = avg_wc_fst, color = pair, fill = pair)) +
  geom_boxplot(alpha = 0.25, linewidth = 0.6) +
  labs(
    x = "FST in 50mb Windows/for each chromosome",
    color = "pair",
    fill  = "pair"
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

fst_plot <- fst_all_50 %>%
  filter(!is.na(avg_wc_fst)) %>%
  mutate(pair = factor(pair)) %>%
  ggplot(aes(x = avg_wc_fst, color = pair, fill = pair)) +
  geom_density(alpha = 0.25, linewidth = 0.6) +
  labs(
    x = "FST in 50mb Windows/for each chromosome",
    color = "pair",
    fill  = "pair"
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave(fst_plot, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/FST_dist_bypair_wneg_50mb.png", device = "png", height = 6, width = 8, units = "in")

pair_means <- fst_all_50 %>%
  filter(!is.na(avg_wc_fst)) %>%
  group_by(pair) %>%
  summarise(mean_fst = mean(avg_wc_fst, na.rm = TRUE), .groups = "drop")

#get overall mean across all windows and pairs
mean(fst_all_50$avg_wc_fst, na.rm = T)
```

#Pixy from dif pairs- plot scaffs from dif pairs in same plot
```{r}
plot_pairs <- function(ext){
  # Provide path to input. Can be pi or Dxy.
  # NOTE: this is the only line you should have to edit to run this code:
  file <- c(paste(ext, "pixy_dxy_", sep = ""), paste(ext, "pixy_pi_", sep = ""))
  pidif <- tibble("Site" = numeric(), "Greater" = numeric(), "Less" = numeric(), "CHR" = factor(), "Center" = numeric())
  diftable <- tibble("window_pos_1_ag" = numeric(), "window_pos_1_nat" = numeric(), "avg_pi_ag" = numeric(), "avg_pi_nat" = numeric(), "Pair" = character(), "Pair_nat" = character(), "avg_pi_dif" = numeric(), "chrOrder" = factor(), "windowcenter" = numeric(), "P" = numeric())
  
  
  for(i in 1:2){
    for(s in 1:16){
      inp <- data.frame()
      for(p in 1:17){
        file_p <- paste(file[i], "pair_", p, "_scaf_", s, ".txt", sep ="")
        inpadd <- read.table(file_p,sep="\t",header=T)
        inpadd <- add_column(inpadd, "Pair" = paste("P_", p, sep =""))
        inp <- rbind(inp, inpadd)
      }
      startext <- paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/pixy/", "scaf_", s, sep = "")
      #startext <- paste(ext, "scaf_", s, sep = "")
      # Find the chromosome names and order them: first numerical order, then any non-numerical chromosomes
      #   e.g., chr1, chr2, chr22, chrX
      chroms <- unique(inp$chromosome)
      chrOrder <- sort(chroms)
      inp$chrOrder <- factor(inp$chromosome,levels=chrOrder)
      
      #switching from all bon p < .05 sites to just clump sites
      #sig <- filter(dcgsigsites, dcgsigsites$V1 == chroms[1])
      sig <- filter(dcgclumpsites, dcgclumpsites$Scaffold == chroms[1])
      
      # Plot pi for each population found in the input file
      # Saves a copy of each plot in the working directory
      if("avg_pi" %in% colnames(inp)){
        pops <- unique(inp$pop)
        for (po in 1:length(pops)){
          thisPop <- subset(inp, pop == pops[po])
          # Plot stats along all chromosomes:
          popPlot <- ggplot(thisPop, aes(window_pos_1, avg_pi, color= Pair)) +
            geom_line()+
            facet_grid(. ~ chrOrder)+
            labs(title=paste("Pi for population", pops[po]))+
            labs(x="Position of window start", y="Pi")+
            theme_bw()+
            theme(legend.position = "right") + geom_vline(data = sig, aes(xintercept = BP), color = "gray")
          ggsave(paste(startext, "piplot_", pops[po],".png", sep=""), plot = popPlot, device = "png", dpi = 300)
        }
        
        #plot diff in pi between pops
        agPop <- subset(inp, pop == "AG")
        natPop <- subset(inp, pop == "NAT")
        
        difPop <- as.tibble(agPop$window_pos_1)
        difPop <- difPop %>%
          rename("window_pos_1_ag" = value) %>%
          add_column(window_pos_1_nat = natPop$window_pos_1, avg_pi_ag = agPop$avg_pi, avg_pi_nat = natPop$avg_pi, Pair = agPop$Pair, Pair_nat = natPop$Pair, avg_pi_dif = (agPop$avg_pi - natPop$avg_pi), chrOrder = agPop$chrOrder, windowcenter = natPop$window_pos_1 + 5000) %>% 
          filter(windowcenter %in% sig$BP) 
        
        difPop$P <- rep(NA, length(difPop$window_pos_1_ag))

        sites <- unique(difPop$windowcenter)
        for(l in 1:length(sites)){
          rows <- which(difPop$windowcenter == sites[l])
          difPop$P[rows] <- sig$P[which(sig$BP == sites[l])]
        }

        
        #add to all scaffold data frame
        diftable <- add_row(diftable, difPop)
        
        #Add consistent differences into 
        consistent <- difPop %>%
          mutate("Site" = window_pos_1_ag + 5000) %>%
          group_by(Site) %>%
          summarize("Greater" = sum(avg_pi_dif > .005), "Less" = sum(avg_pi_dif < -.005), "CHR" = unique(chrOrder)) %>%
          filter(Site %in% sig$BP)
        pidif <- add_row(pidif, consistent)
        
        # Plot dif along all chromosomes:
        difPlot <- ggplot(difPop, aes(window_pos_1_ag, avg_pi_dif, color= Pair)) +
          geom_line()+
          facet_grid(. ~ chrOrder)+
          labs(title="Difference in Pi")+
          labs(x="Position of window start", y="Ag Pi - Nat Pi")+
          theme_bw()+
          theme(legend.position = "right") + geom_vline(data = sig, aes(xintercept = BP), color = "gray")
          ggsave(paste(startext, "difpiplot.png", sep=""), plot = difPlot, device = "png", dpi = 300) 
          
        } else {
        print("Pi not found in this file")
      }
      
      # Plot Dxy for each combination of populations found in the input file
      # Saves a copy of each plot in the working directory
      if("avg_dxy" %in% colnames(inp)){
          inp$Pair <- factor(inp$Pair, levels = paste0("P_", 1:17))
          # Plot stats along all chromosomes:
          popPlot <- ggplot(inp, aes(window_pos_1, avg_dxy, color= Pair)) +
            geom_line()+
            facet_grid(. ~ chrOrder)+
            labs(title=paste("Dxy for", inp$pop1[[1]], "&", inp$pop2[[1]]))+
            labs(x="Position of window start", y="Dxy")+
            theme(legend.position = "right")+
            theme_bw() + geom_vline(data = sig, aes(xintercept = BP), color = "gray")
          ggsave(paste(startext, "dxyplot_", inp$pop1[[1]], "_", inp$pop2[[1]],".png", sep=""), plot = popPlot, device = "png", dpi = 300)
        }else {
        print("Dxy not found in this file")
        }
    }
  }
  consistentdif <- ggplot(pidif) + geom_point(aes(x = Greater, y = Less), alpha = .2) + ggtitle("Number of Pairs outside of (.005, -.005)") + geom_abline(intercept = 0, slope = 1, color = "red") + theme_bw()
  ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/pixy/consistentdif.png", plot = consistentdif, device = "png", dpi = 300)
  
  return(diftable)
}

dcgclumpsites <- read.table("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clumphits_dcgm.txt", header = T)

dcgclumpsites <- dcgclumpsites %>%
  mutate("Scaffold" = paste("Scaffold_", CHROM, sep = ""))

dif <- plot_pairs("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/")
```

#Pixy dif in pi plots - MUST RUN PREVIOUS CODE CHUNK BEFORE
#pixy_output.R and Pixy_analysis.rmd both have more (not sure why 2?)
```{r}
options(digits = 3)
options(pillar.sigfig = 3)
options(scipen = 0)
#Looking for consistent reduction in diversity
dif$log_pi_ag <- log10(dif$avg_pi_ag)
dif$log_pi_nat <- log10(dif$avg_pi_nat)

#check for normal distribution
pis <- c(dif$avg_pi_ag, dif$avg_pi_nat)
ggplot() + geom_histogram(aes(x = pis)) + ggtitle("Distribution of pis before transformation")
qqnorm(pis)
qqline(pis)

logpis <- c(dif$log_pi_ag, dif$log_pi_nat)
logpis <- logpis[is.finite(logpis)]
ggplot() + geom_histogram(aes(x = logpis)) + ggtitle("Distribution of pis after log10")
qqnorm(logpis)
qqline(logpis)
qqnorm(dif$log_pi_nat)
qqline(dif$log_pi_nat)
qqnorm(dif$log_pi_ag[is.finite(dif$log_pi_ag)])
qqline(dif$log_pi_ag[is.finite(dif$log_pi_ag)])

dif$exp_pi_ag <- (dif$avg_pi_ag)^.5
dif$exp_pi_nat <- (dif$avg_pi_nat)^.5
exppis <- c(dif$exp_pi_nat, dif$exp_pi_ag)
ggplot() + geom_histogram(aes(x = exppis)) + ggtitle("Distribution of pis after exponential transformation")
qqnorm(exppis)
qqline(exppis)

exp_nat_agdif <- dif %>%
  group_by(windowcenter) %>%
  summarise("Mean_nat" = mean(exp_pi_nat), "Mean_ag" = mean(exp_pi_ag), "Median_nat" = median(exp_pi_nat), "Median_ag" = median(exp_pi_ag), "sd_nat" = sd(exp_pi_nat, na.rm = T), "sd_ag" = sd(exp_pi_ag, na.rm = T), "N_pair" = n(), "se_nat" = sd_nat/sqrt(N_pair), "se_ag" = sd_ag/sqrt(N_pair), "CIupper_ag" = 1.96*se_ag + Mean_ag, "CIlower_ag" = Mean_ag - 1.96*se_ag, "CIupper_nat" = Mean_nat + 1.96*se_nat, "CIlower_nat" = Mean_nat - 1.96*se_nat, "CHR" = unique(chrOrder))

ggplot(exp_nat_agdif, aes(x = Mean_ag, y = Mean_nat)) + geom_point() + geom_smooth(method="loess") + geom_abline(intercept = 0, slope = 1, linetype = 2) + geom_errorbar(aes(ymin = (Mean_nat - se_nat), ymax = (Mean_nat + se_nat))) + geom_errorbar(aes(xmin = (Mean_ag - se_ag), xmax = (Mean_ag + se_ag))) + ggtitle("Mean Pi in Ag and Nat by site after exp transformation")

exp_nat_agdif <- exp_nat_agdif %>%
  mutate("y_overlap" = case_when(
    Mean_ag >= CIlower_nat & Mean_ag <= CIupper_nat ~ "Overlap",
    Mean_ag < CIlower_nat & Mean_ag < CIupper_nat ~ "Above",
    Mean_ag > CIlower_nat & Mean_ag > CIupper_nat ~ "Below"), "x_overlap" = case_when(
      Mean_nat >= CIlower_ag & Mean_nat <= CIupper_ag ~ "Overlap",
      Mean_nat < CIlower_ag & Mean_nat < CIupper_ag ~ "Below",
      Mean_nat > CIlower_ag & Mean_nat > CIupper_ag ~ "Above"
    )) %>% mutate("Selected_in" = case_when(
      x_overlap == "Overlap" & y_overlap == "Overlap" ~ "Both",
      x_overlap == "Above" & y_overlap == "Above" ~ "Nat",
      x_overlap == "Below" & y_overlap == "Below" ~ "Ag"
    ))

plot_leg <- exp_nat_agdif %>%
  group_by(Selected_in) %>%
  summarise("sample" = n()) %>%
  mutate("label" = paste(Selected_in, " (n = ", sample, ")", sep = ""))

exp_nat_agdif <- left_join(exp_nat_agdif, plot_leg, by = "Selected_in")

na_all <- ggplot(filter(exp_nat_agdif, is.na(exp_nat_agdif$Selected_in) == F), aes(x = Mean_ag, y = Mean_nat)) + geom_point(data = filter(exp_nat_agdif, exp_nat_agdif$Selected_in != "Both"), aes(color = label)) + geom_point(data = filter(exp_nat_agdif, exp_nat_agdif$Selected_in == "Both"), aes(color = label), alpha = .6) + geom_abline(intercept = 0, slope = 1, linetype = 2) + geom_errorbar(data = filter(exp_nat_agdif, exp_nat_agdif$Selected_in != "Both"), aes(ymin = CIlower_nat, ymax = CIupper_nat, color = label)) + geom_errorbar(data = filter(exp_nat_agdif, exp_nat_agdif$Selected_in != "Both"), aes(xmin = CIlower_ag, xmax = CIupper_ag, color = label)) + ggtitle("Mean Pi in Ag and Nat by site exp scale \nCIs either both overlap or both don't both overlap") + scale_color_manual(values = c("#60CEACFF","#395D9CFF","#382A54FF")) + xlim(0,.4) + ylim(0,.4) + theme_bw()
#+ geom_errorbar(data = filter(exp_nat_agdif, exp_nat_agdif$Selected_in == "Both"), aes(ymin = CIlower_nat, ymax = CIupper_nat, color = label)) + geom_errorbar(data = filter(exp_nat_agdif, exp_nat_agdif$Selected_in == "Both"), aes(xmin = CIlower_ag, xmax = CIupper_ag, color = label)) 
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/pixy/divplot_expmean_byPair_selected_naall.png", plot = na_all, device = "png", dpi = 300, height = 6, width = 10, units = "in")

#plot empty plot for proposal 
empty_na_all <- ggplot(filter(exp_nat_agdif, is.na(exp_nat_agdif$Selected_in) == F), aes(x = Mean_ag, y = Mean_nat)) + geom_abline(intercept = 0, slope = 1, linetype = 2) + xlim(0,.4) + ylim(0,.4) + ggtitle("Mean Pi in Ag and Nat by site exp scale \nCIs either both overlap or both don't both overlap") + scale_color_manual(values = c("#60CEACFF","#395D9CFF","#382A54FF")) + theme_bw()

ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/pixy/divplot_empty_naall.png", plot = empty_na_all, device = "png", dpi = 300, height = 6, width = 10, units = "in")
```

#get 0 and 4 fold degeneracy average pi and theta
```{r}
pi_file_0 <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy_allsites/pixy_pi_", "zerofold", ".txt", sep = ""), sep="\t",header=T)

pi_file_4 <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy_allsites/pixy_pi_", "fourfold", ".txt", sep = ""), sep="\t",header=T)

gen_pi_ag_0 <- sum(pi_file_0$count_diffs[which(pi_file_0$pop == "AG")], na.rm = T)/sum(pi_file_0$count_comparisons[which(pi_file_0$pop == "AG")], na.rm = T)
gen_pi_nat_0 <- sum(pi_file_0$count_diffs[which(pi_file_0$pop == "NAT")], na.rm = T)/sum(pi_file_0$count_comparisons[which(pi_file_0$pop == "NAT")], na.rm = T)

gen_pi_ag_4 <- sum(pi_file_4$count_diffs[which(pi_file_4$pop == "AG")], na.rm = T)/sum(pi_file_4$count_comparisons[which(pi_file_4$pop == "AG")], na.rm = T)
gen_pi_nat_4 <- sum(pi_file_4$count_diffs[which(pi_file_4$pop == "NAT")], na.rm = T)/sum(pi_file_4$count_comparisons[which(pi_file_4$pop == "NAT")], na.rm = T)


theta_file_4 <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy_allsites/pixy_watterson_theta_", "fourfold", ".txt", sep = ""), sep="\t",header=T)
theta_ag_4 <- theta_file_4$avg_watterson_theta[which(theta_file_4$pop == "AG")]
gen_theta_ag_4 <- sum(theta_ag_4, na.rm = T)/length(is.na(theta_ag_4) == F)
theta_nat_4 <- theta_file_4$avg_watterson_theta[which(theta_file_4$pop == "NAT")]
gen_theta_nat_4 <- sum(theta_nat_4, na.rm = T)/length(is.na(theta_nat_4) == F)
```

#selscan output
```{r}
#no window analysis (with normalization run across all chromsomes and not individually)
#xpehh <- tibble("locus_id" = character(), "POS" = integer(), "CHROM" = integer(), "gpos" = numeric(), "popA_freq" = numeric(), "ihhA" = numeric(), "popB_freq"= numeric(), "ihhB" = numeric(), "raw_xpehh"= numeric(), "norm_xpehh"= numeric(), "crit" = integer(), "P" = numeric())

#for(s in 1:16){
#  tmp <- fread(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/sel_scan/Scaffold_", s,".xpehh.out.norm", sep = ""))
#  colnames(tmp) <- c("locus_id", "POS", "gpos", "popA_freq", "ihhA", "popB_freq", "ihhB", "raw_xpehh", "norm_xpehh", "crit")
#  tmp$CHROM <- s
  
#  xpehh <- add_row(xpehh, tmp)
#}

#write.table(xpehh, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/sel_scan/Aggregated_xpehh.txt", row.names = F, quote = F)

xpehh <- read.table("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/sel_scan/Aggregated_xpehh.txt", header = T)

#ggplot(xpehh) + geom_histogram(aes(norm_xpehh))

#plot selscan
xpehh_file <- xpehh %>%
  mutate("direction" = case_when(
    norm_xpehh > 0 ~ "Sel_in_Ag",
    norm_xpehh < 0 ~ "Sel_in_Nat"))
    
xpehh_file_cum <- xpehh_file %>%
    group_by(CHROM) %>%
    summarise(chr_len = max(POS)) %>%
    arrange(as.integer(CHROM)) %>%
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    dplyr::select(-chr_len) %>%
    left_join(xpehh_file, ., by = c("CHROM" = "CHROM")) %>%
    arrange(CHROM, POS) %>%
    mutate(BPcum = POS + tot) %>% 
    select(CHROM, POS, BPcum, locus_id, norm_xpehh, crit)
    
axisdf <- xpehh_file_cum %>%
      group_by(CHROM) %>%
      summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

clump_loci <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clumphits_dcgm.txt")
clump_loci$SNP <- as.character(paste(clump_loci$CHROM, ":", clump_loci$BP, sep = ""))

dcgclump_add <- clump_loci %>%
  rowwise() %>%
  mutate("BPcum" = xpehh_file_cum$BPcum[which(xpehh_file_cum$locus_id == SNP)], "norm_xpehh" = xpehh_file_cum$norm_xpehh[which(xpehh_file_cum$locus_id == SNP)])

man_p_ggplo <- ggplot(xpehh_file_cum, aes(x = BPcum, y = norm_xpehh)) + geom_point(aes(color=as.factor(CHROM)), alpha=0.3, size=1.3) + scale_color_manual(values = rep(c("#3497A9FF", "#A0DFB9FF"), 8)) + scale_x_continuous(label = axisdf$CHROM, breaks= axisdf$center, expand = c(0.02, 0.02)) + ggtitle( "XPEHH across the genome") + labs(x = "Scaffold") + geom_hline(aes(yintercept = 5), linetype = "dashed", color="#0B0405FF") + geom_hline(aes(yintercept = -5), linetype = "dashed", color="#0B0405FF") + geom_point(data = dcgclump_add, aes(x = BPcum, y = norm_xpehh), color="#0B0405FF", alpha=.3 , size=1.3) + theme_bw() + theme(axis.text = element_text(size = 15), legend.position = "none") + ylim(-10,10)
    #man_p_ggplo
    
ggsave(man_p_ggplo, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/colornew_selscan_16.png", sep = ""), device = "png", height = 6, width = 16, units = "in")

#look for intercept of selscan and cmh hits - currently looking at the LD thinned hits
#crit is 1 or -1 if extreme, 0 if not extreme
critical_xpehh <- xpehh %>%
  filter(norm_xpehh > 5 | norm_xpehh < -5)
  #filter(crit == 1 | crit == -1) 
  
overlap <- merge(clump_loci, critical_xpehh, by.x = "SNP", by.y = "locus_id")

#checking for PPO
dcgm_clump_int_ID <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/dcgm_clump_bon.05_bedtools_ID.txt")
#searched POX2 in that table
PPO <- c("8:10446220", "8:10444061", "8:10438798", "8:10448825")
ppo_selscan_hits <- PPO %in% overlap$SNP
#"8:10444061", "8:10438798" in the hits

#window analysis - not going to analyze this
#xpehh_win <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/sel_scan/Scaffold_8.xpehh.out.norm.100kb.windows")
#colnames(xpehh_win) <- c("start", "end", "num_snps", "fraction_extremesnps", "percentile")
```

#Look at how XPEHH correlates with recomb rate
``` {r}
xpehh <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/sel_scan/Aggregated_xpehh.txt", header = T)
map <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm_cM.map")
colnames(map) <- c("Scaffold", "SNP", "cM", "Locus")
map$chr <- as.numeric(gsub("Scaffold_([0-9]+).*", "\\1", map$Scaffold))

#covert cM to r (cM in rhotocM_forlibby_new.R is just r*100)
map$r <- map$cM/100

#make recomb rate bins
bins <- quantile(map$r)

#assign recomb rates to these bins and set window positions
rec_bins <- map %>%
  group_by(chr) %>%
  mutate("win_pos_1" = case_when(
    Locus == min(Locus) ~ 0, 
    Locus > min(Locus) ~ lag(Locus))) %>%
  mutate("win_pos_2" = Locus) %>%
  select(chr, SNP, win_pos_1, win_pos_2, r) %>%
  ungroup()

#get rec rate for each of our xpehh scores
#rec_bins_wr <- xpehh %>% - this didn't work so ran next chunk that does same thing
#  rowwise() %>%
#  mutate("rec_rate" = rec_bins$r[which(rec_bins$chr == CHROM & rec_bins$win_pos_1 <= POS & rec_bins$win_pos_2 >= POS)[1]])

rec_bins_wr <- xpehh %>%
  left_join(
    rec_bins %>% select(chr, win_pos_1, win_pos_2, r),
    join_by(CHROM == chr, POS >= win_pos_1, POS <= win_pos_2)
  ) %>%
  rename(rec_rate = r) %>%
  slice_head(n = 1, by = locus_id)

#assign these to their rec rate bin
rec_bins_wr <- rec_bins_wr %>%
  mutate("bin" = case_when(
    rec_rate <= bins[1] ~ "1",
    rec_rate > bins[1] & rec_rate <= bins[2] ~ "2",
    rec_rate > bins[2] & rec_rate <= bins[3] ~ "3",
    rec_rate > bins[3] & rec_rate <= bins[4] ~ "4",
    rec_rate > bins[4] ~ "5"
  ))

norm_rec <- ggplot(rec_bins_wr) + geom_boxplot(aes(y = norm_xpehh, group = bin)) + xlab("Recombination rate") + labs(title = "Recombination rate binned by quantile \n with normalized xpehh values in that rec rate")

raw_rec <- ggplot(rec_bins_wr) + geom_boxplot(aes(y = raw_xpehh, group = bin)) + xlab("Recombination rate") + labs(title = "Recombination rate binned by quantile \n with raw xpehh values in that rec rate")

ggsave(norm_rec, filename = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/recomb_rate_wnormxpehh.png")
ggsave(raw_rec, filename = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/recomb_rate_wrawxpehh.png")
```
#selscan BY PAIR - run on randi because too much data to download
```{r}
#module load gcc/12.1.0
#module load R (version 4.2.1)
#R
library(data.table)
library(tidyverse)
xpehh <- tibble("pair" = numeric(), "locus_id" = character(), "POS" = integer(), "CHROM" = integer(), "gpos" = numeric(), "popA_freq" = numeric(), "ihhA" = numeric(), "popB_freq"= numeric(), "ihhB" = numeric(), "raw_xpehh"= numeric(), "norm_xpehh"= numeric(), "crit" = integer(), "P" = numeric())

for(p in 1:17){
 for(s in 1:16){
  tmp <- fread(paste("/scratch/espolston/sel_scan_pair/Scaffold_", s,"_pair_", p, ".xpehh.out.norm", sep = ""))
  colnames(tmp) <- c("locus_id", "POS", "gpos", "popA_freq", "ihhA", "popB_freq", "ihhB", "raw_xpehh", "norm_xpehh", "crit")
  tmp$CHROM <- s
  tmp$pair <- p
  
  xpehh <- add_row(xpehh, tmp)
 } 
  write.table(tmp, paste("/scratch/espolston/sel_scan_pair/Aggregated_pair_", p, "_xpehh.txt", sep =""), row.names = F, quote = F)
}


clump_loci <- fread("/scratch/espolston/clumphits_dcgm.txt")
clump_loci$SNP <- as.character(paste(clump_loci$CHROM, ":", clump_loci$BP, sep = ""))

for(p in 3:17){
xpehh <- read.table(paste("/scratch/espolston/sel_scan_pair/Aggregated_pair_", p, "_xpehh.txt", sep =""), header = T)

#plot selscan
min(xpehh$norm_xpehh)
max(xpehh$norm_xpehh)

xpehh_file_cum <- xpehh %>%
    group_by(CHROM) %>%
    summarise(chr_len = max(POS)) %>%
    arrange(as.integer(CHROM)) %>%
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    dplyr::select(-chr_len) %>%
    left_join(xpehh, ., by = c("CHROM" = "CHROM")) %>%
    arrange(CHROM, POS) %>%
    mutate(BPcum = POS + tot) %>% 
    select(CHROM, POS, BPcum, locus_id, norm_xpehh, crit, pair)

axisdf <- xpehh_file_cum %>%
      group_by(CHROM) %>%
      summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

dcgclump_add <- clump_loci %>%
  rowwise() %>%
  mutate("BPcum" = xpehh_file_cum$BPcum[which(xpehh_file_cum$locus_id == SNP)], "norm_xpehh" = xpehh_file_cum$norm_xpehh[which(xpehh_file_cum$locus_id == SNP)])  

xpehh_file_cum <- bind_rows(
  xpehh_file_cum %>% filter(norm_xpehh <= -5 | norm_xpehh >= 5),
  xpehh_file_cum %>% filter(norm_xpehh > -5 & norm_xpehh < 5) %>% slice_sample(prop = 0.1)
)


  man_p_ggplo <- ggplot(xpehh_file_cum %>% filter(pair == p), aes(x = BPcum, y = norm_xpehh)) + geom_point(aes(color=as.factor(CHROM)), alpha=0.3, size=1.3) + scale_color_manual(values = rep(c("#3497A9FF", "#A0DFB9FF"), 8)) + scale_x_continuous(label = axisdf$CHROM, breaks= axisdf$center, expand = c(0.02, 0.02)) + ggtitle( "XPEHH across the genome") + labs(x = "Scaffold") + geom_hline(aes(yintercept = 5), linetype = "dashed", color="#0B0405FF") + geom_hline(aes(yintercept = -5), linetype = "dashed", color="#0B0405FF") + theme_bw() + theme(axis.text = element_text(size = 15), legend.position = "none") + ylim(-20,20)
    
  ggsave(man_p_ggplo, file = paste("/scratch/espolston/sel_scan_pair/colornew_selscan_16_pair_", p, "_noclump.png", sep = ""), device = "png", height = 6, width = 16, units = "in")

  man_p_ggplo <- ggplot(xpehh_file_cum %>% filter(pair == p), aes(x = BPcum, y = norm_xpehh)) + geom_point(aes(color=as.factor(CHROM)), alpha=0.3, size=1.3) + scale_color_manual(values = rep(c("#3497A9FF", "#A0DFB9FF"), 8)) + scale_x_continuous(label = axisdf$CHROM, breaks= axisdf$center, expand = c(0.02, 0.02)) + ggtitle( "XPEHH across the genome") + labs(x = "Scaffold") + geom_hline(aes(yintercept = 5), linetype = "dashed", color="#0B0405FF") + geom_hline(aes(yintercept = -5), linetype = "dashed", color="#0B0405FF") + geom_point(data = dcgclump_add, aes(x = BPcum, y = norm_xpehh), color="#0B0405FF", alpha=.3 , size=1.3) + theme_bw() + theme(axis.text = element_text(size = 15), legend.position = "none") + ylim(-20,20)
    #man_p_ggplo
    
  ggsave(man_p_ggplo, file = paste("/scratch/espolston/sel_scan_pair/colornew_selscan_16_pair_", p, ".png", sep = ""), device = "png", height = 6, width = 16, units = "in")
}
```

#look at H12 overlaps - after controlling for recomb rate (dropping everything in first quartile of recomb rate)
```{r}
#Move on to this after running peaks after removing the low recomb regions
#REDO THIS TO CONSIDER TOP 50 PEAKS GENOME WIDE
#also check if these peaks are in areas of possible SVs or missing data- do I already control for missing data enough??
H12Scan_all <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12_H2H1/dcgm_hap_stats_Scaffold_", 1,"_recombfiltered.txt", sep = ""))
colnames(H12Scan_all) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2/H1", "H123")
peaks_all <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12_H2H1/dcgm_hap_peaks_Scaffold_", 1,"_recombfiltered.txt", sep = ""))
colnames(peaks_all) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2/H1", "H123","smallest_edge_coordofpeak", "smallest_edge_coordofpeak")

#get the top 10 peaks for each chrom
peaks_all$top_peak <- c(rep(T, 10), rep(F, length(peaks_all$win_center) - 10))
peaks_all$chr <- rep(1, length(peaks_all$win_center))
H12Scan_all$chr <- rep(1, length(H12Scan_all$win_center))
H12Scan_all$bpcum <- H12Scan_all$win_center
peaks_all$bpcum <- peaks_all$win_center

for(s in 2:16){
  H12Scan_all_in = read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12_H2H1/dcgm_hap_stats_Scaffold_", s,"_recombfiltered.txt", sep = ""))
  peaks_all_in = read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12_H2H1/dcgm_hap_peaks_Scaffold_", s,"_recombfiltered.txt", sep = ""))
  colnames(H12Scan_all_in) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2/H1", "H123")
  colnames(peaks_all_in) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2/H1", "H123","smallest_edge_coordofpeak", "smallest_edge_coordofpeak")
  
  peaks_all_in$top_peak <- c(rep(T, 10), rep(F, length(peaks_all_in$win_center) - 10))
  peaks_all_in$chr <- rep(s, length(peaks_all_in$win_center))
  H12Scan_all_in$chr <- rep(s, length(H12Scan_all_in$win_center))
  
  #make bpcum
  prelength <- max(H12Scan_all$bpcum)
  peaks_all_in$bpcum <- peaks_all_in$win_center + prelength
  H12Scan_all_in$bpcum <- H12Scan_all_in$win_center + prelength
  
  H12Scan_all <- rbind(H12Scan_all, H12Scan_all_in)
  peaks_all <- rbind(peaks_all, peaks_all_in)
}

#load in ag run
peaks_ag <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12_H2H1/dcgm_hap_peaks_Scaffold_", 1,"_recombfiltered_ag.txt", sep = ""))
H12Scan_ag <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12_H2H1/dcgm_hap_stats_Scaffold_", 1,"_recombfiltered_ag.txt", sep = ""))
colnames(H12Scan_ag) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2/H1", "H123")
colnames(peaks_ag) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2/H1", "H123","smallest_edge_coordofpeak", "smallest_edge_coordofpeak")

#get the top 10 peaks for each chrom
peaks_ag$top_peak <- c(rep(T, 10), rep(F, length(peaks_ag$win_center) - 10))
peaks_ag$chr <- rep(1, length(peaks_ag$win_center))
H12Scan_ag$chr <- rep(1, length(H12Scan_ag$win_center))
H12Scan_ag$bpcum <- H12Scan_ag$win_center
peaks_ag$bpcum <- peaks_ag$win_center

for(s in 2:16){
  H12Scan_ag_in = read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12_H2H1/dcgm_hap_stats_Scaffold_", s,"_recombfiltered_ag.txt", sep = ""))
  peaks_ag_in = read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12_H2H1/dcgm_hap_peaks_Scaffold_", s,"_recombfiltered_ag.txt", sep = ""))
  colnames(H12Scan_ag_in) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2/H1", "H123")
  colnames(peaks_ag_in) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2/H1", "H123","smallest_edge_coordofpeak", "smallest_edge_coordofpeak")
  
  peaks_ag_in$top_peak <- c(rep(T, 10), rep(F, length(peaks_ag_in$win_center) - 10))
  peaks_ag_in$chr <- rep(s, length(peaks_ag_in$win_center))
  H12Scan_ag_in$chr <- rep(s, length(H12Scan_ag_in$win_center))
  
  #make bpcum
  prelength <- max(H12Scan_ag$bpcum)
  peaks_ag_in$bpcum <- peaks_ag_in$win_center + prelength
  H12Scan_ag_in$bpcum <- H12Scan_ag_in$win_center + prelength
  
  H12Scan_ag <- rbind(H12Scan_ag, H12Scan_ag_in)
  peaks_ag <- rbind(peaks_ag, peaks_ag_in)
}

#load in nat run
peaks_nat <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12_H2H1/dcgm_hap_peaks_Scaffold_", 1,"_recombfiltered_nat.txt", sep = ""))
H12Scan_nat <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12_H2H1/dcgm_hap_stats_Scaffold_", 1,"_recombfiltered_nat.txt", sep = ""))
colnames(H12Scan_nat) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2/H1", "H123")
colnames(peaks_nat) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2/H1", "H123","smallest_edge_coordofpeak", "smallest_edge_coordofpeak")

#get the top 10 peaks for each chrom
peaks_nat$top_peak <- c(rep(T, 10), rep(F, length(peaks_nat$win_center) - 10))
peaks_nat$chr <- rep(1, length(peaks_nat$win_center))
H12Scan_nat$chr <- rep(1, length(H12Scan_nat$win_center))
H12Scan_nat$bpcum <- H12Scan_nat$win_center
peaks_nat$bpcum <- peaks_nat$win_center

for(s in 2:16){
  H12Scan_nat_in = read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12_H2H1/dcgm_hap_stats_Scaffold_", s,"_recombfiltered_nat.txt", sep = ""))
  peaks_nat_in = read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12_H2H1/dcgm_hap_peaks_Scaffold_", s,"_recombfiltered_nat.txt", sep = ""))
  colnames(H12Scan_nat_in) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2/H1", "H123")
  colnames(peaks_nat_in) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2/H1", "H123","smallest_edge_coordofpeak", "smallest_edge_coordofpeak")
  
  peaks_nat_in$top_peak <- c(rep(T, 10), rep(F, length(peaks_nat_in$win_center) - 10))
  peaks_nat_in$chr <- rep(s, length(peaks_nat_in$win_center))
  H12Scan_nat_in$chr <- rep(s, length(H12Scan_nat_in$win_center))
  
  #make bpcum
  prelength <- max(H12Scan_nat$bpcum)
  peaks_nat_in$bpcum <- peaks_nat_in$win_center + prelength
  H12Scan_nat_in$bpcum <- H12Scan_nat_in$win_center + prelength
  
  H12Scan_nat <- rbind(H12Scan_nat, H12Scan_nat_in)
  peaks_nat <- rbind(peaks_nat, peaks_nat_in)
}

rm(H12Scan_ag_in, H12Scan_nat_in, H12Scan_all_in, peaks_ag_in, peaks_nat_in, peaks_all_in)

#now compare peaks
#add snp
peaks_ag$snp <- paste(peaks_ag$chr, peaks_ag$win_center, sep = ":")
peaks_all$snp <- paste(peaks_all$chr, peaks_all$win_center, sep = ":")
peaks_nat$snp <- paste(peaks_nat$chr, peaks_nat$win_center, sep = ":")

ag_all <- which(peaks_ag$snp %in% peaks_all$snp)
nat_all <- which(peaks_nat$snp %in% peaks_all$snp)
ag_nat <- which(peaks_ag$snp %in% peaks_nat$snp)
ag_nat_all <- which(peaks_ag$snp[ag_nat] %in% peaks_all$snp)

#NOTE!!!!!!!
##h12 need to redo top peaks- look at top 50, do they overlap w inversions, missing data, low recomb regions
```

#AFvapeR to look at amount of parallelism - write here and then move to bbedit and then run on server 
``` {r}
#from: https://github.com/JimWhiting91/afvaper
# install.packages("remotes",repos = "http://cran.us.r-project.org")
# remotes::install_github("JimWhiting91/afvaper")
library(afvaper,verbose = F)
library(vcfR,verbose = F)
vcf_in <- read.vcfR("/Users/libbypolston/Desktop/phased_dcgm_Scaffold_16.vcf.gz",verbose = F)

#making popmap
#awk 'NR > 1 {print $1,$3,$5}' /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/Chp1_DroughtandCommonGardenInfo_Merged-MetaData.tsv > /Users/libbypolston/Desktop/popmap.txt

popmap <- read.table("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm.fam")
popmap <- popmap %>%
  mutate("indiv" = case_when(
    grepl("_", V1) == T ~ paste(V1, "_", V2, sep =""),
    grepl("_", V1) == F ~ V1), 
    "env" = case_when(
      V6 == 1 ~ "NAT",
      V6 == 2 ~ "AG"), "contrast" = paste(V3, "_", env, sep ="")) %>%
  select(indiv, contrast)

#just for this one bc it has AT
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

vector_list <- list(c("12_AG","12_NAT"),
                    c("15_AG","15_NAT"),
                    c("8_AG","8_NAT"),
                    c("5_AG","5_NAT"))
# Name vectors
names(vector_list) <- c("contrast1","contrast2","contrast3","contrast4")

# Set our window size
window_snps = 200

# Calculate Allele Frequency Change Vector Matrices
AF_input <- calc_AF_vectors(vcf = vcf_in,
                            window_size = window_snps,
                            popmap = popmap,
                            vectors = vector_list,
                            n_cores = 4,
                            data_type = "vcf")
print(paste0("Number of windows = ",length(AF_input)))
print(paste0("Number of SNPs per window = ",ncol(AF_input[[1]])))
print(paste0("Number of vectors per window = ",nrow(AF_input[[1]])))


# How many permutations do we want in total?
#total_perms <- 10000

# Imagine a hypothetical fasta index with three chromosomes
#genome_fai <- data.frame(chr=c("Scaffold_1","Scaffold_2","Scaffold_3", "Scaffold_4","Scaffold_5","Scaffold_6","Scaffold_7","Scaffold_8","Scaffold_9","Scaffold_10","Scaffold_11","Scaffold_12","Scaffold_13","Scaffold_14","Scaffold_15","Scaffold_16"),
#                         size=c(2e7,4e7,3e7)) #STILL NEED TO CHANGE SIZES!!!

# Fetch proportional size of all chromosomes
#chr_props <- genome_fai$size/sum(genome_fai$size)
#chr_perms <- data.frame(chr=genome_fai$chr,
#                        perms=round(chr_props * total_perms))

# This gives us approximately 10000 null perms in total, distributed across the genome according to relative size of chromosomes...
#chr_perms
# Calculate Allele Frequency Change Vector Matrices - do this for each chr
#n <- 1 (which chrom)
#null_input <- calc_AF_vectors(vcf = vcf_in_n,
#                              window_size = window_snps,
#                              popmap = popmap,
#                              vectors = vector_list,
#                              n_cores = 4,
#                              null_perms = chr_perms[n],
#                              data_type = "vcf")

# Show features of input...
#print(paste0("Number of null windows = ",length(null_input)))
#print(paste0("Number of SNPs per window = ",ncol(null_input[[1]])))
#print(paste0("Number of vectors per window = ",nrow(null_input[[1]])))

#next 2 steps will need to be replaced by above commented once i get sizes 
null_perm_N = 1000

# Calculate Allele Frequency Change Vector Matrices
null_input <- calc_AF_vectors(vcf = vcf_in,
                              window_size = window_snps,
                              popmap = popmap,
                              vectors = vector_list,
                              n_cores = 4,
                              null_perms = null_perm_N,
                              data_type = "vcf")

###eigen analysis over AF matrices
# Perform eigen analysis
eigen_res <- lapply(AF_input,eigen_analyse_vectors)
# View chromosomal regions:
head(names(eigen_res))
# View eigenvalue distribution of first matrix
eigen_res[[1]]$eigenvals
# View eigenvector loadings of first matrix
eigen_res[[1]]$eigenvecs
# View head of SNP scores
head(eigen_res[[1]]$A_matrix)
# Get cutoffs for 95%, 99% and 99.9%
null_cutoffs <- find_null_cutoff(null_input,cutoffs = c(0.95,0.99,0.999))
null_cutoffs


# Calculate p-vals
pvals <- eigen_pvals(eigen_res,null_input)

# Showpvals
head(pvals)

# Plot the raw eigenvalues, and visualise the cutoff of 99%
all_plots <- eigenval_plot(eigen_res,cutoffs = null_cutoffs[,"99%"])

# Show the plots for eigenvalue 1
all_plots[[1]]
# Plot empirical p-values, -log10(p) of 2 ~ p=0.01, 3 ~ p=0.001 etc.
all_plots_p <- eigenval_plot(eigen_res,null_vectors = null_input,plot.pvalues = T)
# Show the plots for eigenvalue 1
all_plots_p[[1]]
# Plot empirical p-values, -log10(p) of 2 ~ p=0.01, 3 ~ p=0.001 etc.
chr1_windows <- grep("Scaffold_16",names(eigen_res))
all_plots_p_chr1 <- eigenval_plot(eigen_res[chr1_windows],null_vectors = null_input,plot.pvalues = T)
# Show the plots for eigenvalue 1
all_plots_p_chr1[[1]]
library(ggplot2,verbose = F)

# Pull the figure 
eig1_pval_fig <- all_plots_p[[1]]

# Edit
eig1_pval_fig + theme(title = element_blank()) + geom_step(colour="red2")

# Recall the use of find_null_cutoffs() to fetch a matrix of cutoffs...
# null_cutoffs

# Find significant windows above 99.9% null permutation
significant_windows <- signif_eigen_windows(eigen_res,null_cutoffs[,"99.9%"])

# Display 'outliers'
significant_windows

# Summarise parallel evolution in windows that are significant on eigenvector 4
eig4_parallel <- summarise_window_parallelism(window_id = significant_windows[[4]],
                                              eigen_res = eigen_res,
                                              loading_cutoff = 0.3,
                                              eigenvector = 1)
# Show results
head(eig4_parallel)

# Fetch an A matrix
A_mat <- eigen_res[[4]]$A_matrix
head(A_mat)

to_plot <- data.frame(snp=rownames(A_mat),
                      eig4_score=A_mat[,1])

to_plot <- to_plot %>% tidyr::separate("snp",into=c("sc","scaf","pos"),sep="_")
to_plot$pos <- as.integer(to_plot$pos)

ggplot(to_plot,aes(x=pos,y=abs(eig4_score)))+
  geom_point()+
  labs(y="Eig4 Score",x="Pos (bp)")

#this all works but figure out what these plots mean
```

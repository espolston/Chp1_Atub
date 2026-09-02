---
title: "Chp1_mergedvcf"
author: "Libby Polston"
date: "2025-12-19"
output: pdf_document
editor_options: 
  markdown: 
    wrap: 72
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

#MDS and PCA for pre relatedness removal (this fam, bim etc files
overwritten for 440 indiv but mds and pca files not rerun)

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

cgd <- read.table("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/Chp1_DroughtandCommonGardenInfo_MergedMetaData.tsv", sep="\t",na.strings = c("","NA"),header=T)
cols <- c(1,6,7)
cgd <- cgd[,cols]

pca_winfo <- merge(pca_dat, cgd, by.x = "lab", by.y = "Sample_Name", all= T)

#PCA for dataset
library(ggrepel)
library(viridis)
pca_plot <- ggplot(pca_winfo, aes(x = PC1, y = PC2)) + geom_point(aes(color = Pair, shape = Enviro)) + ggtitle("PCA with .01 MAF .05 missingness") + theme_bw() + theme(axis.text = element_text(size = 15)) + scale_color_viridis_d(option = "mako")
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/PCA_dcgm.png", plot = pca_plot, device = "png", height = 6, width = 12, units = "in") 

#repeating for Latitude
pca_plotlat <- ggplot(pca_winfo, aes(x = PC1, y = PC2)) + geom_point(aes(color = Lat, shape = Enviro)) + ggtitle("PCA with .01 MAF .05 missingness") + theme_bw() + theme(axis.text = element_text(size = 15))
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/PCA_dcgm_lat.png", plot = pca_plotlat, device = "png", height = 6, width = 12, units = "in") 

model <- lm(PC1 ~ Dataset + Enviro + Lat + Long, pca_winfo)
summary(model)

model2 <- lm(PC2 ~ Dataset + Enviro + Lat + Long, pca_winfo)
summary(model2)

#checking which indivs are in the clump
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

#CMH scan after related individuals removed #p.fdr set at 10% (confirmed
same threshold at Science 2022 paper)

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
library(ggrepel)
man_plot_clump_gg <- function(cmh_file, clump_file, name){
    p_bon_cut <- cmh_file %>% #look for highest raw pval that corresponds to .05 Bon
    filter(Bon_p < .05) 
    p_bc <- max(p_bon_cut$P)
  
    cmh_file <- cmh_file %>%
    mutate("Gene_Symbol" = case_when(
      Gene_Symbol == "Protein of unknown function;" ~ NA,
      Gene_Symbol != "Protein of unknown function;" ~ Gene_Symbol
    ))
    
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
        select(CHR, POS, BPcum, SNP, P, Bon_p, Bonp_rank, Gene_Symbol)
    
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
    
    #filter data set to ones I want to label (top 10)
    gene_ID <- cmh_file_cum %>%
      filter(Bonp_rank < 188)
    
    man_p_ggplo <- ggplot(cmh_file_cum, aes(x = BPcum, y = -log(P, base = 10))) + geom_point(aes(color=as.factor(CHR)), alpha=0.3, size=1.3) + scale_color_manual(values = rep(c("#3497A9FF", "#A0DFB9FF"), 8)) + geom_point(data = clump, aes(x = BPcum, y = -log(P, base = 10)), color= "#0B0405FF", alpha=.3 , size=1.3) + geom_text_repel(
      data = gene_ID,
      aes(x = BPcum, y = -log(P, base = 10), label = Gene_Symbol),
      size = 3,
      color = "#0B0405FF",
      max.overlaps = Inf,
      segment.size = 0.3,
      box.padding = 0.4,
      min.segment.length = 0
    ) + scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center, expand = c(0.02, 0.02)) + scale_y_continuous(expand = c(0,0)) + ggtitle( paste("Manhattan Plot with clump sites", name, sep = "")) + labs(x = "Scaffold") + geom_hline(aes(yintercept = -log(p_bc, base = 10)), linetype = "dashed", color="#0B0405FF") + theme_bw() + theme(axis.text = element_text(size = 15), legend.position = "none")
    
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
#dcgm <- fread(file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/FDR_dcgm")
dcgm_clump_ID_rank <- fread(file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/dcgm_.1ID")

dcgm_rank <- dcgm_clump_ID_rank %>%
  arrange(Bon_p) %>%
  mutate(Bonp_rank = dense_rank(Bon_p)) %>%
  separate(Info, into = c("Info","Note"), sep = ";Note=") %>%
  separate(Note, into = c("Gene_Symbol", "Info2"), sep = ":") %>%
  select(CHR, SNP, BP, P, Bon_p, Bonp_rank, Gene_Symbol, Dataset)

#Add gene ID


#Make Manhattan Plots
#man_plot_clump(dcgm, dcgm_clump, "dcgm")
#man_plot_clump_gg(dcgm, clump_loci, "dcgm")
man_plot_clump_gg(dcgm_rank, clump_loci, "dcgm")
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
  ggsave(man_p, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ManhattanP_", name, ".png", sep = ""), device = "png", height = 3, width = 6, units = "in")
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

#repeating for xpehh
xpehh_clumpsites <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/xpehh_sigsites.txt")
xpehh_bed <- xpehh_clumpsites %>%
  mutate("START" = POS - 5000, "STOP" = POS + 5000) %>%
  rowwise() %>%
  arrange(CHROM) %>%
  select(CHROM, START, STOP)

write.table(xpehh_bed, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/dcgm_sig_xpehh.bed", sep = "\t", row.names = F, col.names = F, quote = F)
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
        
        #now looking at the xpehh outlier sites
        xpehh_add <- xpehh_clumpsites %>%
          rowwise() %>%
          mutate("POS_end" = pix_file$BPcum[which(pix_file$window_pos_1 < POS & pix_file$window_pos_2 > POS & pix_file$Scaffold == CHROM)], "FST" = pix_file$FST[which(pix_file$window_pos_1 < POS & pix_file$window_pos_2 > POS & pix_file$Scaffold == CHROM)])
        man_p_fst_xpehh <- ggplot(pix_file, aes(x = BPcum, y = FST)) + geom_point(aes(color=as.factor(Scaffold)), alpha=0.3, size=1.3) + scale_color_manual(values = rep(c("#3497A9FF", "#A0DFB9FF"), 8)) + geom_point(data = xpehh_add, aes(x = POS_end, y = FST), color="#0B0405FF", alpha=.3 , size=1.3) + scale_x_continuous(label = axisdf$Scaffold, breaks= axisdf$center, expand = c(0.02, 0.02)) + scale_y_continuous(expand = c(0,0)) + ggtitle( paste("FST across the genome ", name, sep = "")) + labs(x = "Scaffold") + theme_bw() + theme(axis.text = element_text(size = 15), legend.position = "none") #+ theme(legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank())
        ggsave(man_p_fst_xpehh, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/xpehhhigh_ManhattanP_FST_", name, ".png", sep = ""), device = "png", height = 4, width = 18, units = "in")
        
        #having both highlighted
        man_p_fst_both <- ggplot(pix_file, aes(x = BPcum, y = FST)) + geom_point(aes(color=as.factor(Scaffold)), alpha=0.3, size=1.3) + scale_color_manual(values = rep(c("#3497A9FF", "#A0DFB9FF"), 8)) + geom_point(data = xpehh_add, aes(x = POS_end, y = FST), color="purple", alpha=.3 , size=1.3) + geom_point(data = dcgfclump_add, aes(x = POS, y = FST), color="#0B0405FF", alpha=.3 , size=1.3) + scale_x_continuous(label = axisdf$Scaffold, breaks= axisdf$center, expand = c(0.02, 0.02)) + scale_y_continuous(expand = c(0,0)) + ggtitle( paste("FST across the genome ", name, sep = "")) + labs(x = "Scaffold") + theme_bw() + theme(axis.text = element_text(size = 15), legend.position = "none") #+ theme(legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank())
        ggsave(man_p_fst_both, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/both_ManhattanP_FST_", name, ".png", sep = ""), device = "png", height = 4, width = 18, units = "in")
          
        #ok also looking at the xpehh sites shared across 2 or more pairs
        hitlist <- readRDS("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/sitesshared_perpairxpehh.RDS")

        pairadd <- hitlist %>%
          filter(totalpairs_sigforsite > 1) %>%
          select(locus_id, totalpairs_sigforsite) %>%
          rowwise() %>%
          mutate("CHROM" = as.numeric(str_split(locus_id, ":")[[1]][1]), "POS" = as.numeric(str_split(locus_id, ":")[[1]][2])) 
        pairadd <- pairadd %>%
          mutate("POS_end" = pix_file$BPcum[which(pix_file$window_pos_1 < POS & pix_file$window_pos_2 > POS & pix_file$Scaffold == CHROM)], "FST" = pix_file$FST[which(pix_file$window_pos_1 < POS & pix_file$window_pos_2 > POS & pix_file$Scaffold == CHROM)])
        man_p_fst_perpairxpehh <- ggplot(pix_file, aes(x = BPcum, y = FST)) + geom_point(aes(color=as.factor(Scaffold)), alpha=0.3, size=1.3) + scale_color_manual(values = rep(c("#3497A9FF", "#A0DFB9FF"), 8)) + geom_point(data = pairadd, aes(x = POS_end, y = FST), color="#0B0405FF", alpha=.3 , size=1.3) + scale_x_continuous(label = axisdf$Scaffold, breaks= axisdf$center, expand = c(0.02, 0.02)) + scale_y_continuous(expand = c(0,0)) + ggtitle( paste("FST across the genome ", name, sep = "")) + labs(x = "Scaffold") + theme_bw() + theme(axis.text = element_text(size = 15), legend.position = "none") #+ theme(legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank())
        ggsave(man_p_fst_perpairxpehh, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/perpairxpehh_ManhattanP_FST_", name, ".png", sep = ""), device = "png", height = 4, width = 18, units = "in")
        
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

        man_p_pi_ag_noc <- ggplot(pix_ag, aes(x = BPcum, y = avg_pi)) + geom_point(aes(color=as.factor(Scaffold)), alpha=0.3, size=1.3) + scale_color_manual(values = rep(c("#3497A9FF", "#A0DFB9FF"), 8)) + ggtitle(paste("Pi across the genome ", name, sep = "")) + scale_x_continuous(label = axisdf_ag$Scaffold, breaks= axisdf_ag$center, expand = c(0.02, 0.02)) + scale_y_continuous(expand = c(0,0)) + labs(x = "Scaffold", y = "pi") + theme_bw() + theme(axis.text = element_text(size = 15), legend.position = "none")
        man_p_pi_nat_noc <- ggplot(pix_nat, aes(x = BPcum, y = avg_pi)) + geom_point(aes(color=as.factor(Scaffold)), alpha=0.3, size=1.3) + scale_color_manual(values = rep(c("#3497A9FF", "#A0DFB9FF"), 8)) + scale_x_continuous(label = axisdf_nat$Scaffold, breaks= axisdf_nat$center, expand = c(0.02, 0.02)) + scale_y_continuous(expand = c(0,0)) + ggtitle(paste("Pi across the genome ", name, sep = "")) + labs(x = "Scaffold", y = "pi") + theme_bw() + theme(axis.text = element_text(size = 15), legend.position = "none")
        ggsave(plot_grid(man_p_pi_ag_noc, man_p_pi_nat_noc, labels = c("AG", "NAT"), nrow = 2), file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Noclump_ManhattanP_Pi_", name, ".png", sep = ""), device = "png", height = 4, width = 18, units = "in")
        
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

xpehh_clumpsites <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/xpehh_sigsites.txt")
xpehh_clumpsites <- xpehh_clumpsites %>%
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

#repeat with negative fsts
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
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 15))


ggsave(fst_plot, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/FST_dist_bypair_wneg_50mb.png", device = "png", height = 6, width = 8, units = "in")

pair_means <- fst_all_50 %>%
  filter(!is.na(avg_wc_fst)) %>%
  group_by(pair) %>%
  summarise(mean_fst = mean(avg_wc_fst, na.rm = TRUE), .groups = "drop")

#get overall mean across all windows and pairs
mean(pair_means$mean_fst, na.rm = T)
mean(fst_genomewide$genomewide_fst, na.rm = T)

fst_plot <- fst_genomewide %>%
  ggplot(aes(x = genomewide_fst, fill = T)) +
  geom_density(alpha = 0.25, linewidth = 0.6) + scale_fill_manual(values = c("gray")) +
  labs(x = "Mean FST") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 15))

ggsave(fst_plot, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/FST_genomewideavg_bypair.png", device = "png", height = 6, width = 8, units = "in")

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

write.table(dif, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy_pidif.txt", quote = F, row.names = F)
#just a note- this is the avg pi from the 10 kb window around the clumped site so the window pos 1 is actually 5 kb downstream of the clumped site
```

#Pixy dif in pi plots - MUST RUN PREVIOUS CODE CHUNK BEFORE
#pixy_output.R and Pixy_analysis.rmd both have more (not sure why 2?)
#just a note- this is the avg pi from the 10 kb window around the
clumped site so the window pos 1 is actually 5 kb downstream of the
clumped site

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
      x_overlap == "Overlap" & y_overlap == "Overlap" ~ "Neither",
      x_overlap == "Above" & y_overlap == "Above" ~ "Nat",
      x_overlap == "Below" & y_overlap == "Below" ~ "Ag"
    ))

plot_leg <- exp_nat_agdif %>%
  group_by(Selected_in) %>%
  summarise("sample" = n()) %>%
  mutate("label" = paste(Selected_in, " (n = ", sample, ")", sep = ""))

exp_nat_agdif <- left_join(exp_nat_agdif, plot_leg, by = "Selected_in")

na_all <- ggplot(filter(exp_nat_agdif, is.na(exp_nat_agdif$Selected_in) == F), aes(x = Mean_ag, y = Mean_nat)) + geom_point(data = filter(exp_nat_agdif, exp_nat_agdif$Selected_in != "Neither"), aes(color = label)) + geom_point(data = filter(exp_nat_agdif, exp_nat_agdif$Selected_in == "Neither"), aes(color = label), alpha = .6) + geom_abline(intercept = 0, slope = 1, linetype = 2) + geom_errorbar(data = filter(exp_nat_agdif, exp_nat_agdif$Selected_in != "Neither"), aes(ymin = CIlower_nat, ymax = CIupper_nat, color = label)) + geom_errorbar(data = filter(exp_nat_agdif, exp_nat_agdif$Selected_in != "Neither"), aes(xmin = CIlower_ag, xmax = CIupper_ag, color = label)) + ggtitle("Mean Pi in Ag and Nat by site exp scale \nCIs either both overlap or both don't both overlap") + scale_color_manual(values = c("#60CEACFF","#382A54FF","#395D9CFF")) + xlim(0,.4) + ylim(0,.4) + theme_bw() + theme(axis.text = element_text(size = 15))
#+ geom_errorbar(data = filter(exp_nat_agdif, exp_nat_agdif$Selected_in == "Both"), aes(ymin = CIlower_nat, ymax = CIupper_nat, color = label)) + geom_errorbar(data = filter(exp_nat_agdif, exp_nat_agdif$Selected_in == "Both"), aes(xmin = CIlower_ag, xmax = CIupper_ag, color = label)) 
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/pixy/divplot_expmean_byPair_selected_naall.png", plot = na_all, device = "png", dpi = 300, height = 6, width = 10, units = "in") #I replaced here below w better genomewide perms and estimates

#plot empty plot for proposal 
empty_na_all <- ggplot(filter(exp_nat_agdif, is.na(exp_nat_agdif$Selected_in) == F), aes(x = Mean_ag, y = Mean_nat)) + geom_abline(intercept = 0, slope = 1, linetype = 2) + xlim(0,.4) + ylim(0,.4) + ggtitle("Mean Pi in Ag and Nat by site exp scale \nCIs either both overlap or both don't both overlap") + scale_color_manual(values = c("#60CEACFF","#395D9CFF","#382A54FF")) + theme_bw()

ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/pixy/divplot_empty_naall.png", plot = empty_na_all, device = "png", dpi = 300, height = 6, width = 10, units = "in")
```

#Making distribution of slopes based on permuted nat v ag labels for a
locus #Have to run code chunk with plot pairs function first to get dif
#just a note- this is the avg pi from the 10 kb window around the
clumped site so the window pos 1 is actually 5 kb downstream of the
clumped site

```{r}
dif <- read.table("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy_pidif.txt", header = T)

permed_global <- function(x){
  perm <- dif
  perm_site <- perm %>%
    group_by(windowcenter) %>%
    summarize(
      all_list = list(c(avg_pi_ag, avg_pi_nat)),
      chrom = unique(chrOrder),
      .groups = "drop") %>%
    mutate(
      ag_draw  = map(all_list, ~ sample(.x, 17, replace = F)),
      nat_draw = map2(all_list, ag_draw, ~ .x[!.x %in% .y])
    ) %>%
    rowwise() %>%
    mutate("Mean_nat" = mean(unlist(nat_draw)), "Mean_ag" = mean(unlist(ag_draw)), "Median_nat" = median(unlist(nat_draw)), "Median_ag" = median(unlist(ag_draw)), "sd_nat" = sd(unlist(nat_draw), na.rm = T), "sd_ag" = sd(unlist(ag_draw), na.rm = T), "N_pair" = 17, "se_nat" = sd_nat/sqrt(N_pair), "se_ag" = sd_ag/sqrt(N_pair))
  
  #Get best fit line to compare
  line_fit <- lm(Mean_nat ~ Mean_ag, data = perm_site)
  perm_lines <- c(as.numeric(coef(line_fit)["Mean_ag"]), as.numeric(coef(line_fit)["(Intercept)"]))
  return(perm_lines)
  
  #ggplot(perm_site, aes(x = Mean_ag, y = Mean_nat)) + geom_point() + geom_smooth(method="loess") + geom_abline(intercept = 0, slope = 1, linetype = 2) + ggtitle("Mean Pi in Ag and Nat by site") + theme_bw()
}

results <- replicate(1000, permed_global(x), simplify = T)
rownames(results) <- c("Slope", "Intercept")
results_global <- t(results)

natvagdif <- dif %>%
  group_by(windowcenter) %>%
  summarise("Mean_nat" = mean(avg_pi_nat), "Mean_ag" = mean(avg_pi_ag), "Median_nat" = median(avg_pi_nat), "Median_ag" = median(avg_pi_ag), "sd_nat" = sd(avg_pi_nat, na.rm = T), "sd_ag" = sd(avg_pi_ag, na.rm = T), "N_pair" = n(), "se_nat" = sd_nat/sqrt(N_pair), "se_ag" = sd_ag/sqrt(N_pair), "CHR" = unique(chrOrder))

obs_fit <- lm(Mean_nat ~ Mean_ag, data = natvagdif)
obs_line <- tibble("Slope" = as.numeric(coef(obs_fit)["Mean_ag"]), "Intercept" = as.numeric(coef(obs_fit)["(Intercept)"]))

permed_global_plot <- ggplot() + geom_histogram(data = results_global, aes(x = Slope)) + geom_vline(xintercept = obs_line$Slope, color = "blue") + ggtitle("Slope of best fit line for 1000 replicates of pi means \nfor permuted ag and nat label within any pair") + theme_bw()
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/pixy/divplot_nonpair_permuted_slopes.png", plot = permed_global_plot, device = "png", dpi = 300, height = 3.19, width = 6.93, units = "in")
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
critical_xpehh <- xpehh_file_cum %>%
  filter(norm_xpehh > 5 | norm_xpehh < -5) %>%
  mutate(dataset_x = "XPEHH")
  
overlap <- merge(clump_loci, critical_xpehh, by.x = "SNP", by.y = "locus_id")
#16 overlap - but xpehh hits not LD thinned

#checking non clump loci
nonclump_cmhsig <- read.table("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/P.05FDR_dcgm", header = T)

nonclump_cmhsig <- nonclump_cmhsig %>%
  filter(Bon_p < .05) %>%
  mutate(dataset_c = "CMH") %>%
  rowwise() %>%
  mutate("chr" = as.numeric(str_split(SNP, ":")[[1]][1]))

overlap_nonclump <- merge(nonclump_cmhsig, critical_xpehh, by.x = "SNP", by.y = "locus_id", all = T)
#93 loci overlap

overlap_nonclump <- overlap_nonclump %>%
  mutate("dataset" = case_when(
    dataset_c == "CMH" & dataset_x == "XPEHH" ~ "both",
    is.na(dataset_c) == T & dataset_x == "XPEHH" ~ "XPEHH",
    dataset_c == "CMH" & is.na(dataset_x) == T ~ "CMH"), "CHROMOSOME" = case_when(
      dataset_c == "CMH" ~ chr,
      dataset_x == "XPEHH" ~ CHROM
    ))

xpehh_ref <- overlap_nonclump %>%
  filter(dataset_x == "XPEHH") %>%
  select(CHROMOSOME, POS, BPcum)

overlap_nonclump <- overlap_nonclump %>%
  rowwise() %>%
  mutate(BPcum_all = if (dataset == "XPEHH" | dataset == "both") {
    BPcum
  } else {
      ref <- xpehh_ref[xpehh_ref$CHROMOSOME == CHROMOSOME & xpehh_ref$POS <= BP, ]
    if (nrow(ref) != 0){
      as.numeric(ref$BPcum[which.max(ref$POS)] + (BP - ref$POS[which.max(ref$POS)]))
    } else if (CHROMOSOME > 1) {
      ref <- xpehh_ref[xpehh_ref$CHROMOSOME == CHROMOSOME - 1, ]
      as.numeric(ref$BPcum[which.max(ref$POS)] + BP)
    } else {
      BP
    }}) %>%
  ungroup()

man_p_ggplo <- ggplot() + geom_point(data = overlap_nonclump %>% filter(dataset == "XPEHH"), aes(x = BPcum_all, y = 1, color=as.factor(CHROMOSOME)), alpha=0.3, size=1.3) + geom_point(data = overlap_nonclump %>% filter(dataset == "CMH"), aes(x = BPcum_all, y = 0, color=as.factor(CHROMOSOME)), alpha=0.3, size=1.3) + geom_point(data = overlap_nonclump %>% filter(dataset == "both"), aes(x = BPcum_all, y = .5, color=as.factor(CHROMOSOME)), alpha=0.3, size=1.3) + scale_color_manual(values = rep(c("#3497A9FF", "#A0DFB9FF"), 8)) + scale_x_continuous(label = axisdf$CHROM, breaks= axisdf$center, expand = c(0.02, 0.02)) + ggtitle( "Signficant XPEHH and CMH Loci (not LD thinned)") + labs(x = "Scaffold") + theme_bw() + theme(axis.text = element_text(size = 15), legend.position = "none")
man_p_ggplo
ggsave(man_p_ggplo, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/comp_xpehh_cmh.png", sep = ""), device = "png", height = 6, width = 32, units = "in")

#plot by scaffold to see better
for(s in 1:16){
  man_p_ggplo <- ggplot() + geom_point(data = overlap_nonclump %>% filter(dataset == "XPEHH", CHROMOSOME == s), aes(x = BPcum_all, y = 1, color=as.factor(CHROMOSOME)), alpha=0.3, size=1.3) + geom_point(data = overlap_nonclump %>% filter(dataset == "CMH", CHROMOSOME == s), aes(x = BPcum_all, y = 0, color=as.factor(CHROMOSOME)), alpha=0.3, size=1.3) + geom_point(data = overlap_nonclump %>% filter(dataset == "both", CHROMOSOME == s), aes(x = BPcum_all, y = .5, color=as.factor(CHROMOSOME)), alpha=0.3, size=1.3) + scale_color_manual(values = rep(c("#3497A9FF", "#A0DFB9FF"), 8)) + scale_x_continuous(label = axisdf$CHROM, breaks= axisdf$center, expand = c(0.02, 0.02)) + ggtitle( "Signficant XPEHH and CMH Loci (not LD thinned)") + labs(x = "Scaffold") + theme_bw() + theme(axis.text = element_text(size = 15), legend.position = "none")
ggsave(man_p_ggplo, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/comp_xpehh_cmh", s, ".png", sep = ""), device = "png", height = 6, width = 32, units = "in")
}

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

#Look at how XPEHH and CMH correlates with recomb rate
```{r}
xpehh <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/sel_scan/Aggregated_xpehh.txt", header = T)
cmh <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/FDR_dcgm")
h12_pooled <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/Aggregated_H12scan_all")
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
    rec_rate > bins[4] ~ "5"))

norm_rec <- ggplot(rec_bins_wr) + geom_boxplot(aes(y = norm_xpehh, group = bin)) + xlab("Recombination rate") + labs(title = "Recombination rate binned by quantile \n with normalized xpehh values in that rec rate")

raw_rec <- ggplot(rec_bins_wr) + geom_boxplot(aes(y = raw_xpehh, group = bin)) + xlab("Recombination rate") + labs(title = "Recombination rate binned by quantile \n with raw xpehh values in that rec rate")

ggsave(norm_rec, filename = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/recomb_rate_wnormxpehh.png")
ggsave(raw_rec, filename = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/recomb_rate_wrawxpehh.png")

#repeat for cmh
cmh <- cmh %>%
  rowwise() %>%
  mutate("CHROM" = as.integer(str_split(CHR, "_")[[1]][2])) %>%
  ungroup()

rec_bins_wr_cmh <- cmh %>%
  left_join(
    rec_bins %>% select(chr, win_pos_1, win_pos_2, r),
    join_by(CHROM == chr, BP >= win_pos_1, BP <= win_pos_2)) %>%
  rename(rec_rate = r) %>%
  slice_head(n = 1, by = SNP)

#assign these to their rec rate bin
rec_bins_wr_cmh <- rec_bins_wr_cmh %>%
  mutate("bin" = case_when(
    rec_rate <= bins[1] ~ "1",
    rec_rate > bins[1] & rec_rate <= bins[2] ~ "2",
    rec_rate > bins[2] & rec_rate <= bins[3] ~ "3",
    rec_rate > bins[3] & rec_rate <= bins[4] ~ "4",
    rec_rate > bins[4] ~ "5"))

#site 12:14827 is the only one that has rec bin 1, and it shouldnt so just dropping it bc the rest look good

norm_rec_cmh <- ggplot(rec_bins_wr_cmh) + geom_boxplot(aes(y = FDR_p, group = bin)) + xlab("Recombination rate") + labs(title = "Recombination rate binned by quantile \n with FDR pvalues in that rec rate")

raw_rec_cmh <- ggplot(rec_bins_wr_cmh) + geom_boxplot(aes(y = P, group = bin)) + xlab("Recombination rate") + labs(title = "Recombination rate binned by quantile \n with raw pvalues in that rec rate")

ggsave(norm_rec_cmh, filename = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/recomb_rate_wfdrcmh.png")
ggsave(raw_rec_cmh, filename = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/recomb_rate_wrawcmh.png")

#repeat for h12 pooled ag-nat run
rec_bins_wr_h12 <- h12_pooled %>%
  left_join(
    rec_bins %>% select(chr, win_pos_1, win_pos_2, r),
    join_by(chr == chr, win_center >= win_pos_1, win_center <= win_pos_2)) %>%
  rename(rec_rate = r) %>%
  slice_head(n = 1, by = locusID)

#assign these to their rec rate bin
rec_bins_wr_h12 <- rec_bins_wr_h12 %>%
  mutate("bin" = case_when(
    rec_rate <= bins[1] ~ "1",
    rec_rate > bins[1] & rec_rate <= bins[2] ~ "2",
    rec_rate > bins[2] & rec_rate <= bins[3] ~ "3",
    rec_rate > bins[3] & rec_rate <= bins[4] ~ "4",
    rec_rate > bins[4] ~ "5"))

#site 12:14827 is the only one that has rec bin 1, and it shouldnt so just dropping it bc the rest look good

rec_h12 <- ggplot(rec_bins_wr_h12) + geom_boxplot(aes(y = H12, group = bin)) + xlab("Recombination rate") + labs(title = "Recombination rate binned by quantile \n with H12 values from pooled ag-nat run in that rec rate")

rec_h2h1 <- ggplot(rec_bins_wr_h12) + geom_boxplot(aes(y = H2H1, group = bin)) + xlab("Recombination rate") + labs(title = "Recombination rate binned by quantile \n with H2h1 values from pooled ag-nat run in that rec rate")

ggsave(rec_h12, filename = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/recomb_rate_wh12.png")
ggsave(rec_h2h1, filename = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/recomb_rate_wh2h1.png")

#plot continuous and point

rec_h12 <- ggplot(rec_bins_wr_h12) + geom_point(aes(x = rec_rate, y = H12)) + geom_smooth(method = lm, aes(x = rec_rate, y = H12)) + xlab("Recombination rate") + labs(title = "Recombination rate binned by quantile \n with H12 values from pooled ag-nat run in that rec rate")
ggsave(rec_h12, filename = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/recomb_rate_wh12_pts.png")
rec_h12

rec_h2h1 <- ggplot(rec_bins_wr_h12) + geom_point(aes(x = rec_rate, y = H2H1)) + geom_smooth(method = lm, aes(x = rec_rate, y = H2H1)) + xlab("Recombination rate") + labs(title = "Recombination rate binned by quantile \n with H2h1 values from pooled ag-nat run in that rec rate")
ggsave(rec_h2h1, filename = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/recomb_rate_wh2h1_pts.png")
rec_h2h1

norm_rec_cmh <- ggplot(rec_bins_wr_cmh) + geom_point(aes(y = Bon_p, x = rec_rate)) + geom_smooth(method = lm, aes(y = Bon_p, x = rec_rate)) + xlab("Recombination rate") + labs(title = "Recombination rate binned by quantile \n with Bon pvalues in that rec rate")
norm_rec_cmh
ggsave(norm_rec_cmh, filename = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/recomb_rate_wboncmh_pts.png")

plot_df <- rec_bins_wr[, c("rec_rate", "norm_xpehh")]
norm_rec <- ggplot(plot_df) + geom_point(aes(x = rec_rate, y = norm_xpehh)) + geom_smooth(method = lm, aes(y = norm_xpehh, x = rec_rate)) + xlab("Recombination rate") + labs(title = "Recombination rate binned by quantile \n with normalized xpehh values in that rec rate")
ggsave(norm_rec, filename = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/recomb_rate_wnormxpehh_pts.png")
norm_rec

h12_line <- glm(rec_bins_wr_h12$H12 ~ rec_bins_wr_h12$rec_rate)
h2h1_line <- glm(rec_bins_wr_h12$H2H1 ~ rec_bins_wr_h12$rec_rate)
cmh_line <- glm(rec_bins_wr_cmh$Bon_p ~ rec_bins_wr_cmh$rec_rate)
xpehh_line <- glm(rec_bins_wr$norm_xpehh ~ rec_bins_wr$rec_rate)

summary(h12_line)
summary(h2h1_line)
summary(cmh_line)
summary(xpehh_line)

pseudo_r2_cmh <- 1 - (cmh_line$deviance / cmh_line$null.deviance) #3.14x10^-10
pseudo_r2_xpehh <- 1 - (xpehh_line$deviance / xpehh_line$null.deviance) #5.44x10^-5
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

for(p in 1:17){
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
  mutate("BPcum" = xpehh_file_cum$BPcum[which(xpehh_file_cum$locus_id == SNP)[1]], "norm_xpehh" = xpehh_file_cum$norm_xpehh[which(xpehh_file_cum$locus_id == SNP)[1]])  

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

#comparing selscan pooled to by pair - run on randi bc a lot of data

```{r}
#module load gcc/12.1.0
#module load R (version 4.2.1)
#R
library(data.table)
library(tidyverse)

pooled_run <- fread("/scratch/espolston/Aggregated_xpehh.txt", header = T)
pooled_run <- pooled_run %>%
  filter(norm_xpehh >= 5 | norm_xpehh <= -5) 

summary <- tibble("pair" = numeric(), "num_sigvalsforpair" = numeric(), "num_sitesoverlap" = numeric(), "num_sitesoverlap_wpooled(both_above_5)" = numeric(), "num_sitesoverlap_wpooled(both_below_5)" = numeric(), "num_sitesoverlap_wpooled(signs_dontmatch)" = numeric(), "sites" = list())

for(p in 1:17){
  pair_run <- fread(paste("/scratch/espolston/sel_scan_pair/Aggregated_pair_", p, "_xpehh.txt", sep =""), header = T)
  
  #filter down to just the significant sites, add col of locus_id:norm_xpehh to look at later
  pair_run <- pair_run %>%
    filter(norm_xpehh >= 5 | norm_xpehh <= -5) %>%
    mutate("locus_id:norm_xpehh" = paste(locus_id, ":", norm_xpehh, sep = ""))
  
  #get which sites overlap w pooled run
  overlaps <- pair_run[which(pair_run$locus_id %in% pooled_run$locus_id),]
  
  #see how sign of these compare
  overlaps <- overlaps %>%
    rowwise() %>%
    mutate("sign_match" = case_when(
      norm_xpehh >= 5 & pooled_run$norm_xpehh[which(pooled_run$locus_id == locus_id)] >= 5 ~ "both_pos",
      norm_xpehh <= -5 & pooled_run$norm_xpehh[which(pooled_run$locus_id == locus_id)] <= -5 ~ "both_neg",
      (norm_xpehh >= 5 & pooled_run$norm_xpehh[which(pooled_run$locus_id == locus_id)] <= -5) | (norm_xpehh <= -5 & pooled_run$norm_xpehh[which(pooled_run$locus_id == locus_id)] >= 5) ~ "nomatch"))
  
  #add to summary table of all 17 pairwise comparisons
  summary <- add_row(summary, "pair" = p, "num_sigvalsforpair" = length(pair_run$locus_id), "num_sitesoverlap" = length(overlaps$locus_id), "num_sitesoverlap_wpooled(both_above_5)" = sum(overlaps$sign_match == "both_pos"), "num_sitesoverlap_wpooled(both_below_5)" = sum(overlaps$sign_match == "both_neg"), "num_sitesoverlap_wpooled(signs_dontmatch)" = sum(overlaps$sign_match == "nomatch"), "sites" = list(pair_run$'locus_id:norm_xpehh'))
}

saveRDS(summary, "/scratch/espolston/sel_scan_pair/summary_pooledvpaired.txt")
```

#finishing analyzing the xpehh pool vs pairs local

```{r}
poolvpair <- readRDS("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/summary_pooledvpaired.rds")

nums_poolvpair <- poolvpair[,c(1,2,3,4,5,6)]
colnames(nums_poolvpair) <- c("pair", "num_sig", "num_overlap_total", "num_overlap_bothpos", "num_overlap_bothneg", "num_overlap_sign_nomatch")

#add in sample size in ag + nat
size <- read.csv("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pop_samplesize.csv")

size <- size %>%
  group_by(Pair.) %>%
  summarise("popsize" = sum(Sample.size.after.IBD.and.other.spp.drop))

#fixing pop 7b
size$popsize[15] <- size$popsize[15] + size$popsize[16]
size <- size[-16,]
size$pair <- as.numeric(size$'Pair.')

nums_poolvpair <- merge(nums_poolvpair, size)
nums_poolvpair <- select(nums_poolvpair, pair, num_sig, popsize, num_overlap_total, num_overlap_bothpos, num_overlap_bothneg, num_overlap_sign_nomatch)
```

#looking at selscan pair hits to see if in same regions, just not exact same hit
```{r}
#ran xpehh_pairsoverlap.sbatch on randi
window_xpehh_allpairs <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_xpehh_allpairs.txt")

#confirming worked correctly
num_hits_perpair <- window_xpehh_allpairs %>%
  group_by(pair) %>%
  summarise("total" = sum(n_sig_sites)) %>%
  ungroup()
#all correct except for pair 2 where missing 10 sites but I'm gonna ignore bc all the rest correct?

window_xpehh_allpairs$snp <- paste(window_xpehh_allpairs$chrom, ":", window_xpehh_allpairs$win_start, sep = "")

overlap_sites <- window_xpehh_allpairs %>%
  filter(n_sig_sites > 0) %>%
  group_by(snp) %>%
  summarize("pairs_sig" = list(pair), "nsig_perpair" = list(n_sig_sites), "num_pairoverlap" = n()) %>%
  ungroup() 

length(filter(overlap_sites, overlap_sites$num_pairoverlap > 1)$num_pairoverlap)
#only 1281 windows overlap for 2 or more pairs


#looking at the number of these that overlap with the pooled hits
xpehh <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_xpehh_all.txt")
colnames(xpehh) <- c("chrom", "win_start", "win_end", "n_sig_sites")
xpehh$snp <- paste(xpehh$chrom, ":", xpehh$win_start, sep = "")
xpehh$pair <- "pool"

pool <- rbind(window_xpehh_allpairs, xpehh)

overlap_pool <- pool %>%
  filter(n_sig_sites > 0) %>%
  group_by(snp) %>%
  summarize("pairs_sig" = list(pair), "nsig_perpair" = list(n_sig_sites), "num_pairoverlap" = n(), "pool?" = "pool" %in% pair) %>%
  ungroup() 

overlap_pool_display <- overlap_pool %>%
  mutate(pairs_sig_str = map_chr(pairs_sig, ~ paste(.x, collapse = ", "))) %>%
  select(snp, pairs_sig_str, nsig_perpair, num_pairoverlap, "pool?")

length(filter(overlap_pool, overlap_pool$num_pairoverlap > 1)$num_pairoverlap)
length(filter(overlap_pool, overlap_pool$num_pairoverlap > 1 & overlap_pool$"pool?" == T)$num_pairoverlap)
```

#look at H12 overlaps - after controlling for recomb rate (dropping
everything in first quartile of recomb rate) - not using this

```{r}
#have thinned out low recomb rate regions
#also check if these peaks are in areas of possible SVs or missing data- do I already control for missing data enough??
H12Scan_all <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/dcgm_hap_stats_Scaffold_", 1,"_recombfiltered.txt", sep = ""))
colnames(H12Scan_all) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123")
peaks_all <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/dcgm_hap_peaks_Scaffold_", 1,"_recombfiltered.txt", sep = ""))
colnames(peaks_all) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123","smallest_edge_coordofpeak", "largest_edge_coordofpeak")

#get the top 10 peaks for each chrom
#peaks_all$top_peak <- c(rep(T, 10), rep(F, length(peaks_all$win_center) - 10))

peaks_all$chr <- rep(1, length(peaks_all$win_center))
H12Scan_all$chr <- rep(1, length(H12Scan_all$win_center))
H12Scan_all$bpcum <- H12Scan_all$win_center
peaks_all$bpcum <- peaks_all$win_center

for(s in 2:16){
  H12Scan_all_in = read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/dcgm_hap_stats_Scaffold_", s,"_recombfiltered.txt", sep = ""))
  peaks_all_in = read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/dcgm_hap_peaks_Scaffold_", s,"_recombfiltered.txt", sep = ""))
  colnames(H12Scan_all_in) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123")
  colnames(peaks_all_in) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123","smallest_edge_coordofpeak", "largest_edge_coordofpeak")
  
  #peaks_all_in$top_peak <- c(rep(T, 10), rep(F, length(peaks_all_in$win_center) - 10))
  peaks_all_in$chr <- rep(s, length(peaks_all_in$win_center))
  H12Scan_all_in$chr <- rep(s, length(H12Scan_all_in$win_center))
  
  #make bpcum
  prelength <- max(H12Scan_all$bpcum)
  peaks_all_in$bpcum <- peaks_all_in$win_center + prelength
  H12Scan_all_in$bpcum <- H12Scan_all_in$win_center + prelength
  
  H12Scan_all <- rbind(H12Scan_all, H12Scan_all_in)
  peaks_all <- rbind(peaks_all, peaks_all_in)
}

write.table(H12Scan_all, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/Aggregated_H12scan_all", quote = F, row.names = F)
write.table(peaks_all, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/Aggregated_peaks_all", quote = F, row.names = F)

#get top 50 h12 peaks genomewide
H12Scan_all <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/Aggregated_H12scan_all")

H12Scan_all <- H12Scan_all %>%
  arrange(desc(H12))

H12Scan_all$top50 <- c(seq(1, 50, by = 1), rep(NA, length(H12Scan_all$win_center)-50))
H12Scan_all <- H12Scan_all %>%
  mutate("locusID" = paste(chr, ":", win_center, sep = ""))

#Pull out these top 50 for ag and nat run
H12Scan_ag <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/dcgm_hap_stats_Scaffold_", 1,"_recombfiltered_ag.txt", sep = ""))
colnames(H12Scan_ag) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123")
H12Scan_ag <- H12Scan_ag %>%
  mutate("chr" = rep(1, length(win_center)), "bpcum" = win_center, "locusID" = paste(chr, ":", win_center, sep = ""))

for(s in 2:16){
  H12Scan_ag_in = read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/dcgm_hap_stats_Scaffold_", s,"_recombfiltered_ag.txt", sep = ""))
  colnames(H12Scan_ag_in) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123")
  
  #make bpcum
  prelength <- max(H12Scan_ag$bpcum)
  H12Scan_ag_in$bpcum <- H12Scan_ag_in$win_center + prelength
  
  H12Scan_ag_in <- H12Scan_ag_in %>%
  mutate("chr" = rep(s, length(win_center)), "locusID" = paste(chr, ":", win_center, sep = ""))
  
  H12Scan_ag <- rbind(H12Scan_ag, H12Scan_ag_in)
}

H12Scan_nat <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/dcgm_hap_stats_Scaffold_", 1,"_recombfiltered_nat.txt", sep = ""))
colnames(H12Scan_nat) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123")
H12Scan_nat <- H12Scan_nat %>%
  mutate("chr" = rep(1, length(win_center)), "bpcum" = win_center, "locusID" = paste(chr, ":", win_center, sep = ""))

for(s in 2:16){
  H12Scan_nat_in = read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/dcgm_hap_stats_Scaffold_", s,"_recombfiltered_nat.txt", sep = ""))
  colnames(H12Scan_nat_in) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123")
  
  #make bpcum
  prelength <- max(H12Scan_nat$bpcum)
  H12Scan_nat_in$bpcum <- H12Scan_nat_in$win_center + prelength
  
  H12Scan_nat_in <- H12Scan_nat_in %>%
  mutate("chr" = rep(s, length(win_center)), "locusID" = paste(chr, ":", win_center, sep = ""))
  
  H12Scan_nat <- rbind(H12Scan_nat, H12Scan_nat_in)
}

#filter to top 50 h12 vals genome wide, then pull out the H2H1 vals in ag and nat and then plot these
toph12all_agnatvals <- H12Scan_all %>%
  filter(is.na(top50) == F) %>%
  rowwise() %>%
  mutate(H12_ag   = H12Scan_ag$H12[match(locusID, H12Scan_ag$locusID)],
    H12_nat  = H12Scan_nat$H12[match(locusID, H12Scan_nat$locusID)],
    H2H1_ag  = H12Scan_ag$H2H1[match(locusID, H12Scan_ag$locusID)],
    H2H1_nat = H12Scan_nat$H2H1[match(locusID, H12Scan_nat$locusID)])

#now look at H2/H1 (so H12 is a sweep, H2/H1 is the softness of it)
h2h1plot <- ggplot() + geom_point(data = toph12all_agnatvals, aes(x = H2H1_ag, y = H2H1_nat))  + ggtitle("Values of H2H1 in ag and nat\nfor top H12 values in pooled run") + geom_abline(intercept = 0, slope = 1, linetype = 2) + theme_bw() + theme(axis.text = element_text(size = 15)) + scale_color_manual(values = c("#60CEACFF","#395D9CFF","#382A54FF","gray80")) + ylim(0,.65) + xlim(0,.65)
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/h2h1_agvnat_fortop50.png", plot = h2h1plot, device = "png", dpi = 300, height = 6, width = 10, units = "in")


#Maybe above plot isn’t the right thing to ask- distribution of top 50 h2h1 vals in pooled, distribution of h2h1 vals for cmh/xpehh sites
top50pool <- filter(H12Scan_all, is.na(top50) == F)

cmh_clumpsites <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clumphits_dcgm.txt")
cmh_clumpsites$snp <- paste(cmh_clumpsites$CHROM, ":", cmh_clumpsites$BP, sep ="")

xpehh_clumpsites <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/xpehh_sigsites.txt")

#get sites that the locus ID matches (the cmh site is the window center)
cmh_ag <- H12Scan_ag %>%
  filter(locusID %in% cmh_clumpsites$snp)

cmh_nat <- H12Scan_nat %>%
  filter(locusID %in% cmh_clumpsites$snp)

cmh_pool <- H12Scan_all %>%
  filter(locusID %in% cmh_clumpsites$snp)

xpehh_ag <- H12Scan_ag %>%
  filter(locusID %in% xpehh_clumpsites$locus_id)

xpehh_nat <- H12Scan_nat %>%
  filter(locusID %in% xpehh_clumpsites$locus_id)

xpehh_pool <- H12Scan_all %>%
  filter(locusID %in% xpehh_clumpsites$locus_id)

#for the remaining sites, take top h12 value from the windows that contain the snp
for(site in 1:length(cmh_clumpsites$snp)){
  if(cmh_clumpsites$snp[site] %in% cmh_ag$locusID) {
  }else{
    row_ag <- H12Scan_ag[which(H12Scan_ag$chr == cmh_clumpsites$CHROM[site] & cmh_clumpsites$BP[site] >= H12Scan_ag$win_pos_1 & cmh_clumpsites$BP[site] <= H12Scan_ag$win_pos_2),] 
    row_ag <- arrange(row_ag, desc(H12))
    cmh_ag <- add_row(cmh_ag, row_ag[1,])
    
    row_nat <- H12Scan_nat[which(H12Scan_nat$chr == cmh_clumpsites$CHROM[site] & cmh_clumpsites$BP[site] >= H12Scan_nat$win_pos_1 & cmh_clumpsites$BP[site] <= H12Scan_nat$win_pos_2),] 
    row_nat <- arrange(row_nat, desc(H12))
    cmh_nat <- add_row(cmh_nat, row_nat[1,])
    
    row_pool <- H12Scan_all[which(H12Scan_all$chr == cmh_clumpsites$CHROM[site] & cmh_clumpsites$BP[site] >= H12Scan_all$win_pos_1 & cmh_clumpsites$BP[site] <= H12Scan_all$win_pos_2),] 
    row_pool <- arrange(row_pool, desc(H12))
    cmh_pool <- add_row(cmh_pool, row_pool[1,])
  }
}

for(site in 1:length(xpehh_clumpsites$locus_id)){
  if(xpehh_clumpsites$locus_id[site] %in% xpehh_ag$locusID) {
  }else{
    row_ag <- H12Scan_ag[which(H12Scan_ag$chr == xpehh_clumpsites$CHROM[site] & xpehh_clumpsites$POS[site] >= H12Scan_ag$win_pos_1 & xpehh_clumpsites$POS[site] <= H12Scan_ag$win_pos_2),] 
    row_ag <- arrange(row_ag, desc(H12))
    xpehh_ag <- add_row(xpehh_ag, row_ag[1,])
    
    row_nat <- H12Scan_nat[which(H12Scan_nat$chr == xpehh_clumpsites$CHROM[site] & xpehh_clumpsites$POS[site] >= H12Scan_nat$win_pos_1 & xpehh_clumpsites$POS[site] <= H12Scan_nat$win_pos_2),] 
    row_nat <- arrange(row_nat, desc(H12))
    xpehh_nat <- add_row(xpehh_nat, row_nat[1,])
    
    row_pool <- H12Scan_all[which(H12Scan_all$chr == xpehh_clumpsites$CHROM[site] & xpehh_clumpsites$POS[site] >= H12Scan_all$win_pos_1 & xpehh_clumpsites$POS[site] <= H12Scan_all$win_pos_2),] 
    row_pool <- arrange(row_pool, desc(H12))
    xpehh_pool <- add_row(xpehh_pool, row_pool[1,])
  }
}

#a lot of the cmh and xpehh sites are in low recombination regions so not showing up in the table
#get num of sites shown
length(xpehh_ag$win_center) - sum(is.na(xpehh_ag$win_center)) #830 sites
length(xpehh_nat$win_center) - sum(is.na(xpehh_nat$win_center)) #830 sites

length(cmh_ag$win_center) - sum(is.na(cmh_ag$win_center)) #729
length(cmh_nat$win_center) - sum(is.na(cmh_nat$win_center)) #729

length(cmh_pool$win_center) - sum(is.na(cmh_pool$win_center)) #729
length(xpehh_pool$win_center) - sum(is.na(xpehh_pool$win_center)) #830

#H2h1_dist <- ggplot() + geom_histogram(data = top50pool, aes(x = H2H1, fill = "top50H12"), alpha = 0.5, position = "identity") + geom_histogram(data = cmh_ag, aes(x = H2H1, fill = "CMH_ag"), alpha = 0.5, position = "identity") + geom_histogram(data = cmh_nat, aes(x = H2H1, fill = "CMH_nat"), alpha = 0.5, position = "identity") + geom_histogram(data = xpehh_ag, aes(x = H2H1, fill = "XPEHH_ag"), alpha = 0.5, position = "identity") + geom_histogram(data = xpehh_nat, aes(x = H2H1, fill = "XPEHH_nat"), alpha = 0.5, position = "identity") + labs(title = "H2H1 for: top 50 h12 vals in pooled set, cmh/xpehh loci in ag and nat")

H2h1_dist <- ggplot() + geom_density(data = top50pool, aes(x = H2H1, fill = "top50H12"), alpha = 0.5, position = "identity") + geom_density(data = cmh_ag, aes(x = H2H1, fill = "CMH_ag"), alpha = 0.5, position = "identity") + geom_density(data = cmh_nat, aes(x = H2H1, fill = "CMH_nat"), alpha = 0.5, position = "identity") + geom_density(data = xpehh_ag, aes(x = H2H1, fill = "XPEHH_ag"), alpha = 0.5, position = "identity") + geom_density(data = xpehh_nat, aes(x = H2H1, fill = "XPEHH_nat"), alpha = 0.5, position = "identity") + geom_density(data = cmh_pool, aes(x = H2H1, fill = "CMH_pool"), alpha = 0.5, position = "identity") + geom_density(data = xpehh_pool, aes(x = H2H1, fill = "XPEHH_pool"), alpha = 0.5, position = "identity") + labs(title = "H2H1 for: top 50 h12 vals in pooled set, cmh/xpehh loci in ag and nat")

ggsave(plot = H2h1_dist, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/h2h1_dist.png", device = "png", dpi = 300, height = 6, width = 10, units = "in")

h12_dist <- ggplot() + geom_histogram(data = top50pool, aes(x = H12, fill = "top50H12"), alpha = 0.5, position = "identity") + geom_histogram(data = cmh_ag, aes(x = H12, fill = "CMH_ag"), alpha = 0.5, position = "identity") + geom_histogram(data = cmh_nat, aes(x = H12, fill = "CMH_nat"), alpha = 0.5, position = "identity") + geom_histogram(data = xpehh_ag, aes(x = H12, fill = "XPEHH_ag"), alpha = 0.5, position = "identity") + geom_histogram(data = xpehh_nat, aes(x = H12, fill = "XPEHH_nat"), alpha = 0.5, position = "identity") + geom_histogram(data = cmh_pool, aes(x = H12, fill = "CMH_pool"), alpha = 0.5, position = "identity") + geom_histogram(data = xpehh_pool, aes(x = H12, fill = "XPEHH_pool"), alpha = 0.5, position = "identity") + labs(title = "H12 for: top 50 h12 vals in pooled set, cmh/xpehh loci in ag and nat")

ggsave(plot = h12_dist, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/h12_dist.png", device = "png", dpi = 300, height = 6, width = 10, units = "in")

```

#Now looking at ag and nat separate h12 runs - not using this

```{r}
#load in ag run
peaks_ag <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/dcgm_hap_peaks_Scaffold_", 1,"_recombfiltered_ag.txt", sep = ""))
H12Scan_ag <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/dcgm_hap_stats_Scaffold_", 1,"_recombfiltered_ag.txt", sep = ""))
colnames(H12Scan_ag) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123")
colnames(peaks_ag) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123","smallest_edge_coordofpeak", "largest_edge_coordofpeak")

#get the top 10 peaks for each chrom
peaks_ag$top_peak <- c(rep(T, 1))
peaks_ag$chr <- rep(1, length(peaks_ag$win_center))
H12Scan_ag$chr <- rep(1, length(H12Scan_ag$win_center))
H12Scan_ag$bpcum <- H12Scan_ag$win_center
peaks_ag$bpcum <- peaks_ag$win_center

for(s in 2:16){
  H12Scan_ag_in = read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/dcgm_hap_stats_Scaffold_", s,"_recombfiltered_ag.txt", sep = ""))
  peaks_ag_in = read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/dcgm_hap_peaks_Scaffold_", s,"_recombfiltered_ag.txt", sep = ""))
  colnames(H12Scan_ag_in) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123")
  colnames(peaks_ag_in) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123","smallest_edge_coordofpeak", "largest_edge_coordofpeak")
  
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
peaks_nat <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/dcgm_hap_peaks_Scaffold_", 1,"_recombfiltered_nat.txt", sep = ""))
H12Scan_nat <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/dcgm_hap_stats_Scaffold_", 1,"_recombfiltered_nat.txt", sep = ""))
colnames(H12Scan_nat) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123")
colnames(peaks_nat) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123","smallest_edge_coordofpeak", "largest_edge_coordofpeak")

#get the top 10 peaks for each chrom
peaks_nat$top_peak <- c(rep(T, 10), rep(F, length(peaks_nat$win_center) - 10))
peaks_nat$chr <- rep(1, length(peaks_nat$win_center))
H12Scan_nat$chr <- rep(1, length(H12Scan_nat$win_center))
H12Scan_nat$bpcum <- H12Scan_nat$win_center
peaks_nat$bpcum <- peaks_nat$win_center

for(s in 2:16){
  H12Scan_nat_in = read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/dcgm_hap_stats_Scaffold_", s,"_recombfiltered_nat.txt", sep = ""))
  peaks_nat_in = read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/dcgm_hap_peaks_Scaffold_", s,"_recombfiltered_nat.txt", sep = ""))
  colnames(H12Scan_nat_in) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123")
  colnames(peaks_nat_in) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123","smallest_edge_coordofpeak", "largest_edge_coordofpeak")
  
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

write.table(H12Scan_ag, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/Aggregated_H12scan_ag", quote = F, row.names = F)
write.table(H12Scan_nat, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/Aggregated_H12scan_nat", quote = F, row.names = F)
write.table(peaks_ag, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/Aggregated_peaks_ag", quote = F, row.names = F)
write.table(peaks_nat, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/recombfilter_H12/Aggregated_peaks_nat", quote = F, row.names = F)

#now compare top 50 peaks - this is based on top H12 value
#add snp
peaks_ag$snp <- paste(peaks_ag$chr, peaks_ag$win_center, sep = ":")
peaks_all$snp <- paste(peaks_all$chr, peaks_all$win_center, sep = ":")
peaks_nat$snp <- paste(peaks_nat$chr, peaks_nat$win_center, sep = ":")

peaks_agtop <- peaks_ag %>%
  filter(top_peak == T) %>%
  arrange(desc(H12))
peaks_agtop <- peaks_agtop[1:50,]

peaks_nattop <- peaks_nat %>%
  filter(top_peak == T) %>%
  arrange(desc(H12))
peaks_nattop <- peaks_nattop[1:50,]

peaks_alltop <- peaks_all %>%
  filter(top_peak == T) %>%
  arrange(desc(H12))
peaks_alltop <- peaks_alltop[1:50,]

ag_all <- which(peaks_agtop$snp %in% peaks_alltop$snp)
nat_all <- which(peaks_nattop$snp %in% peaks_alltop$snp)
ag_nat <- which(peaks_agtop$snp %in% peaks_nattop$snp)
ag_nat_all <- which(peaks_agtop$snp[ag_nat] %in% peaks_alltop$snp)

#investigating ag+nat peaks (top 10 for each chr)
ag_vals <- peaks_ag %>%
  mutate("dataset" = "ag", "H1_ag" = H1, "H2_ag" = H2, "H12_ag" = H12, "H2H1_ag" = H2H1, "H123_ag" = H123, "top_peak_ag" = top_peak)  %>%
  select(snp, dataset, top_peak_ag, H1_ag, H2_ag, H12_ag, H2H1_ag, H123_ag)

nat_vals <- peaks_nat %>%
  mutate("dataset" = "nat", "H1_nat" = H1, "H2_nat" = H2, "H12_nat" = H12, "H2H1_nat" = H2H1, "H123_nat" = H123, "top_peak_nat" = top_peak)  %>%
  select(snp, dataset, top_peak_nat, H1_nat, H2_nat, H12_nat, H2H1_nat, H123_nat)

ag_nat_vals <- merge(ag_vals, nat_vals, by = "snp")

ag_nat_vals <- ag_nat_vals  %>%
  mutate("peak_type" = case_when(
    top_peak_ag == T & top_peak_nat == T ~ "Both",
    top_peak_ag == T & top_peak_nat == F ~ "Ag", 
    top_peak_ag == F & top_peak_nat == T ~ "Nat"
  ))

ag_nat_vals$peak_type <- as.factor(ag_nat_vals$peak_type)

h12plot <- ggplot() + geom_point(data = ag_nat_vals %>% filter(is.na(peak_type) == F), aes(x = log(H12_ag), y = log(H12_nat), color = peak_type)) + geom_point(data = ag_nat_vals %>% filter(is.na(peak_type) == T), aes(x = log(H12_ag), y = log(H12_nat), color = peak_type, alpha = .5)) + ggtitle("Values of H12 in ag and nat for window centers that are in both data sets\nAg-nat peak when one of top 10 on each chrom in ag and same in nat") + geom_abline(intercept = 0, slope = 1, linetype = 2) + theme_bw() + theme(axis.text = element_text(size = 15)) + scale_color_manual(values = c("#60CEACFF","#395D9CFF","#382A54FF","gray80"))
h2h1plot <- ggplot(ag_nat_vals, aes(x = H2H1_ag, y = H2H1_nat)) + geom_point(aes(color = peak_type)) + ggtitle("Values of H2/H1 in ag and nat for window centers that are in both data sets\nAg-nat peak when one of top 10 on each chrom in ag and same in nat") + geom_abline(intercept = 0, slope = 1, linetype = 2) + theme_bw()
#+ scale_color_manual(values = c("#60CEACFF","#395D9CFF","#382A54FF"))

ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/h12_agvnat.png", plot = h12plot, device = "png", dpi = 300, height = 6, width = 10, units = "in")
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/h2h1_agvnat.png", plot = h2h1plot, device = "png", dpi = 300, height = 6, width = 10, units = "in")

#look at distribution of H12, h2h1 values in both habitats
h12ag <- ggplot() + geom_boxplot(data = peaks_ag, aes(x = H12, y = 1, color = "All ag")) + geom_boxplot(data = peaks_ag %>% filter(top_peak == T), aes(x = H12, y = 0, color = "Peak ag")) + ggtitle("Boxplot of Ag peaks H12 scores") + xlim(0, .2)
h12nat <- ggplot() + geom_boxplot(data = peaks_nat, aes(x = H12, y = 1, color = "All nat")) + geom_boxplot(data = peaks_nat %>% filter(top_peak == T), aes(x = H12, y = 0, color = "Peak nat")) + ggtitle("Boxplot of Nat peaks H12 scores") + xlim(0, .2)
h2h1ag <- ggplot() + geom_boxplot(data = peaks_ag, aes(x = H2H1, y = 1, color = "All ag")) + geom_boxplot(data = peaks_ag %>% filter(top_peak == T), aes(x = H2H1, y = 0, color = "Peak ag")) + ggtitle("Boxplot of Ag peaks H2H1 scores") + xlim(0, 1)
h2h1nat <- ggplot() + geom_boxplot(data = peaks_nat, aes(x = H2H1, y = 1, color = "All nat")) + geom_boxplot(data = peaks_nat %>% filter(top_peak == T), aes(x = H2H1, y = 0, color = "Peak nat")) + ggtitle("Boxplot of Nat peaks H2H1 scores") + xlim(0, 1)

ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/h12_agbox.png", plot = h12ag, device = "png", dpi = 300, height = 6, width = 10, units = "in")
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/h12_natbox.png", plot = h12nat, device = "png", dpi = 300, height = 6, width = 10, units = "in")
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/h2h1_agbox.png", plot = h2h1ag, device = "png", dpi = 300, height = 6, width = 10, units = "in")
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/h2h1_natbox.png", plot = h2h1nat, device = "png", dpi = 300, height = 6, width = 10, units = "in")
```

#LD decay window size 200 kb

```{r}
#from https://speciationgenomics.github.io/ld_decay/
library(tidyverse)
# read in data
#ld_bins.5 <- read_tsv("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm_lddecay_.5.ld_decay_bins")

ld_bins <- read_tsv("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/dcgm_lddecay.ld_decay_bins")

# plot LD decay
ld_10 <- ggplot(ld_bins, aes(distance, avg_R2)) + geom_line() +
  xlab("Distance (bp)") + ylab(expression(italic(r)^2)) + labs(title = "LD decay in 200kb windows for thinned 10% data")
#ld_50 <- ggplot(ld_bins.5, aes(distance, avg_R2)) + geom_line() +
  xlab("Distance (bp)") + ylab(expression(italic(r)^2)) + labs(title = "LD decay in 200kb windows for thinned 50% data")

ggsave(ld_10, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/LD_decayplot", device = "png", dpi = 300)

#10 kb is good based on Rozenn's analyses
```

#Comparing XPEHH and CMH

```{r}
#Now need to clump XPEHH and add those loci to the cmh clump loci
#creating clump file- set up for p values so using abs of norm_xpehh
xpehh_in <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/sel_scan/Aggregated_xpehh.txt")

plink_xpehh <- xpehh_in %>%
  mutate(zscore = (norm_xpehh - mean(norm_xpehh))/sd(norm_xpehh)) %>%
  select(CHROM, POS, locus_id, norm_xpehh, zscore)
  
plink_xpehh$pval <-2*pnorm(-abs(plink_xpehh$norm_xpehh))

write.table(plink_xpehh, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/xpehh_plink", row.names = F, quote = F)

#plink_xpehh <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/xpehh_plink")

#cutoffs - none == 5 so had to get min norm xpehh that was sig
min(abs(plink_xpehh$norm_xpehh[which(abs(plink_xpehh$norm_xpehh) >= 5)]))
max(plink_xpehh$pval[which(abs(plink_xpehh$norm_xpehh) == 5.00001)]) 
#stringent cutoff = 5.732734e-07 = .0000005732734
min(abs(plink_xpehh$norm_xpehh[which(abs(plink_xpehh$norm_xpehh) >= 2)]))
max(plink_xpehh$pval[which(abs(plink_xpehh$norm_xpehh) == 2)]) 
#less stringent cutoff = 0.04550026

#ran clumping in xpehh_cmh_overlap.sh

#---analysis of overlap results-----------
xpehh <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_xpehh_all.txt")
colnames(xpehh) <- c("CHROM", "POS1", "POS2", "NUM_XPEHH_CRIT")
cmh <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_cmh_all.txt")
colnames(cmh) <- c("CHROM", "POS1", "POS2", "NUM_CMH_CRIT")

hits <- merge(cmh, xpehh)
hits <- hits %>%
  filter(NUM_CMH_CRIT > 0 | NUM_XPEHH_CRIT > 0) #1432

overlaps <- filter(hits, NUM_CMH_CRIT > 0 & NUM_XPEHH_CRIT > 0) #95
cmhonly <- filter(hits, NUM_CMH_CRIT > 0 & NUM_XPEHH_CRIT == 0) #1043
xpehhonly <- filter(hits, NUM_CMH_CRIT == 0 & NUM_XPEHH_CRIT > 0) #294

#check for overlap w clump cmh loci
cmhclump <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_clumpcmh_all.txt")
colnames(cmhclump) <- c("CHROM", "POS1", "POS2", "NUM_CMH_CLUMP")
hits_clump <- merge(hits, cmhclump)
hits_clump <- hits_clump %>%
  filter(NUM_CMH_CRIT > 0 | NUM_XPEHH_CRIT > 0) #1432 like it should be

clump <- filter(hits_clump, NUM_CMH_CLUMP > 0) #720
clumponly <- filter(hits_clump, NUM_CMH_CLUMP > 0 & NUM_XPEHH_CRIT == 0) #655

xpehhclump <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_clumpxpehh_all.txt")
colnames(xpehhclump) <- c("CHROM", "POS1", "POS2", "NUM_XPEHH_CLUMP")
xhits_clump <- merge(hits, xpehhclump)
xhits_clump <- xhits_clump %>%
  filter(NUM_CMH_CRIT > 0 | NUM_XPEHH_CRIT > 0) #1432 like it should be

xclump <- filter(xhits_clump, NUM_XPEHH_CLUMP > 0) #335
xclumponly <- filter(xhits_clump, NUM_XPEHH_CLUMP > 0 & NUM_CMH_CRIT == 0) #264 windows unique to XPEHH

#merging clumps
allclump <- merge(cmhclump, xpehhclump)
allclump <- filter(allclump, NUM_XPEHH_CLUMP > 0 | NUM_CMH_CLUMP > 0) #1008 windows together

#Now in 264 only XPEHH windows, I pulled out the snps in these and then added these to the 865 CMH clump sites
write.table(xclumponly, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/onlyxpehhclump_windows.txt", row.names = F, col.names = F, sep ="\t", quote = F)

xpehh_sitesinwin <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/onlyxpehhclump_sites.txt")
colnames(xpehh_sitesinwin) <- colnames(xpehhclump) <- c("CHROM", "POS1", "POS2", "XPEHH_CLUMP_SITES")
xpehh_sites <- xpehh_sitesinwin %>%
  rowwise() %>%
  mutate(CHR = strsplit(CHROM, "_")[[1]][2]) %>%
  mutate(snps = paste(paste0(CHR, ":", unlist(strsplit(XPEHH_CLUMP_SITES, ","))), collapse = ", "))

#so these are Xpehh sig sites that are in 10kb windows not tagged by CMH
toadd_xpehh_sites <- unlist(strsplit(paste(xpehh_sites$snps, collapse = ", "), ", "))
xpehh_wval <- filter(plink_xpehh, locus_id %in% toadd_xpehh_sites)
write.table(xpehh_wval, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/xpehh_sigsites.txt", quote = F, row.names = F)

#looking to see if 100kb win in clumping fucks with 10kb win above
xpehh_clumpall <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/clumped_dcgm_xpehh.clumped")
xpehh_clumpall <- xpehh_clumpall %>%
  rowwise() %>%
  mutate("snps"= list((unlist(strsplit(SP2,",")[[1]])))) %>%
  mutate("firstsnp" = snps[1], "lastsnp" = snps[length(snps)])
#looks like at least several are 100kb apart

#now looking at parallel/anti parallel in AFVaper
```

#doing diversity permutations for the AF freq matched windows across the
genome- in case nat pi is just globally elevated

```{r}
#load in pixy results for 10kb windows
pi_allpair_10kbwin <- tibble("pop" = character(), "pair" = numeric(), "chromosome" = character(), "window_pos_1" = numeric(), "window_pos_2" = numeric(), "avg_pi" = numeric(), "no_sites" = numeric(), "count_diffs"= numeric(), "count_comparisons" = numeric(), "count_missing" = numeric(), "snp_id" = character())
pi_allpair_cmh <- tibble("pop" = character(), "pair" = numeric(), "chromosome" = character(), "window_pos_1" = numeric(), "window_pos_2" = numeric(), "avg_pi" = numeric(), "no_sites" = numeric(), "count_diffs"= numeric(), "count_comparisons" = numeric(), "count_missing" = numeric(), "snp_id" = character())
pi_allpair_xpehh <- tibble("pop" = character(), "pair" = numeric(), "chromosome" = character(), "window_pos_1" = numeric(), "window_pos_2" = numeric(), "avg_pi" = numeric(), "no_sites" = numeric(), "count_diffs"= numeric(), "count_comparisons" = numeric(), "count_missing" = numeric(), "snp_id" = character())

dif_table_10kbwin <- tibble("chromosome" = character(), "window_pos_1" = numeric(), "window_pos_2" = numeric(), "pair" = numeric(), "avg_pi_ag-nat" = numeric(), "avg_pi_ag" = numeric(), "avg_pi_nat" = numeric(), "snp_id" = character())
dif_table_cmh <- tibble("chromosome" = character(), "window_pos_1" = numeric(), "window_pos_2" = numeric(), "pair" = numeric(), "avg_pi_ag-nat" = numeric(), "avg_pi_ag" = numeric(), "avg_pi_nat" = numeric(), "snp_id" = character())
dif_table_xpehh <- tibble("chromosome" = character(), "window_pos_1" = numeric(), "window_pos_2" = numeric(), "pair" = numeric(), "avg_pi_ag-nat" = numeric(), "avg_pi_ag" = numeric(), "avg_pi_nat" = numeric(), "snp_id" = character())

for(p in 1:17){
  for(s in 1:16){
    file <- fread(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/perm_10kbwin/pixy_pi_pair_", p, "_scaf_", s, "_10kbwin.txt", sep = ""), sep="\t", header=T)
    file <- file %>%
      mutate("pair" = p, "snp_id" = paste(chromosome, ":", window_pos_1, sep = ""))
    
    pi_allpair_10kbwin <- add_row(pi_allpair_10kbwin, "pop" = file$pop, "pair" = file$pair, "chromosome" = file$chromosome, "window_pos_1" = file$window_pos_1, "window_pos_2" = file$window_pos_2, "avg_pi" = file$avg_pi, "no_sites" = file$no_sites, "count_diffs"= file$count_diffs, "count_comparisons" = file$count_comparisons, "count_missing" = file$count_missing, "snp_id" = file$snp_id)
    
    file_cmh <- fread(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/perm_cmh/pixy_pi_pair_", p, "_scaf_", s, "_cmh.txt", sep = ""), sep="\t", header=T)
    file_cmh <- file_cmh %>%
      mutate("pair" = p, "snp_id" = paste(chromosome, ":", window_pos_1, sep = ""))
    
    pi_allpair_cmh <- add_row(pi_allpair_cmh, "pop" = file_cmh$pop, "pair" = file_cmh$pair, "chromosome" = file_cmh$chromosome, "window_pos_1" = file_cmh$window_pos_1, "window_pos_2" = file_cmh$window_pos_2, "avg_pi" = file_cmh$avg_pi, "no_sites" = file_cmh$no_sites, "count_diffs"= file_cmh$count_diffs, "count_comparisons" = file_cmh$count_comparisons, "count_missing" = file_cmh$count_missing, "snp_id" = file_cmh$snp_id)
    
    file_xpehh <- fread(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/perm_xpehh/pixy_pi_pair_", p, "_scaf_", s, "_xpehh.txt", sep = ""), sep="\t", header=T)
    file_xpehh <- file_xpehh %>%
      mutate("pair" = p, "snp_id" = paste(chromosome, ":", window_pos_1, sep = ""))
    
    pi_allpair_xpehh <- add_row(pi_allpair_xpehh, "pop" = file_xpehh$pop, "pair" = file_xpehh$pair, "chromosome" = file_xpehh$chromosome, "window_pos_1" = file_xpehh$window_pos_1, "window_pos_2" = file_xpehh$window_pos_2, "avg_pi" = file_xpehh$avg_pi, "no_sites" = file_xpehh$no_sites, "count_diffs"= file_xpehh$count_diffs, "count_comparisons" = file_xpehh$count_comparisons, "count_missing" = file_xpehh$count_missing, "snp_id" = file_xpehh$snp_id)
    
    #get table of the pi dif for each pair for s errors
    dif_table_10kbwin_ag <- file %>%
      filter(pop == "AG") %>%
      mutate("avg_pi_ag" = avg_pi)
    
    dif_table_10kbwin_nat <- file %>%
      filter(pop == "NAT") %>%
      mutate("avg_pi_nat" = avg_pi)
    
    dif_table_10kbwin_file <- merge(dif_table_10kbwin_ag, dif_table_10kbwin_nat, by = "snp_id")
    
     dif_table_10kbwin <- add_row(dif_table_10kbwin, "chromosome" = dif_table_10kbwin_file$chromosome.x, "window_pos_1" = dif_table_10kbwin_file$window_pos_1.x, "window_pos_2" = dif_table_10kbwin_file$window_pos_2.x, "pair" = dif_table_10kbwin_file$pair.x, "avg_pi_ag-nat" = dif_table_10kbwin_file$avg_pi_ag - dif_table_10kbwin_file$avg_pi_nat, "avg_pi_ag" = dif_table_10kbwin_file$avg_pi_ag, "avg_pi_nat" = dif_table_10kbwin_file$avg_pi_nat, "snp_id" = dif_table_10kbwin_file$snp_id)
    
    dif_table_cmh_ag <- file_cmh %>%
      filter(pop == "AG") %>%
      mutate("avg_pi_ag" = avg_pi)
    
    dif_table_cmh_nat <- file_cmh %>%
      filter(pop == "NAT") %>%
      mutate("avg_pi_nat" = avg_pi)
    
    dif_table_cmh_file <- merge(dif_table_cmh_ag, dif_table_cmh_nat, by = "snp_id")
    
     dif_table_cmh <- add_row(dif_table_cmh, "chromosome" = dif_table_cmh_file$chromosome.x, "window_pos_1" = dif_table_cmh_file$window_pos_1.x, "window_pos_2" = dif_table_cmh_file$window_pos_2.x, "pair" = dif_table_cmh_file$pair.x, "avg_pi_ag-nat" = dif_table_cmh_file$avg_pi_ag - dif_table_cmh_file$avg_pi_nat, "avg_pi_ag" = dif_table_cmh_file$avg_pi_ag, "avg_pi_nat" = dif_table_cmh_file$avg_pi_nat, "snp_id" = dif_table_cmh_file$snp_id)
     
     dif_table_xpehh_ag <- file_xpehh %>%
      filter(pop == "AG") %>%
      mutate("avg_pi_ag" = avg_pi)
    
    dif_table_xpehh_nat <- file_xpehh %>%
      filter(pop == "NAT") %>%
      mutate("avg_pi_nat" = avg_pi)
    
    dif_table_xpehh_file <- merge(dif_table_xpehh_ag, dif_table_xpehh_nat, by = "snp_id")
    
     dif_table_xpehh <- add_row(dif_table_xpehh, "chromosome" = dif_table_xpehh_file$chromosome.x, "window_pos_1" = dif_table_xpehh_file$window_pos_1.x, "window_pos_2" = dif_table_xpehh_file$window_pos_2.x, "pair" = dif_table_xpehh_file$pair.x, "avg_pi_ag-nat" = dif_table_xpehh_file$avg_pi_ag - dif_table_xpehh_file$avg_pi_nat, "avg_pi_ag" = dif_table_xpehh_file$avg_pi_ag, "avg_pi_nat" = dif_table_xpehh_file$avg_pi_nat, "snp_id" = dif_table_xpehh_file$snp_id)
  }
}
#write.table(pi_allpair_10kbwin, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/pi_allpair_10kbwin.txt", quote = F, row.names = F)
#write.table(pi_allpair_cmh, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/pi_allpair_cmh.txt", quote = F, row.names = F)
write.table(dif_table_cmh, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/diftable_cmh.txt", quote = F, row.names = F)
write.table(dif_table_10kbwin, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/diftable_10kbwin.txt", quote = F, row.names = F)
write.table(dif_table_xpehh, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/diftable_xpehh.txt", quote = F, row.names = F)

#get average pi across pairs for each site
avgpi_10kbwin <- pi_allpair_10kbwin %>%
  group_by(chromosome, window_pos_1, window_pos_2, pop, snp_id) %>%
  summarise(pi_overpair = sum(count_diffs) / sum(count_comparisons))

avgpi_cmh <- pi_allpair_cmh %>%
  group_by(chromosome, window_pos_1, window_pos_2, pop, snp_id) %>%
  summarise(pi_overpair = sum(count_diffs) / sum(count_comparisons))

avgpi_xpehh <- pi_allpair_xpehh %>%
  group_by(chromosome, window_pos_1, window_pos_2, pop, snp_id) %>%
  summarise(pi_overpair = sum(count_diffs) / sum(count_comparisons))

write.table(avgpi_10kbwin, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/avgpi_10kbwins", quote = F, row.names = F)
write.table(avgpi_cmh, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/avgpi_cmh", quote = F, row.names = F)
write.table(avgpi_xpehh, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/avgpi_xpehh", quote = F, row.names = F)

#these are before avg over pairs
#pi_allpair_10kbwin <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/pi_allpair_10kbwin.txt")
#pi_allpair_cmh <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/pi_allpair_cmh.txt")
dif_table_cmh <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/diftable_cmh.txt")
dif_table_10kbwin <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/diftable_10kbwin.txt")
dif_table_xpehh <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/diftable_xpehh.txt")

#avg pi over pairs
avgpi_10kbwin <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/avgpi_10kbwins")
avgpi_cmh <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/avgpi_cmh")
avgpi_xpehh <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/avgpi_xpehh")

#get set of windows that fall into each AF bin
#NA when no snps in that bin
win_10kb_AF <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/dcgm_10kbsites_all_AF_final.txt")
colnames(win_10kb_AF) <- c("Chr", "Pos1", "Pos2", "AF")

#get list of windows to choose from for the permutation AF matching (and remove the ones not in avg pi tables)
#window is from .010 to .019
valid_wins <- unique(avgpi_10kbwin$snp_id)
binned_10kb_AF <- win_10kb_AF %>%
  filter(!is.na(AF)) %>%
  mutate(snp_id = paste(Chr, ":", Pos1, sep = ""), AF_bin = floor(AF * 100) / 100) %>%
  filter(snp_id %in% valid_wins) %>%
  group_by(AF_bin) %>%
  summarize(bins = list(snp_id), .groups = "drop")

win_cmh_AF <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/dcgm_sig_all_AF_final.txt")
colnames(win_cmh_AF) <- c("Chr", "Pos1", "Pos2", "AF")

#get af bins of the cmh 10kb win to match bins in binned_10kb_AF to 
win_cmh_AF <- win_cmh_AF %>%
  mutate(AF_bin = floor(AF * 100) / 100, snp_id = paste(Chr, ":", Pos1, sep = ""))

win_xpehh_AF <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/dcgm_sig_xpehh_all_AF_final.txt")
colnames(win_xpehh_AF) <- c("Chr", "Pos1", "Pos2", "AF")
win_xpehh_AF <- win_xpehh_AF %>%
  mutate(AF_bin = floor(AF * 100) / 100, snp_id = paste(Chr, ":", Pos1, sep = ""))

#Make lookup tables for perms - first is vector of values that map to the vector of names in the second vector so avg pis = snp ids
ag_pi_lookup  <- setNames(
  avgpi_10kbwin$pi_overpair[avgpi_10kbwin$pop == "AG"],
  avgpi_10kbwin$snp_id[avgpi_10kbwin$pop == "AG"]
)
nat_pi_lookup <- setNames(
  avgpi_10kbwin$pi_overpair[avgpi_10kbwin$pop == "NAT"],
  avgpi_10kbwin$snp_id[avgpi_10kbwin$pop == "NAT"]
)

#gets the row of that AF bin and then draws from the list of windows that fall into that AF bin 
perms <- function(AFtomatch, potential_wins){
  #set.seed(4599)
  #for allele frequencies for our cmh sites, choose one randomly drawn window from the list of windows in that AF bin
  wintouse <- AFtomatch %>%
    rowwise() %>%
    mutate("chosenwin_matchedAF" = sample(unlist(potential_wins$bins[which(potential_wins$AF_bin == AF_bin)]), size = 1)) 
  
  #look up wintouse in the pixy results now and get avg pi ag and nat from that perm win
  #snp id is the cmh snp ide and chosenwin_matchedAF is the perm win to use
  perm_pis <- wintouse %>%
    rowwise() %>%
    mutate("perm_avgpi_AG" = ag_pi_lookup[[chosenwin_matchedAF]], "perm_avgpi_NAT" = nat_pi_lookup[[chosenwin_matchedAF]])
  
  perm_pis$perm_avgpi_AG_exp <- perm_pis$perm_avgpi_AG^.5
  perm_pis$perm_avgpi_NAT_exp <- perm_pis$perm_avgpi_NAT^.5
  
  #calculate pi ag - pi nat (this is from avg of window, averaged over so not as accurate but this is going to be avgs and is good enough)
  perm_pis$'pi_ag-nat' <- perm_pis$perm_avgpi_AG - perm_pis$perm_avgpi_NAT
  perm_pis$'pi_ag_exp_-pi_nat_exp' <- perm_pis$perm_avgpi_AG_exp - perm_pis$perm_avgpi_NAT_exp  
  
  #calculate the slope for exp 
  obs_fit <- lm(perm_avgpi_NAT_exp ~ perm_avgpi_AG_exp, data = perm_pis)
  
  obs <- tibble("Slope" = as.numeric(coef(obs_fit)["perm_avgpi_AG_exp"]), "Intercept" = as.numeric(coef(obs_fit)["(Intercept)"]), "Avg_pi_ag_minus_nat" = mean(perm_pis$`pi_ag-nat`), "Avg_pi_ag_exp_minus_nat_exp" = mean(perm_pis$`pi_ag_exp_-pi_nat_exp`))
  
  return(obs)
}

results <- replicate(1000, perms(win_cmh_AF, binned_10kb_AF), simplify = T)
results <- t(results)
results <- unlist(results)
results_global <- tibble("Slope" = results[1:1000], "Intercept" = results[1001:2000], "Avg_pi_ag_minus_nat" = results[2001:3000], "Avg_pi_ag_exp_minus_nat_exp" = results[3001:4000])

results_x <- replicate(1000, perms(win_xpehh_AF, binned_10kb_AF), simplify = T)
results_x <- t(results_x)
results_x <- unlist(results_x)
results_global_x <- tibble("Slope" = results_x[1:1000], "Intercept" = results_x[1001:2000], "Avg_pi_ag_minus_nat" = results_x[2001:3000], "Avg_pi_ag_exp_minus_nat_exp" = results_x[3001:4000])

write.table(results_global, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/perm_cmh_sites.txt", quote = F, row.names = F)
write.table(results_global_x, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/perm_xpehh_sites.txt", quote = F, row.names = F)



#RUN FROM HERE NOW----------------------------
results_global <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/perm_cmh_sites.txt")
results_global_x <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/perm_xpehh_sites.txt")
#avg pi over pairs
avgpi_10kbwin <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/avgpi_10kbwins")
avgpi_cmh <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/avgpi_cmh")
avgpi_xpehh <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/avgpi_xpehh")

#get avg pi from the 10kb cmh windows
avgpi_cmh$snp_id <- paste(avgpi_cmh$chromosome, ":", avgpi_cmh$window_pos_1, sep = "")
obs_ag_pi_lookup  <- setNames(
    avgpi_cmh$pi_overpair[avgpi_cmh$pop == "AG"],
    avgpi_cmh$snp_id[avgpi_cmh$pop == "AG"])
obs_nat_pi_lookup <- setNames(
    avgpi_cmh$pi_overpair[avgpi_cmh$pop == "NAT"],
    avgpi_cmh$snp_id[avgpi_cmh$pop == "NAT"])
obs_pis_cmh <- avgpi_cmh %>%
  rowwise() %>%
  mutate("obs_avgpi_AG" = obs_ag_pi_lookup[[snp_id]], "obs_avgpi_NAT" = obs_nat_pi_lookup[[snp_id]]) %>%
  filter(pop == "NAT") %>% #filter to just one obs a site
  ungroup()

#and for xpehh sites
avgpi_xpehh$snp_id <- paste(avgpi_xpehh$chromosome, ":", avgpi_xpehh$window_pos_1, sep = "")
obs_ag_pi_lookup  <- setNames(
    avgpi_xpehh$pi_overpair[avgpi_xpehh$pop == "AG"],
    avgpi_xpehh$snp_id[avgpi_xpehh$pop == "AG"])
obs_nat_pi_lookup <- setNames(
    avgpi_xpehh$pi_overpair[avgpi_xpehh$pop == "NAT"],
    avgpi_xpehh$snp_id[avgpi_xpehh$pop == "NAT"])
obs_pis_xpehh <- avgpi_xpehh %>%
  rowwise() %>%
  mutate("obs_avgpi_AG" = obs_ag_pi_lookup[[snp_id]], "obs_avgpi_NAT" = obs_nat_pi_lookup[[snp_id]]) %>%
  filter(pop == "NAT") %>% #filter to just one obs a site
  ungroup()

#get exp dist
obs_pis_cmh$perm_avgpi_AG_exp <- obs_pis_cmh$obs_avgpi_AG^.5
obs_pis_cmh$perm_avgpi_NAT_exp <- obs_pis_cmh$obs_avgpi_NAT^.5

obs_pis_xpehh$perm_avgpi_AG_exp <- obs_pis_xpehh$obs_avgpi_AG^.5
obs_pis_xpehh$perm_avgpi_NAT_exp <- obs_pis_xpehh$obs_avgpi_NAT^.5

#get average ag-nat difference
obs_pis_cmh$perm_pi_AG_minus_NAT <- obs_pis_cmh$obs_avgpi_AG - obs_pis_cmh$obs_avgpi_NAT
obs_pis_xpehh$perm_pi_AG_minus_NAT <- obs_pis_xpehh$obs_avgpi_AG - obs_pis_xpehh$obs_avgpi_NAT

obs_pis_cmh$perm_pi_per_dif <- obs_pis_cmh$perm_pi_AG_minus_NAT/obs_pis_cmh$obs_avgpi_NAT
obs_pis_xpehh$perm_pi_per_dif <- obs_pis_xpehh$perm_pi_AG_minus_NAT/obs_pis_xpehh$obs_avgpi_NAT

#see how much difference is percentage of avg pi nat
pidifper <- ggplot() + geom_point(data = obs_pis_cmh, aes(x = perm_pi_AG_minus_NAT, y = perm_pi_per_dif, color = "CMH", alpha = .5)) + geom_point(data = obs_pis_xpehh, aes(x = perm_pi_AG_minus_NAT, y = perm_pi_per_dif, color = "XPEHH", alpha = .5)) + labs(title = "Avg pi ag - nat for real cmh/xpehh sites compared to its percentage of total pi", x = "Avg pi ag - nat", y = "pi_dif/nat_pi") + theme_bw() + scale_color_manual(values = c("#395D9CFF", "#8AD9B1FF")) + theme_bw() + theme(axis.text = element_text(size = 15))
ggsave(pidifper, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/pixy/ag-natpi_over_natpi.png", device = "png", dpi = 300, height = 6, width = 10, units = "in")

obs_fit <- lm(perm_avgpi_NAT_exp ~ perm_avgpi_AG_exp, data = obs_pis_cmh)
obs <- tibble("Slope" = as.numeric(coef(obs_fit)["perm_avgpi_AG_exp"]), "Intercept" = as.numeric(coef(obs_fit)["(Intercept)"]), "Avg_pi_ag_minus_nat" = mean(obs_pis_cmh$perm_pi_AG_minus_NAT))

obs_fit_x <- lm(perm_avgpi_NAT_exp ~ perm_avgpi_AG_exp, data = obs_pis_xpehh)
obs_x <- tibble("Slope" = as.numeric(coef(obs_fit_x)["perm_avgpi_AG_exp"]), "Intercept" = as.numeric(coef(obs_fit_x)["(Intercept)"]), "Avg_pi_ag_minus_nat" = mean(obs_pis_xpehh$perm_pi_AG_minus_NAT))

slopedist <- ggplot() + geom_histogram(data = results_global, aes(x = Slope)) + geom_vline(xintercept = obs$Slope) + labs(title = "Slopes of genomewide permuted 10kb windows (AF matched to cmh AFs) \ncompared to slope of observed. All exponentially distributed") + theme_bw() + theme(axis.text = element_text(size = 15))
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/pixy/genomewideperm_exptransform_slopeofagtonatpi.png", plot = slopedist, device = "png", dpi = 300, height = 6, width = 10, units = "in")

slopedist <- ggplot() + geom_histogram(data = results_global_x, aes(x = Slope)) + geom_vline(xintercept = obs_x$Slope) + labs(title = "Slopes of genomewide permuted 10kb windows (AF matched to xpehh AFs) \ncompared to slope of observed. All exponentially distributed") + theme_bw() + theme(axis.text = element_text(size = 15))
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/pixy/genomewideperm_exptransform_slopeofagtonatpi_xpehh.png", plot = slopedist, device = "png", dpi = 300, height = 6, width = 10, units = "in")

ggplot() + geom_point(data = obs_pis_cmh, aes(x = perm_avgpi_AG_exp, y = perm_avgpi_NAT_exp)) + geom_abline(xintercept = obs$Intercept, slope = obs$Slope) + geom_abline(xintercept = 0, slope = 1, linetype = "dashed") + labs(title = "Avg_pi^.5 and slope of 10kb win around cmh site")

ggplot() + geom_point(data = obs_pis_xpehh, aes(x = perm_avgpi_AG_exp, y = perm_avgpi_NAT_exp)) + geom_abline(xintercept = obs$Intercept, slope = obs$Slope) + geom_abline(xintercept = 0, slope = 1, linetype = "dashed") + labs(title = "Avg_pi^.5 and slope of 10kb win around xpehh site")

#figuring out how to plot slope and intercept and avg dif
slopesw_intercept <- ggplot() + geom_abline(slope = results_global$Slope, intercept = results_global$Intercept, alpha = .3) + geom_abline(slope = obs$Slope, intercept = obs$Intercept, color = "red") + xlim(0, .4) + ylim(0,.4) + labs(title = "Slopes of genomewide permuted 10kb windows (AF matched to cmh AFs) \ncompared to slope of observed. All exponentially distributed\nHave slope and intercept of inferred line")
ggsave(slopesw_intercept, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/pixy/genomewideperm_exptransform_slopeand_intofagtonatpi.png", device = "png", dpi = 300, height = 6, width = 10, units = "in")

ggplot() + geom_histogram(data = results_global, aes(x = Intercept)) + geom_vline(xintercept = obs$Intercept) + labs(title = "Intercepts of genomewide permuted 10kb windows (AF matched to cmh AFs) \ncompared to intercept of observed. All exponentially distributed")

pi_dif <- ggplot() + geom_histogram(data = results_global, aes(x = Avg_pi_ag_minus_nat)) + geom_vline(xintercept = obs$Avg_pi_ag_minus_nat, color = "blue") + labs(title = "ag pi ag - nat of genomewide permuted 10kb windows (AF matched to cmh AFs) \ncompared to intercept of observed. real pis (no transform)") + theme_bw() + theme(axis.text = element_text(size = 20))
ggsave(pi_dif, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/pixy/genomewideperm_ag-natpi.png", device = "png", dpi = 300, height = 4, width = 8, units = "in")

#get pvals of significantly different these values are
p_cmh_onesided <- (sum(results_global$Avg_pi_ag_minus_nat <= obs$Avg_pi_ag_minus_nat) + 1) /
                   (nrow(results_global) + 1)

#for non exponential distribution
obs_fit_notransform <- lm(obs_avgpi_NAT ~ obs_avgpi_AG, data = obs_pis_cmh)
obs_line_notransform <- tibble("Slope" = as.numeric(coef(obs_fit_notransform)["obs_avgpi_AG"]), "Intercept" = as.numeric(coef(obs_fit_notransform)["(Intercept)"]))
nonexp_trans_div <- ggplot() + geom_point(data = obs_pis_cmh, aes(x = obs_avgpi_AG, y = obs_avgpi_NAT)) + geom_abline(xintercept = obs$Intercept, slope = obs$Slope) + geom_abline(xintercept = 0, slope = 1, linetype = "dashed") + labs(title = "Slope of non exponentially transformed avg pi ag to avg pi nat")
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/pixy/genomewideperm_slopeofagtonatpi.png", plot = nonexp_trans_div, device = "png", dpi = 300, height = 6, width = 10, units = "in")

#looking at means difs - adapted from section around 1050
#calculate mean pi difference for ag and nat, standard error and CIs (from exp dist)
mean_difs_cmh <- dif_table_cmh %>%
    group_by(snp_id) %>%
    summarise("mean_pi_dif_frompairs" = mean(`avg_pi_ag-nat`, na.rm = T), "se_pi_dif" = sd(`avg_pi_ag-nat`, na.rm = T)/sqrt(n()), "mean_pi_dif_exp" = mean(`avg_pi_ag-nat`^.5, na.rm = T), "se_pi_dif_exp" = sd(`avg_pi_ag-nat`^.5, na.rm = T)/sqrt(n()), "se_pi_ag_exp" = sd(`avg_pi_ag`^.5, na.rm = T)/sqrt(n()), "mean_pi_ag_exp" = mean(`avg_pi_ag`^.5, na.rm = T), "se_pi_nat_exp" = sd(`avg_pi_nat`^.5, na.rm = T)/sqrt(n()), "mean_pi_nat_exp" = mean(`avg_pi_nat`^.5, na.rm = T), "CIupper_ag" = 1.96*se_pi_ag_exp + mean_pi_ag_exp, "CIlower_ag" =  mean_pi_ag_exp - 1.96*se_pi_ag_exp, "CIupper_nat" = 1.96*se_pi_nat_exp + mean_pi_nat_exp, "CIlower_nat" = mean_pi_nat_exp - 1.96*se_pi_nat_exp)

#get if the CIs cross 1:1 line
mean_difs_cmh <- mean_difs_cmh %>%
  mutate("y_overlap" = case_when(
    mean_pi_ag_exp >= CIlower_nat & mean_pi_ag_exp <= CIupper_nat ~ "Overlap",
    mean_pi_ag_exp < CIlower_nat & mean_pi_ag_exp < CIupper_nat ~ "Above",
    mean_pi_ag_exp > CIlower_nat & mean_pi_ag_exp > CIupper_nat ~ "Below"), "x_overlap" = case_when(
      mean_pi_nat_exp >= CIlower_ag & mean_pi_nat_exp <= CIupper_ag ~ "Overlap",
      mean_pi_nat_exp < CIlower_ag & mean_pi_nat_exp < CIupper_ag ~ "Below",
      mean_pi_nat_exp > CIlower_ag & mean_pi_nat_exp > CIupper_ag ~ "Above"
    )) %>% mutate("Selected_in" = case_when(
      x_overlap == "Overlap" & y_overlap == "Overlap" ~ "Both",
      x_overlap == "Above" & y_overlap == "Above" ~ "Nat",
      x_overlap == "Below" & y_overlap == "Below" ~ "Ag"
    ))

plot_leg <- mean_difs_cmh %>%
  group_by(Selected_in) %>%
  summarise("sample" = n()) %>%
  mutate("label" = paste(Selected_in, " (n = ", sample, ")", sep = ""))

mean_difs_cmh <- left_join(mean_difs_cmh, plot_leg, by = "Selected_in")

na_all <- ggplot(filter(mean_difs_cmh, is.na(mean_difs_cmh$Selected_in) == F), aes(x = mean_pi_ag_exp, y = mean_pi_nat_exp)) + geom_point(data = filter(mean_difs_cmh, mean_difs_cmh$Selected_in == "Both"), aes(color = label), alpha = .4) + geom_point(data = filter(mean_difs_cmh, mean_difs_cmh$Selected_in != "Both"), aes(color = label)) + geom_abline(intercept = 0, slope = 1, linetype = 2) +geom_abline(xintercept = obs$Intercept, slope = obs$Slope) + geom_errorbar(data = filter(mean_difs_cmh, mean_difs_cmh$Selected_in != "Both"), aes(ymin = CIlower_nat, ymax = CIupper_nat, color = label)) + geom_errorbar(data = filter(mean_difs_cmh, mean_difs_cmh$Selected_in != "Both"), aes(xmin = CIlower_ag, xmax = CIupper_ag, color = label)) + ggtitle("Mean Pi in Ag and Nat in 10kb win at cmh sites, exponential transformation") + scale_color_manual(values = c("#60CEACFF","#395D9CFF", "#382A54FF")) + xlim(0,.4) + ylim(0,.4) + theme_bw() + theme(axis.text = element_text(size = 15))
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/pixy/divplot_expmean_byPair_selected_naall.png", plot = na_all, device = "png", dpi = 300, height = 6, width = 10, units = "in") #replacing old plot

#repeat plots for xpehh----------------------------
slopesw_intercept <- ggplot() + geom_abline(slope = results_global_x$Slope, intercept = results_global_x$Intercept, alpha = .3) + geom_abline(slope = obs_x$Slope, intercept = obs_x$Intercept, color = "red") + xlim(0, .4) + ylim(0,.4) + labs(title = "Slopes of genomewide permuted 10kb windows (AF matched to xpehh AFs) \ncompared to slope of observed. All exponentially distributed\nHave slope and intercept of inferred line")
ggsave(slopesw_intercept, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/pixy/genomewideperm_exptransform_slopeand_intofagtonatpi_xpehh.png", device = "png", dpi = 300, height = 6, width = 10, units = "in")

ggplot() + geom_histogram(data = results_global_x, aes(x = Intercept)) + geom_vline(xintercept = obs_x$Intercept) + labs(title = "Intercepts of genomewide permuted 10kb windows (AF matched to xpehh AFs) \ncompared to intercept of observed. All exponentially distributed")

pi_dif <- ggplot() + geom_histogram(data = results_global_x, aes(x = Avg_pi_ag_minus_nat)) + geom_vline(xintercept = obs_x$Avg_pi_ag_minus_nat, color = "blue") + labs(title = "ag pi ag - nat of genomewide permuted 10kb windows (AF matched to xpehh AFs) \ncompared to intercept of observed. real pis (no transform)") + theme_bw() + theme(axis.text = element_text(size = 20))
ggsave(pi_dif, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/pixy/genomewideperm_ag-natpi_xpehh.png", device = "png", dpi = 300, height = 4, width = 8, units = "in")

#get pvals of significantly different these values are
p_xpehh_onesided <- (sum(results_global_x$Avg_pi_ag_minus_nat <= obs_x$Avg_pi_ag_minus_nat) + 1) /
                      (nrow(results_global_x) + 1)

#for non exponential distribution
obs_x_fit_notransform <- lm(obs_avgpi_NAT ~ obs_avgpi_AG, data = obs_pis_xpehh)
obs_x_line_notransform <- tibble("Slope" = as.numeric(coef(obs_x_fit_notransform)["obs_avgpi_AG"]), "Intercept" = as.numeric(coef(obs_x_fit_notransform)["(Intercept)"]))
nonexp_trans_div <- ggplot() + geom_point(data = obs_pis_xpehh, aes(x = obs_avgpi_AG, y = obs_avgpi_NAT)) + geom_abline(xintercept = obs_x$Intercept, slope = obs_x$Slope) + geom_abline(xintercept = 0, slope = 1, linetype = "dashed") + labs(title = "Slope of non exponentially transformed avg pi ag to avg pi nat")
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/pixy/genomewideperm_slopeofagtonatpi_xpehh.png", plot = nonexp_trans_div, device = "png", dpi = 300, height = 6, width = 10, units = "in")

#looking at means difs - adapted from section around 1050
#calculate mean pi difference for ag and nat, standard error and CIs (from exp dist)
mean_difs_xpehh <- dif_table_xpehh %>%
    group_by(snp_id) %>%
    summarise("mean_pi_dif_frompairs" = mean(`avg_pi_ag-nat`, na.rm = T), "se_pi_dif" = sd(`avg_pi_ag-nat`, na.rm = T)/sqrt(n()), "mean_pi_dif_exp" = mean(`avg_pi_ag-nat`^.5, na.rm = T), "se_pi_dif_exp" = sd(`avg_pi_ag-nat`^.5, na.rm = T)/sqrt(n()), "se_pi_ag_exp" = sd(`avg_pi_ag`^.5, na.rm = T)/sqrt(n()), "mean_pi_ag_exp" = mean(`avg_pi_ag`^.5, na.rm = T), "se_pi_nat_exp" = sd(`avg_pi_nat`^.5, na.rm = T)/sqrt(n()), "mean_pi_nat_exp" = mean(`avg_pi_nat`^.5, na.rm = T), "CIupper_ag" = 1.96*se_pi_ag_exp + mean_pi_ag_exp, "CIlower_ag" =  mean_pi_ag_exp - 1.96*se_pi_ag_exp, "CIupper_nat" = 1.96*se_pi_nat_exp + mean_pi_nat_exp, "CIlower_nat" = mean_pi_nat_exp - 1.96*se_pi_nat_exp)

#get if the CIs cross 1:1 line
mean_difs_xpehh <- mean_difs_xpehh %>%
  mutate("y_overlap" = case_when(
    mean_pi_ag_exp >= CIlower_nat & mean_pi_ag_exp <= CIupper_nat ~ "Overlap",
    mean_pi_ag_exp < CIlower_nat & mean_pi_ag_exp < CIupper_nat ~ "Above",
    mean_pi_ag_exp > CIlower_nat & mean_pi_ag_exp > CIupper_nat ~ "Below"), "x_overlap" = case_when(
      mean_pi_nat_exp >= CIlower_ag & mean_pi_nat_exp <= CIupper_ag ~ "Overlap",
      mean_pi_nat_exp < CIlower_ag & mean_pi_nat_exp < CIupper_ag ~ "Below",
      mean_pi_nat_exp > CIlower_ag & mean_pi_nat_exp > CIupper_ag ~ "Above"
    )) %>% mutate("Selected_in" = case_when(
      x_overlap == "Overlap" & y_overlap == "Overlap" ~ "Both",
      x_overlap == "Above" & y_overlap == "Above" ~ "Nat",
      x_overlap == "Below" & y_overlap == "Below" ~ "Ag"
    ))

plot_leg <- mean_difs_xpehh %>%
  group_by(Selected_in) %>%
  summarise("sample" = n()) %>%
  mutate("label" = paste(Selected_in, " (n = ", sample, ")", sep = ""))

mean_difs_xpehh <- left_join(mean_difs_xpehh, plot_leg, by = "Selected_in")

na_all <- ggplot(filter(mean_difs_xpehh, is.na(mean_difs_xpehh$Selected_in) == F), aes(x = mean_pi_ag_exp, y = mean_pi_nat_exp)) + geom_point(data = filter(mean_difs_xpehh, mean_difs_xpehh$Selected_in == "Both"), aes(color = label), alpha = .4) + geom_point(data = filter(mean_difs_xpehh, mean_difs_xpehh$Selected_in != "Both"), aes(color = label)) + geom_abline(intercept = 0, slope = 1, linetype = 2) +geom_abline(xintercept = obs_x$Intercept, slope = obs_x$Slope) + geom_errorbar(data = filter(mean_difs_xpehh, mean_difs_xpehh$Selected_in != "Both"), aes(ymin = CIlower_nat, ymax = CIupper_nat, color = label)) + geom_errorbar(data = filter(mean_difs_xpehh, mean_difs_xpehh$Selected_in != "Both"), aes(xmin = CIlower_ag, xmax = CIupper_ag, color = label)) + ggtitle("Mean Pi in Ag and Nat in 10kb win at xpehh sites, exponential transformation") + scale_color_manual(values = c("#60CEACFF","#395D9CFF", "#382A54FF")) + xlim(0,.4) + ylim(0,.4) + theme_bw() + theme(axis.text = element_text(size = 15))
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/pixy/divplot_expmean_byPair_selected_naall_xpehh.png", plot = na_all, device = "png", dpi = 300, height = 6, width = 10, units = "in") #replacing old plot
```

#Comparing within pair to pooled all pairs xpehh results #run on randi
because too much data to download

```{r}
#module load gcc/12.1.0
#module load R (version 4.2.1)
#R
library(data.table)
library(tidyverse)
xpehh <- tibble("pair" = numeric(), "num_sigxpehhvals" = numeric(), "num_sigxpehhvals_belowneg5" = numeric(), "num_sigxpehhvals_above5"= numeric(), "positions" = integer(), "norm_xpehh"= numeric())

for(p in 1:17){
  tmp <- fread(paste("/scratch/espolston/sel_scan_pair/Aggregated_pair_", p, "_xpehh.txt", sep =""))
  
  tmp_final <- tmp %>%
    filter(norm_xpehh >= 5 | norm_xpehh <= -5) %>%
    summarise("num_sigxpehhvals" = n(), "num_sigxpehhvals_belowneg5" = sum(norm_xpehh <= -5), "num_sigxpehhvals_above5" = sum(norm_xpehh >= 5), "positions" = list(locus_id), "norm_xpehh" = list(norm_xpehh), "pair" = p)
  
  xpehh <- rbind(xpehh, tmp_final)
}

#output <- select(xpehh, num_sigxpehhvals, num_sigxpehhvals_belowneg5, num_sigxpehhvals_above5, pair)
#write.table(output, file = "/scratch/espolston/sel_scan_pair/counts_sigxpehh_perpair.txt", quote = F, row.names = F)

saveRDS(xpehh, "/scratch/espolston/sel_scan_pair/counts_sigxpehh_perpair.txt")

sig_perpair <- readRDS("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/sel_scan/counts_sigxpehh_perpair.txt")

#looking at the shared hits between pairs
#make list of all possible hits
hitlist <- tibble("locus_id" = unique(unlist(sig_perpair$positions)), "totalpairs_sigforsite" = 0, "pairs_sigforsite" = vector("list", length(unique(unlist(sig_perpair$positions)))))

#go though each pairs list and add to the cumulative total of pairs and the list of which pairs it is in
for(p in 1:17){
  row <- which(sig_perpair$pair == p)
  pair_hits <- unlist(sig_perpair$positions[[row]])
  
  # find where these loci are in hitlist
  idx <- match(pair_hits, hitlist$locus_id)
  
  # increment the total count of pairs this locus is significant for
  hitlist$totalpairs_sigforsite[idx] <- hitlist$totalpairs_sigforsite[idx] + 1
  
  # append this pair number to each matching locus's list
  for(i in idx){
    hitlist$pairs_sigforsite[[i]] <- c(hitlist$pairs_sigforsite[[i]], p)
  }
}

saveRDS(hitlist, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/sitesshared_perpairxpehh.RDS")

numshared <- hitlist %>%
  filter(totalpairs_sigforsite > 1) #4694 sites are shared between 2 or more pairs

#----------- NOW RUNNING SAME THING LOCALLY FOR POOLED XPEHH----------------------
xpehh_pool <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/sel_scan/Aggregated_xpehh.txt")
xpehh_pool <- xpehh_pool %>%
    filter(norm_xpehh >= 5 | norm_xpehh <= -5) %>%
    summarise("num_sigxpehhvals" = n(), "num_sigxpehhvals_belowneg5" = sum(norm_xpehh <= -5), "num_sigxpehhvals_above5" = sum(norm_xpehh >= 5),  "pair" = "pooled") #"positions" = list(POS), "norm_xpehh" = list(norm_xpehh),

#compare pool and by pair distributions
bypair_distrib <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/sel_scan/counts_sigxpehh_perpair.txt")

dist_xpehh <- rbind(xpehh_pool, bypair_distrib)

#----------- NOW looking at dist of xpehh and H2H1 at cmh sites in POOLED XPEHH----------------------
xpehh_pool <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/sel_scan/Aggregated_xpehh.txt")
cmh_clumpsites <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clumphits_dcgm.txt")
cmh_clumpsites$locusid <- paste(cmh_clumpsites$CHROM, ":", cmh_clumpsites$BP, sep = "")

xpehh_pool <- filter(xpehh_pool, locus_id %in% cmh_clumpsites$locusid)
sum(xpehh_pool$norm_xpehh < -2)
sum(xpehh_pool$norm_xpehh > 2)

xpehh_cmh_dist <- ggplot(xpehh_pool) + geom_histogram(aes(x = norm_xpehh)) + labs(title = "Distribution of norm xpehh scores for cmh sites in pooled xpehh set")
ggsave(xpehh_cmh_dist, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/xpehhscores_forcmh.png", device = "png", dpi = 300, height = 6, width = 10, units = "in")

H12Scan_all <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/Aggregated_H12scan_all")
H12Scan_ag <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/Aggregated_H12scan_ag")
H12Scan_nat <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/Aggregated_H12scan_nat")

H12Scan_all <- mutate(H12Scan_all, "snp" = paste(chr, ":", win_center, sep =""))
H12all_cmh <- filter(H12Scan_all, snp %in% cmh_clumpsites$locusid)
H12Scan_ag <- mutate(H12Scan_ag, "snp" = paste(chr, ":", win_center, sep =""))
H12ag_cmh <- filter(H12Scan_ag, snp %in% cmh_clumpsites$locusid)
H12Scan_nat <- mutate(H12Scan_nat, "snp" = paste(chr, ":", win_center, sep =""))
H12nat_cmh <- filter(H12Scan_nat, snp %in% cmh_clumpsites$locusid)

#for the remaining sites, take top h12 value from the windows that contain the snp
for(site in 1:length(cmh_clumpsites$locusid)){
  if(cmh_clumpsites$locusid[site] %in% H12all_cmh$snp) {
  }else{
    row_all <- H12Scan_all[which(H12Scan_all$chr == cmh_clumpsites$CHROM[site] & cmh_clumpsites$BP[site] >= H12Scan_all$win_pos_1 & cmh_clumpsites$BP[site] <= H12Scan_all$win_pos_2),] 
    row_all <- arrange(row_all, desc(H12))
    H12all_cmh <- add_row(H12all_cmh, row_all[1,])
    
    row_ag <- H12Scan_ag[which(H12Scan_ag$chr == cmh_clumpsites$CHROM[site] & cmh_clumpsites$BP[site] >= H12Scan_ag$win_pos_1 & cmh_clumpsites$BP[site] <= H12Scan_ag$win_pos_2),] 
    row_ag <- arrange(row_ag, desc(H12))
    H12ag_cmh <- add_row(H12ag_cmh, row_ag[1,])
    
    row_nat <- H12Scan_nat[which(H12Scan_nat$chr == cmh_clumpsites$CHROM[site] & cmh_clumpsites$BP[site] >= H12Scan_nat$win_pos_1 & cmh_clumpsites$BP[site] <= H12Scan_nat$win_pos_2),] 
    row_nat <- arrange(row_nat, desc(H12))
    H12nat_cmh <- add_row(H12nat_cmh, row_nat[1,])
  }
}

sum(H12Scan_all$H12 > .1)
sum(H12Scan_ag$H12 > .1)
sum(H12Scan_nat$H12 > .1)

h12_cmh_dist_ag <- ggplot(H12Scan_ag) + geom_histogram(aes(x = H12)) + labs(title = "Distribution of H12 scores for cmh sites in ag only")
ggsave(h12_cmh_dist_ag, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/h12_ag_scores_forcmh.png", device = "png", dpi = 300, height = 6, width = 10, units = "in")
h12_cmh_dist_nat <- ggplot(H12Scan_nat) + geom_histogram(aes(x = H12)) + labs(title = "Distribution of H12 scores for cmh sites in nat only")
ggsave(h12_cmh_dist_nat, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/h12_nat_scores_forcmh.png", device = "png", dpi = 300, height = 6, width = 10, units = "in")
h12_cmh_dist_all <- ggplot(H12Scan_all) + geom_histogram(aes(x = H12)) + labs(title = "Distribution of H12 scores for cmh sites in pooled ag and nat run")
ggsave(h12_cmh_dist_all, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/h12_all_scores_forcmh.png", device = "png", dpi = 300, height = 6, width = 10, units = "in")

H2H1_cmh_dist_ag <- ggplot(H12Scan_ag) + geom_histogram(aes(x = H2H1)) + labs(title = "Distribution of H2H1 scores for cmh sites in ag only")
ggsave(H2H1_cmh_dist_ag, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/H2H1_ag_scores_forcmh.png", device = "png", dpi = 300, height = 6, width = 10, units = "in")
H2H1_cmh_dist_nat <- ggplot(H12Scan_nat) + geom_histogram(aes(x = H2H1)) + labs(title = "Distribution of H2H1 scores for cmh sites in nat only")
ggsave(H2H1_cmh_dist_nat, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/H2H1_nat_scores_forcmh.png", device = "png", dpi = 300, height = 6, width = 10, units = "in")
H2H1_cmh_dist_all <- ggplot(H12Scan_all) + geom_histogram(aes(x = H2H1)) + labs(title = "Distribution of H2H1 scores for cmh sites in pooled ag and nat run")
ggsave(H2H1_cmh_dist_all, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/H2H1_all_scores_forcmh.png", device = "png", dpi = 300, height = 6, width = 10, units = "in")

#----------- NOW looking at dist of cmh and h12 at xpehh clumped sites ----------------------
xpehh_clumpsites <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/xpehh_sigsites.txt")
dcgm <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/FDR_dcgm")

#switched to FDR instead of Bon bc only 2 bon p were < 1
row_lookup  <- setNames(dcgm_.05$FDR_p, dcgm_.05$SNP)
cmhbonp_forxpehhclumpsnps <- xpehh_clumpsites %>%
  rowwise() %>%
  mutate(row_lookup[[locus_id]])

sum(cmhbonp_forxpehhclumpsnps$`row_lookup[[locus_id]]` < .1)
plot <- ggplot(cmhbonp_forxpehhclumpsnps) + geom_histogram(aes(x = `row_lookup[[locus_id]]`)) + labs(title = "Distribution of FDR p for xpehh sites")
ggsave(plot, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/cmhscores_forxpehh.png", device = "png", dpi = 300, height = 6, width = 10, units = "in")

#h12 - some no matches bc in recomb filtered regions
xpehh_h12 <- filter(xpehh_clumpsites, locus_id %in% H12Scan_all$snp)
for(site in 1:length(xpehh_clumpsites$locus_id)){
  if(xpehh_clumpsites$locus_id[site] %in% xpehh_h12$locus_id) {
  }else{
    row_all <- H12Scan_all[which(H12Scan_all$chr == xpehh_clumpsites$CHROM[site] & xpehh_clumpsites$POS[site] >= H12Scan_all$win_pos_1 & xpehh_clumpsites$POS[site] <= H12Scan_all$win_pos_2),] 
    if(length(row_all$win_center) > 1){
      row_all <- arrange(row_all, desc(H12))
      H12all_cmh <- add_row(H12all_cmh, row_all[1,]) 
    }
  }
} #some reason this got no additional matches so I am giving up for now
```

#AFvapeR analysis

```{r}
#look at PR for every significant window to determine parallel, mulit and divergent cutoffs
eigen_res <- readRDS("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/afvapr/Afvaper_eigenresiduals.txt")
sigwin <- readRDS("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/afvapr/afvapr_sigwindows.txt")
alleigen_sig_windows <- as.vector(unlist(sigwin, use.names = FALSE))
alleigen_sig_windows <- as.data.table(alleigen_sig_windows)
colnames(alleigen_sig_windows) <- "window_id"

#add marker for significant on which
alleigen_sig_windows$sig_oneig <- 0
for(e in 1:17){
  alleigen_sig_windows$sig_oneig[which(alleigen_sig_windows$window_id %in% sigwin[[e]])] <- e
}
alleigen_sig_windows$sig_oneig <- factor(alleigen_sig_windows$sig_oneig)

alleigen_sig_windows <- add_column(alleigen_sig_windows, "Eigenvector" = rep(0, length(alleigen_sig_windows$window_id)),"Eigenvalue" = rep(0, length(alleigen_sig_windows$window_id)))

#AFVaper mean eigenvalue across sig windows in each eigenvector
for(e in 1:17){
  eigen_x <- filter(alleigen_sig_windows, sig_oneig == e)
  for(w in 1:length(eigen_x$window_id)){
    #get window we are looking at
    row <- which(alleigen_sig_windows$window_id == eigen_x$window_id[w])
    win <- alleigen_sig_windows$window_id[row] 
    
    #get eigenvalues across all eigenvectors for that window
    eigvals <- eigen_res[[win]]$eigenvals
    
    #update these values in the alleigen table
    alleigen_sig_windows$Eigenvector[row] <- 1
    alleigen_sig_windows$Eigenvalue[row] <- eigvals[1]
    alleigen_sig_windows <- add_row(alleigen_sig_windows, window_id = rep(alleigen_sig_windows$window_id[row], 16), sig_oneig = rep(alleigen_sig_windows$sig_oneig[row],16), Eigenvector = 2:17, Eigenvalue = eigvals[2:17])
  }
}

write.table(alleigen_sig_windows, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/afvapr/alleigen_sigwindows.txt", quote = F, row.names = F)

eig_means <- alleigen_sig_windows %>%
  group_by(sig_oneig, Eigenvector) %>%
  summarise("Eigenvalue_Mean" = mean(Eigenvalue)) %>%
  ungroup()

ggplot(eig_means[18:68,], aes(x = Eigenvector, y = Eigenvalue_Mean, color = sig_oneig, group = sig_oneig)) + geom_point() + geom_line()
ggplot(eig_means[69:289,], aes(x = Eigenvector, y = Eigenvalue_Mean, color = sig_oneig, group = sig_oneig)) + geom_point() + geom_line()
ggplot(eig_means, aes(x = Eigenvector, y = Eigenvalue_Mean, color = sig_oneig, group = sig_oneig)) + geom_point() + geom_line()

#parallel = eig 1
#multiparallel = eig 2-4
#divergent = eig 5-17


#------- Repeating for 95% and 99% cutoff-----
afvap_analysis <- function(cutoff){
  sigwin <- readRDS(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/afvapr/afvapr_sigwindows_", cutoff, ".txt", sep = ""))
  alleigen_sig_windows <- as.vector(unlist(sigwin, use.names = FALSE))
  alleigen_sig_windows <- as.data.table(alleigen_sig_windows)
  colnames(alleigen_sig_windows) <- "window_id"
  
  alleigen_sig_windows <- filter(alleigen_sig_windows, is.na(window_id) == F)

  #add marker for significant on which
  alleigen_sig_windows$sig_oneig <- 0
  for(e in 1:17){
    alleigen_sig_windows$sig_oneig[which(alleigen_sig_windows$window_id %in% sigwin[[e]])] <- e
  }
  alleigen_sig_windows$sig_oneig <- factor(alleigen_sig_windows$sig_oneig)

  alleigen_sig_windows <- add_column(alleigen_sig_windows, "Eigenvector" = rep(0, length(alleigen_sig_windows$window_id)),"Eigenvalue" = rep(0, length(alleigen_sig_windows$window_id)))

#AFVaper mean eigenvalue across sig windows in each eigenvector
  for(e in 1:17){
    eigen_x <- filter(alleigen_sig_windows, sig_oneig == e)
    if(length(eigen_x$window_id) > 0){
      for(w in 1:length(eigen_x$window_id)){
    #get window we are looking at
      row <- which(alleigen_sig_windows$window_id == eigen_x$window_id[w])
      win <- alleigen_sig_windows$window_id[row] 
    
    #get eigenvalues across all eigenvectors for that window
      eigvals <- eigen_res[[win]]$eigenvals
    
    #update these values in the alleigen table
      alleigen_sig_windows$Eigenvector[row] <- 1
      alleigen_sig_windows$Eigenvalue[row] <- eigvals[1]
      alleigen_sig_windows <- add_row(alleigen_sig_windows, window_id = rep(alleigen_sig_windows$window_id[row], 16), sig_oneig = rep(alleigen_sig_windows$sig_oneig[row],16), Eigenvector = 2:17, Eigenvalue = eigvals[2:17])
    }
    }
  }

  write.table(alleigen_sig_windows, paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/afvapr/alleigen_sigwindows_", cutoff, ".txt", sep = ""), quote = F, row.names = F)

  eig_means <- alleigen_sig_windows %>%
    group_by(sig_oneig, Eigenvector) %>%
    summarise("Eigenvalue_Mean" = mean(Eigenvalue)) %>%
    ungroup()
  
  write.table(eig_means, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/afvapr/mean_eig", cutoff, ".txt", sep = ""), quote = F, row.names = F)
}

ggplot(eig_means_95[1:68,], aes(x = Eigenvector, y = Eigenvalue_Mean, color = sig_oneig, group = sig_oneig)) + geom_point() + geom_line()
ggplot(eig_means_95[18:289,], aes(x = Eigenvector, y = Eigenvalue_Mean, color = sig_oneig, group = sig_oneig)) + geom_point() + geom_line()
#for 95% cutoff: 
#parallel: 1
#multiparallel: 2-6
#divergent: 7-17

ggplot(eig_means_99[1:68,], aes(x = Eigenvector, y = Eigenvalue_Mean, color = sig_oneig, group = sig_oneig)) + geom_point() + geom_line()
ggplot(eig_means_99[69:221,], aes(x = Eigenvector, y = Eigenvalue_Mean, color = sig_oneig, group = sig_oneig)) + geom_point() + geom_line()
#for 99% cutoff: 
#parallel: 1
#multiparallel: 2-4
#divergent: 8-17
```

#now looking into afvaper parallel candidates are some pairs found
together more often than others

```{r}
#--------------------------------ON EIGENVECTOR 1--------------------------------
afvaperres <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/afvapr/afvapr_eig1_parallel.txt", fill = T)

#look at number of times pair 1 is in a set, num times pair 1 and pair 2 are in a set
#do this for all combinations for up to 3 in a group
pairs <- paste("pair_", seq(1,17,by=1), sep = "")
pairs_of2 <- as.data.frame(t(combn(pairs, 2)))
pairs_of3 <- as.data.frame(t(combn(pairs, 3)))
pairs_of2$count <- 0
pairs_of3$count <- 0
pairs_of2$V1 <- factor(pairs_of2$V1, c("pair_1","pair_2",  "pair_3",  "pair_4",  "pair_5",  "pair_6",  "pair_7",  "pair_8",  "pair_9",  "pair_10", "pair_11", "pair_12", "pair_13","pair_14", "pair_15", "pair_16", "pair_17"))
pairs_of2$V2 <- factor(pairs_of2$V2, c("pair_1","pair_2",  "pair_3",  "pair_4",  "pair_5",  "pair_6",  "pair_7",  "pair_8",  "pair_9",  "pair_10", "pair_11", "pair_12", "pair_13","pair_14", "pair_15", "pair_16", "pair_17"))

#go over groups of 2 - check for in parallel or antiparallel
for(p in 1:length(pairs_of2$V1)){
  matches <- afvaperres %>%
    rowwise() %>%
    filter((pairs_of2$V1[p] %in% str_split(parallel_pops, ",")[[1]] & pairs_of2$V2[p] %in% str_split(parallel_pops, ",")[[1]]) | (pairs_of2$V1[p] %in% str_split(antiparallel_pops, ",")[[1]] & pairs_of2$V2[p] %in% str_split(antiparallel_pops, ",")[[1]])) %>%
    ungroup()
  
  pairs_of2$count[p] <- length(matches$window_id)
}

overlaps <- ggplot(pairs_of2, aes(x = V1, y = V2, fill = count)) + geom_tile(color = "white", linewidth = 0.5) + scale_fill_viridis_c(option = "magma") + coord_fixed() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = "Parallel eig 1 count of cooccuring pairs")
ggsave(overlaps, file= "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Afvaper_parallel_pairsfoundtogether.png", device = "png", height = 6, width = 12, units = "in")

#go over groups of 3
for(p in 1:length(pairs_of3$V1)){
  matches <- afvaperres %>%
    rowwise() %>%
    filter((pairs_of3$V1[p] %in% str_split(parallel_pops, ",")[[1]] & pairs_of3$V2[p] %in% str_split(parallel_pops, ",")[[1]] & pairs_of3$V3[p] %in% str_split(parallel_pops, ",")[[1]]) | (pairs_of3$V1[p] %in% str_split(antiparallel_pops, ",")[[1]] & pairs_of3$V2[p] %in% str_split(antiparallel_pops, ",")[[1]] & pairs_of3$V3[p] %in% str_split(antiparallel_pops, ",")[[1]])) %>%
    ungroup()
  
  pairs_of3$count[p] <- length(matches$window_id)
}

#go over groups of 4
pairs_of4 <- as.data.frame(t(combn(pairs, 4)))
pairs_of4$count <- 0
pairs_of4$V1 <- factor(pairs_of4$V1, c("pair_1","pair_2",  "pair_3",  "pair_4",  "pair_5",  "pair_6",  "pair_7",  "pair_8",  "pair_9",  "pair_10", "pair_11", "pair_12", "pair_13","pair_14", "pair_15", "pair_16", "pair_17"))
for(p in 1:length(pairs_of4$V1)){
  matches <- afvaperres %>%
    rowwise() %>%
    filter((pairs_of4$V1[p] %in% str_split(parallel_pops, ",")[[1]] & pairs_of4$V2[p] %in% str_split(parallel_pops, ",")[[1]] & pairs_of4$V3[p] %in% str_split(parallel_pops, ",")[[1]] & pairs_of4$V4[p] %in% str_split(parallel_pops, ",")[[1]]) | (pairs_of4$V1[p] %in% str_split(antiparallel_pops, ",")[[1]] & pairs_of4$V2[p] %in% str_split(antiparallel_pops, ",")[[1]] & pairs_of4$V3[p] %in% str_split(antiparallel_pops, ",")[[1]] & pairs_of4$V4[p] %in% str_split(antiparallel_pops, ",")[[1]])) %>%
    ungroup()
  
  pairs_of4$count[p] <- length(matches$window_id)
}
write.table(pairs_of4, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Afvaper_parallel_pairsfoundtogether_p4.txt", quote =F, row.names = F)


#--------------------------------ON EIGENVECTORS 1-4 (MULTIPARALLEL)-------------------------------------
afvaperres <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/afvapr/afvapr_eig1_parallel.txt", fill = T)
afvaperres <- filter(afvaperres, eigenvector == "Eig1") #filter down just to the vector we are looking at
afvaperres$eigenvalue_sum <- afvaperres$eigenvalue

for(e in 2:4){
  afvaperres_in <- fread(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/afvapr/afvapr_eig",e, "_parallel.txt", sep = ""), fill = T)
  afvaperres_in <- filter(afvaperres_in, eigenvector == paste("Eig", e, sep =""))
  
  afvaperres <- rbind(afvaperres, afvaperres_in)
}

pairs <- paste("pair_", seq(1,17,by=1), sep = "")
pairs_of2 <- as.data.frame(t(combn(pairs, 2)))
pairs_of2$count <- 0
pairs_of2$V1 <- factor(pairs_of2$V1, c("pair_1","pair_2",  "pair_3",  "pair_4",  "pair_5",  "pair_6",  "pair_7",  "pair_8",  "pair_9",  "pair_10", "pair_11", "pair_12", "pair_13","pair_14", "pair_15", "pair_16", "pair_17"))
pairs_of2$V2 <- factor(pairs_of2$V2, c("pair_1","pair_2",  "pair_3",  "pair_4",  "pair_5",  "pair_6",  "pair_7",  "pair_8",  "pair_9",  "pair_10", "pair_11", "pair_12", "pair_13","pair_14", "pair_15", "pair_16", "pair_17"))

#go over groups of 2 - check for in parallel 
#for(p in 1:length(pairs_of2$V1)){
#  matches <- afvaperres %>%
#    rowwise() %>%
#    filter((pairs_of2$V1[p] %in% str_split(parallel_pops, ",")[[1]] & pairs_of2$V2[p] %in% str_split(parallel_pops, ",")[[1]]) | (pairs_of2$V1[p] %in% str_split(antiparallel_pops, ",")[[1]] & pairs_of2$V2[p] %in% str_split(antiparallel_pops, ",")[[1]])) %>%
#    ungroup()
  
#  pairs_of2$count[p] <- length(matches$window_id)
#}

#overlaps <- ggplot(pairs_of2, aes(x = V1, y = V2, fill = count)) + geom_tile(color = "white", linewidth = 0.5) + scale_fill_viridis_c(option = "magma") + coord_fixed() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = "Multiparallel eig 2-4 count of cooccuring pairs")
#ggsave(overlaps, file= "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Afvaper_multiparallel_pairsfoundtogether.png", device = "png", height = 6, width = 12, units = "in")

for(p in 1:length(pairs_of2$V1)){
  matches <- afvaperres %>%
    rowwise() %>%
    filter((pairs_of2$V1[p] %in% str_split(parallel_pops, ",")[[1]] & pairs_of2$V2[p] %in% str_split(parallel_pops, ",")[[1]])) %>%
    ungroup()
  
  pairs_of2$count_onlyparallel[p] <- length(matches$window_id)
}

write.table(pairs_of2, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Afvaper_pairsfoundtogether.txt", quote = F, row.names = F)

p2gether<- ggplot(pairs_of2, aes(x = V1, y = V2, fill = count_onlyparallel)) + geom_tile(color = "white", linewidth = 0.5) + scale_fill_viridis_c(option = "mako") + coord_fixed() + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = "Multiparallel eig 1-4 count of cooccuring pairs\nonly parallel (not antiparallel)")
ggsave(p2gether, file= "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Afvaper_pairsfoundtogether.png", device = "png", height = 6, width = 12, units = "in")

#Now looking at count vs IBD
dist <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Pairwise_distances.txt")
dist <- dist %>%
  rowwise() %>%
  mutate("p1" = paste("pair_", str_split(pop1, "_")[[1]][2], sep = ""), "e1" = str_split(pop1, "_")[[1]][3], "p2" = paste("pair_", str_split(pop2, "_")[[1]][2], sep =""),"e2" = str_split(pop2, "_")[[1]][3], "pop_contrast" = case_when(
    e1 == e2 ~ "same_env",
    e1 != e2 ~ "ag-nat"
  )) %>%
  ungroup()
dist$p1 <- factor(dist$p1, c("pair_1","pair_2",  "pair_3",  "pair_4",  "pair_5",  "pair_6",  "pair_7",  "pair_8",  "pair_9",  "pair_10", "pair_11", "pair_12", "pair_13","pair_14", "pair_15", "pair_16", "pair_17"))
dist$p2 <- factor(dist$p2, c("pair_1","pair_2",  "pair_3",  "pair_4",  "pair_5",  "pair_6",  "pair_7",  "pair_8",  "pair_9",  "pair_10", "pair_11", "pair_12", "pair_13","pair_14", "pair_15", "pair_16", "pair_17"))

for(p in 1:length(pairs_of2$V1)){
  ind <- which((dist$p1 == pairs_of2$V1[p] & dist$p2 == pairs_of2$V2[p]) | (dist$p2 == pairs_of2$V1[p] & dist$p1 == pairs_of2$V2[p]))
  ind_use <- filter(dist[ind,])
  #have to average over ag and nat 
  pairs_of2$mean_pairdist[p] <- sum(ind_use$Geo_km)/length(ind_use$pop1)
}

library(ggpmisc)
paronly <- ggplot(pairs_of2, aes(x= mean_pairdist, y = count_onlyparallel)) + geom_point() + geom_smooth(method = lm) + labs(title = "count of pair cooccurance in parallel category vs population distance") + stat_poly_eq(formula = y ~ x, aes(label = after_stat(eq.label)), parse = TRUE) + theme_bw() + theme(axis.text = element_text(size = 15))
ggsave(paronly, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Afvaper_parallel_pairstogether_distance.png", device = "png", height = 6, width = 12, units = "in")

paronlypres <- ggplot(pairs_of2, aes(x= mean_pairdist, y = count_onlyparallel)) + geom_point() + geom_smooth(method = lm) + labs(title = "count of pair cooccurance in parallel category vs population distance") + theme_bw() + theme(axis.text = element_text(size = 20))
ggsave(paronlypres, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Afvaper_parallel_pairstogether_dist.png", device = "png", height = 6, width = 12, units = "in")

#ggplot(pairs_of2, aes(x= mean_pairdist, y = count)) + geom_point() + geom_smooth(method = lm) + labs(title = "count of pair cooccurance in parallel or antiparallel vs population distance") + stat_poly_eq(formula = y ~ x, aes(label = after_stat(eq.label)), parse = TRUE)

fit <- glm(pairs_of2$count_onlyparallel ~ pairs_of2$mean_pairdist)
summary(fit)
pseudo_r2 <- 1 - (fit$deviance / fit$null.deviance)

#------------Repeating ON EIGENVECTORS -4 (MULTIPARALLEL) for 95 and 99% cutoff --------------------
pairtogether <- function(cutoff){
  afvaperres <- fread(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/afvapr/afvapr_eig2_parallel", cutoff, ".txt", sep = ""), fill = T)
afvaperres <- filter(afvaperres, eigenvector == "Eig2") #filter down just to the vector we are looking at

for(e in 3:4){
  afvaperres_in <- fread(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/afvapr/afvapr_eig", e, "_parallel", cutoff, ".txt", sep = ""), fill = T)
  afvaperres_in <- filter(afvaperres_in, eigenvector == paste("Eig", e, sep =""))
  
  afvaperres <- rbind(afvaperres, afvaperres_in)
}

pairs <- paste("pair_", seq(1,17,by=1), sep = "")
pairs_of2 <- as.data.frame(t(combn(pairs, 2)))
pairs_of2$count <- 0
pairs_of2$V1 <- factor(pairs_of2$V1, c("pair_1","pair_2",  "pair_3",  "pair_4",  "pair_5",  "pair_6",  "pair_7",  "pair_8",  "pair_9",  "pair_10", "pair_11", "pair_12", "pair_13","pair_14", "pair_15", "pair_16", "pair_17"))
pairs_of2$V2 <- factor(pairs_of2$V2, c("pair_1","pair_2",  "pair_3",  "pair_4",  "pair_5",  "pair_6",  "pair_7",  "pair_8",  "pair_9",  "pair_10", "pair_11", "pair_12", "pair_13","pair_14", "pair_15", "pair_16", "pair_17"))

#go over groups of 2 - check for in parallel or antiparallel
for(p in 1:length(pairs_of2$V1)){
  matches <- afvaperres %>%
    rowwise() %>%
    filter((pairs_of2$V1[p] %in% str_split(parallel_pops, ",")[[1]] & pairs_of2$V2[p] %in% str_split(parallel_pops, ",")[[1]]) | (pairs_of2$V1[p] %in% str_split(antiparallel_pops, ",")[[1]] & pairs_of2$V2[p] %in% str_split(antiparallel_pops, ",")[[1]])) %>%
    ungroup()
  
  pairs_of2$count[p] <- length(matches$window_id)
}

overlaps <- ggplot(pairs_of2, aes(x = V1, y = V2, fill = count)) + geom_tile(color = "white", linewidth = 0.5) + scale_fill_viridis_c(option = "magma") + coord_fixed() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = paste("Multiparallel eig 2-4 count of cooccuring pairs at null ", cutoff))
ggsave(overlaps, file= paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Afvaper_multiparallel_pairsfoundtogether", cutoff, ".png"), device = "png", height = 6, width = 12, units = "in")

for(p in 1:length(pairs_of2$V1)){
  matches <- afvaperres %>%
    rowwise() %>%
    filter((pairs_of2$V1[p] %in% str_split(parallel_pops, ",")[[1]] & pairs_of2$V2[p] %in% str_split(parallel_pops, ",")[[1]])) %>%
    ungroup()
  
  pairs_of2$count_onlyparallel[p] <- length(matches$window_id)
}
ggplot(pairs_of2, aes(x = V1, y = V2, fill = count_onlyparallel)) + geom_tile(color = "white", linewidth = 0.5) + scale_fill_viridis_c(option = "magma") + coord_fixed() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = "Multiparallel eig 2-4 count of cooccuring pairs only parallel (not antiparallel)")

#Now looking at count vs IBD
dist <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Pairwise_distances.txt")
dist <- dist %>%
  rowwise() %>%
  mutate("p1" = paste("pair_", str_split(pop1, "_")[[1]][2], sep = ""), "e1" = str_split(pop1, "_")[[1]][3], "p2" = paste("pair_", str_split(pop2, "_")[[1]][2], sep =""),"e2" = str_split(pop2, "_")[[1]][3], "pop_contrast" = case_when(
    e1 == e2 ~ "same_env",
    e1 != e2 ~ "ag-nat"
  )) %>%
  ungroup()
dist$p1 <- factor(dist$p1, c("pair_1","pair_2",  "pair_3",  "pair_4",  "pair_5",  "pair_6",  "pair_7",  "pair_8",  "pair_9",  "pair_10", "pair_11", "pair_12", "pair_13","pair_14", "pair_15", "pair_16", "pair_17"))
dist$p2 <- factor(dist$p2, c("pair_1","pair_2",  "pair_3",  "pair_4",  "pair_5",  "pair_6",  "pair_7",  "pair_8",  "pair_9",  "pair_10", "pair_11", "pair_12", "pair_13","pair_14", "pair_15", "pair_16", "pair_17"))

for(p in 1:length(pairs_of2$V1)){
  ind <- which((dist$p1 == pairs_of2$V1[p] & dist$p2 == pairs_of2$V2[p]) | (dist$p2 == pairs_of2$V1[p] & dist$p1 == pairs_of2$V2[p]))
  ind_use <- filter(dist[ind,])
  #ind_use <- filter(dist[ind,], pop_contrast == "ag-nat")
  #have to average over ag and nat 
  pairs_of2$mean_pairdist[p] <- sum(ind_use$Geo_km)/length(ind_use$pop1)
}

paronly <- ggplot(pairs_of2, aes(x= mean_pairdist, y = count_onlyparallel)) + geom_point() + geom_smooth(method =lm) + labs(title = paste("count of pair cooccurance in parallel category vs population distance", cutoff, "%")) + stat_poly_eq(formula = y ~ x, aes(label = after_stat(eq.label)), parse = TRUE)
ggsave(paronly, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Afvaper_parallel_pairstogether_distance", cutoff, ".png", sep = ""), device = "png", height = 6, width = 12, units = "in")
fit <- glm(pairs_of2$count_onlyparallel ~ pairs_of2$mean_pairdist)
#ggplot(pairs_of2, aes(x= mean_pairdist, y = count)) + geom_point() + geom_smooth() + labs(title = "count of pair cooccurance in parallel or antiparallel vs population distance")  
}
```

#NOT THIS ONE- OLD
#looking at pairs together with environmental variables
```{r}
library(raster)
library(geodata)
library(corrplot)
library(fields)
library(tidyverse)

#https://www.stat.cmu.edu/~cshalizi/350/2008/lectures/14/lecture-14.pdf

#popcoords <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/Popcoordinates.tsv")
#popcoords$Pair_name <- paste("pair", popcoords$Pair, sep = "_")

#https://www.worldclim.org/data/bioclim.html
#deleted to save space
#temp <- worldclim_global(var = 'bio', res = .5, path = "/Users/libbypolston/Desktop/Climate/", lon = popcoords$Long, lat = popcoords$Lat) #download raster of climate data from worlclim

#climdataOG <- raster::extract(temp, popcoords[, c("Long", "Lat")], method = 'bilinear') #save data in in a separate frame
#colnames(climdataOG) <- c("ID", "B1_Ann_T","B2_Diurnal_Range","B3_Isothermality","B4_T_Seasonality","B5_MaxT_Wrmst_Month","B6_MinT_Cldst_Month","B7_T_Ann_Range","B8_T_Wettest_Qtr","B9_T_Driest_Qtr","B10_T_Wrmst_Qtr","B11_T_Cldst_Qtr","B12_Ann_Precip","B13_P_Wettest_Month","B14_P_Driest_Month","B15_P_Seasonality","B16_P_Wettest_Qtr","B17_P_Driest_Qtr","B18_P_Wrmst_Qtr","B19_P_Cldst_Qtr")
#clim <- cbind(popcoords, climdataOG)

#write.table(clim, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/bioclim_pairstogetherAFvaper.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#run from here ---------------------------------------
clim <- read.table(file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/bioclim_pairstogetherAFvaper.txt", header = T)

#get mean of paired sites variables
clim_acrosshab <- clim %>%
  group_by(Pair_name) %>%
  summarise("B1_Ann_T" = mean(B1_Ann_T),"B2_Diurnal_Range" = mean(B2_Diurnal_Range),"B3_Isothermality"= mean(B3_Isothermality),"B4_T_Seasonality"= mean(B4_T_Seasonality),"B5_MaxT_Wrmst_Month"= mean(B5_MaxT_Wrmst_Month),"B6_MinT_Cldst_Month"= mean(B6_MinT_Cldst_Month),"B7_T_Ann_Range"= mean(B7_T_Ann_Range),"B8_T_Wettest_Qtr"= mean(B8_T_Wettest_Qtr),"B9_T_Driest_Qtr"= mean(B9_T_Driest_Qtr),"B10_T_Wrmst_Qtr"= mean(B10_T_Wrmst_Qtr),"B11_T_Cldst_Qtr"= mean(B11_T_Cldst_Qtr),"B12_Ann_Precip"= mean(B12_Ann_Precip),"B13_P_Wettest_Month"= mean(B13_P_Wettest_Month),"B14_P_Driest_Month"= mean(B14_P_Driest_Month),"B15_P_Seasonality"= mean(B15_P_Seasonality),"B16_P_Wettest_Qtr"= mean(B16_P_Wettest_Qtr),"B17_P_Driest_Qtr"= mean(B17_P_Driest_Qtr),"B18_P_Wrmst_Qtr"= mean(B18_P_Wrmst_Qtr),"B19_P_Cldst_Qtr"= mean(B19_P_Cldst_Qtr))

#plot the correlation of the climate variables
png("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clim_afvaper/clim_corr.png")
print(corrplot(cor(clim_acrosshab[,c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)]), method = 'color'))
dev.off()
  
#plot clim var with pair count
pairs_of2 <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Afvaper_pairsfoundtogether.txt")

rsq_table <- tibble("variable" = "B0", "rsq" = 0, "pval" = NA_real_)

allvals_pairsof2 <- pairs_of2

for(c in 2:20){
  #dif between pair i and pair j
  varname <- paste0("dif_", colnames(clim_acrosshab)[c])  # e.g. "var_B1_Ann_T"
  varname1 <- paste0("i_", colnames(clim_acrosshab)[c])  # e.g. "i_B1_Ann_T"
  varname2 <- paste0("j_", colnames(clim_acrosshab)[c])  # e.g. "j_B1_Ann_T"
  pairs_of2 <- pairs_of2 %>%
    rowwise() %>%
    mutate(!!varname1 := clim_acrosshab[[which(clim_acrosshab$Pair_name == V1), c]], !!varname2 := clim_acrosshab[[which(clim_acrosshab$Pair_name == V2), c]], !!varname := !!sym(varname1)- !!sym(varname2)) %>%
    ungroup() %>%
    select(-!!varname1, -!!varname2)
  
  #keeping the values of each pair
  allvals_pairsof2 <- allvals_pairsof2 %>%
    rowwise() %>%
    mutate(!!varname1 := clim_acrosshab[[which(clim_acrosshab$Pair_name == V1), c]], !!varname2 := clim_acrosshab[[which(clim_acrosshab$Pair_name == V2), c]], !!varname := !!sym(varname1)- !!sym(varname2)) %>%
    ungroup()
  
  fit <- glm(pairs_of2$count_onlyparallel ~ pairs_of2[[varname]])
  pseudo_r2 <- 1 - (fit$deviance / fit$null.deviance)
  rsq_table <- add_row(rsq_table, "variable" = varname, "rsq" = pseudo_r2, "pval" = summary(fit)$coefficients[2, 4])
  
  climplot <- ggplot(pairs_of2, aes(x = .data[[varname]], y = count_onlyparallel)) + geom_point() + geom_smooth(method = lm) + labs(title = paste("count of pair cooccurance in parallel category vs difference between the 2 pairs", varname,"mean(ag+nat)")) + theme_bw() + theme(axis.text = element_text(size = 20))
ggsave(climplot, file = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clim_afvaper/Afvaper_parallel_pairstogether_clim_b", (c-1), ".png", sep = ""), device = "png", height = 6, width = 12, units = "in")
}

# remove placeholder row and apply multiple-test correction
rsq_table <- rsq_table %>%
  filter(variable != "B0") %>%
  mutate(
    padj_BH = p.adjust(pval, method = "BH"),          # Benjamini-Hochberg (FDR), recommended default
    padj_bonferroni = p.adjust(pval, method = "bonferroni"),  # stricter, controls FWER
    sig_BH = padj_BH < 0.05,
    sig_bonferroni = padj_bonferroni < 0.05
  )

#PCA of clim vars 
clim_acrosshab_pca <- clim_acrosshab %>%
  column_to_rownames(var = "Pair_name")

png(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clim_afvaper/climPCA_dif.png", sep = ""), width = 30, height = 30, units = "cm", res = 120)
print(biplot(prcomp(clim_acrosshab_pca, scale = T), cex = c(0.5,0.4)))
dev.off()

pca <- prcomp(clim_acrosshab_pca, scale = T)
summary(pca)

#now converting this into dif between pair i and j
pca_scores <- as.data.frame(pca$x)
pca_scores$pair_name <- rownames(pca_scores)

p2_contrasts <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Afvaper_pairsfoundtogether.txt") #reentering so not mixed w the single variable stuff
for(c in 1:17){
  #dif between pair i and pair j
  varname <- paste0("dif_PC", c)

  p2_contrasts <- p2_contrasts %>%
    rowwise() %>%
    mutate(val1 = pca_scores[[which(pca_scores$pair_name == V1), c]], val2 = pca_scores[[which(pca_scores$pair_name == V2), c]], !!varname := val1 - val2) %>%
    ungroup() %>%
    select(-val1, -val2)
}

#now adding in distance bw pairs and redoing pca
#have to get average distance between pair i ag,nat and pair j ag,nat
dist <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Pairwise_distances.txt")
dist <- dist %>%
  rowwise() %>%
  mutate("p1" = paste("pair_", str_split(pop1, "_")[[1]][2], sep = ""), "e1" = str_split(pop1, "_")[[1]][3], "p2" = paste("pair_", str_split(pop2, "_")[[1]][2], sep =""),"e2" = str_split(pop2, "_")[[1]][3], "pop_contrast" = case_when(
    e1 == e2 ~ "same_env",
    e1 != e2 ~ "ag-nat"
  )) %>%
  ungroup()
dist$p1 <- factor(dist$p1, c("pair_1","pair_2",  "pair_3",  "pair_4",  "pair_5",  "pair_6",  "pair_7",  "pair_8",  "pair_9",  "pair_10", "pair_11", "pair_12", "pair_13","pair_14", "pair_15", "pair_16", "pair_17"))
dist$p2 <- factor(dist$p2, c("pair_1","pair_2",  "pair_3",  "pair_4",  "pair_5",  "pair_6",  "pair_7",  "pair_8",  "pair_9",  "pair_10", "pair_11", "pair_12", "pair_13","pair_14", "pair_15", "pair_16", "pair_17"))
p2_contrasts$mean_pairdist <- 0

for(p in 1:length(p2_contrasts$V1)){
  ind <- which((dist$p1 == p2_contrasts$V1[p] & dist$p2 == p2_contrasts$V2[p]) | (dist$p2 == p2_contrasts$V1[p] & dist$p1 == p2_contrasts$V2[p]))
  ind_use <- filter(dist[ind,])
  #have to average over ag and nat 
  p2_contrasts$mean_pairdist[p] <- sum(ind_use$Geo_km)/length(ind_use$pop1)
}

p2_contrasts$Contrast <- paste(p2_contrasts$V1, p2_contrasts$V2, sep = ":")

#mean_pairdist not in the pca so adding here
dim_red_clim_dist <- glm(count_onlyparallel ~ dif_PC1 + dif_PC2 + mean_pairdist, data = p2_contrasts)
summary(dim_red_clim_dist)
dim_red_clim <- glm(count_onlyparallel ~ dif_PC1 + dif_PC2, data = p2_contrasts)
summary(dim_red_clim)
r2 <- 1 - (dim_red_clim_dist$deviance / dim_red_clim_dist$null.deviance) #38.5%
r2 <- 1 - (dim_red_clim$deviance / dim_red_clim$null.deviance) #38.5%


#now adding in the land and herbicide use-------------
#/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/amaranthus_predictors_for_libby.xlsx has another tab w info on the columns
other_env <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/amaranthus_predictors_for_libby.csv")
other_env$Pair_name <- paste("pair_", other_env$Pair, sep = "")

#get the difference in the proportion of ag land within 1000 km of the site (mean ag + nat at a site, site i - j)
p2_contrasts <- p2_contrasts %>%
  rowwise() %>%
  mutate("dif_proag2019" = mean(other_env$Ag_1km_2019[which(other_env$Pair_name == V1)]) - mean(other_env$Ag_1km_2019[which(other_env$Pair_name == V2)]), "dif_proag2016" = mean(other_env$Ag_1km_2016[which(other_env$Pair_name == V1)]) - mean(other_env$Ag_1km_2016[which(other_env$Pair_name == V2)])) %>%
  ungroup()

clim_dist_land <- glm(count_onlyparallel ~ dif_PC1 + dif_PC2 + mean_pairdist + dif_proag2019, data = p2_contrasts)
summary(clim_dist_land)
clim_land <- glm(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_proag2019, data = p2_contrasts)
summary(clim_land)

#herbicide use
#mean(ag + nat, pair i) - mean(ag + nat, pair j)
p2_contrasts <- p2_contrasts %>%
  rowwise() %>%
  mutate("dif_gly_kgha" = mean(other_env$glyphosate_kg_per_ha[which(other_env$Pair_name == V1)]) - mean(other_env$glyphosate_kg_per_ha[which(other_env$Pair_name == V2)]), "dif_atrazine_kgha" = mean(other_env$atrazine_kg_per_ha[which(other_env$Pair_name == V1)]) - mean(other_env$atrazine_kg_per_ha[which(other_env$Pair_name == V2)]), "dif_fomesafen_kgha" = mean(other_env$fomesafen_kg_per_ha[which(other_env$Pair_name == V1)]) - mean(other_env$fomesafen_kg_per_ha[which(other_env$Pair_name == V2)])) %>%
  ungroup()

clim_dist_land_allherb <- glm(count_onlyparallel ~ dif_PC1*dif_proag2019 + dif_PC2*dif_proag2019 + mean_pairdist + dif_proag2019 + dif_proag2019:dif_gly_kgha + dif_proag2019:dif_atrazine_kgha + dif_proag2019:dif_fomesafen_kgha, data = p2_contrasts)
summary(clim_dist_land_allherb)
clim_dist_allherb <- glm(count_onlyparallel ~ dif_PC1 + dif_PC2 + mean_pairdist + dif_gly_kgha + dif_atrazine_kgha + dif_fomesafen_kgha, data = p2_contrasts)
summary(clim_dist_allherb)
clim_allherb <- glm(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_gly_kgha + dif_atrazine_kgha + dif_fomesafen_kgha, data = p2_contrasts)
summary(clim_allherb)
#atrazine is the only sig herb? and mean pair dist and land use are not significant

#dim reducing all but distance (still not bc it is lat/long, maybe i should add in lat/long instead?) - mean(ag + nat)
clim_env_pca <- clim_acrosshab %>%
  rowwise() %>%
  mutate("pro_ag2019" = mean(other_env$Ag_1km_2019[which(other_env$Pair_name == Pair_name)]), "glyphosphate_kgha" = mean(other_env$glyphosate_kg_per_ha[which(other_env$Pair_name == Pair_name)]), "atrazine_kgha" = mean(other_env$atrazine_kg_per_ha[which(other_env$Pair_name == Pair_name)]), "fomesafen_kgha" = mean(other_env$fomesafen_kg_per_ha[which(other_env$Pair_name == Pair_name)])) %>%
  ungroup()

#dropping pairname
clim_env_pca <- clim_env_pca %>%
  column_to_rownames(var = "Pair_name")

pca_clim_land <- prcomp(clim_env_pca, scale = T)
summary(pca_clim_land)

#now converting this into dif between pair i and j
pca_scores <- as.data.frame(pca_clim_land$x)
pca_scores$pair_name <- rownames(pca_scores)

p2_alldimred <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Afvaper_pairsfoundtogether.txt") #reentering so not mixed
for(c in 1:17){
  #dif between pair i and pair j
  varname <- paste0("dif_PC", c)

  p2_alldimred <- p2_alldimred %>%
    rowwise() %>%
    mutate(val1 = pca_scores[[which(pca_scores$pair_name == V1), c]], val2 = pca_scores[[which(pca_scores$pair_name == V2), c]], !!varname := val1 - val2) %>%
    ungroup() %>%
    select(-val1, -val2)
}
dist <- dist %>%
  rowwise() %>%
  mutate("p1" = paste("pair_", str_split(pop1, "_")[[1]][2], sep = ""), "e1" = str_split(pop1, "_")[[1]][3], "p2" = paste("pair_", str_split(pop2, "_")[[1]][2], sep =""),"e2" = str_split(pop2, "_")[[1]][3], "pop_contrast" = case_when(
    e1 == e2 ~ "same_env",
    e1 != e2 ~ "ag-nat"
  )) %>%
  ungroup()
dist$p1 <- factor(dist$p1, c("pair_1","pair_2",  "pair_3",  "pair_4",  "pair_5",  "pair_6",  "pair_7",  "pair_8",  "pair_9",  "pair_10", "pair_11", "pair_12", "pair_13","pair_14", "pair_15", "pair_16", "pair_17"))
dist$p2 <- factor(dist$p2, c("pair_1","pair_2",  "pair_3",  "pair_4",  "pair_5",  "pair_6",  "pair_7",  "pair_8",  "pair_9",  "pair_10", "pair_11", "pair_12", "pair_13","pair_14", "pair_15", "pair_16", "pair_17"))
p2_alldimred$mean_pairdist <- 0

for(p in 1:length(p2_alldimred$V1)){
  ind <- which((dist$p1 == p2_alldimred$V1[p] & dist$p2 == p2_alldimred$V2[p]) | (dist$p2 == p2_alldimred$V1[p] & dist$p1 == p2_alldimred$V2[p]))
  ind_use <- filter(dist[ind,])
  #have to average over ag and nat 
  p2_alldimred$mean_pairdist[p] <- sum(ind_use$Geo_km)/length(ind_use$pop1)
}

p2_alldimred$Contrast <- paste(p2_alldimred$V1, p2_alldimred$V2, sep = ":")

#now looking at models
png(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clim_afvaper/clim_envPCA.png", sep = ""), width = 30, height = 30, units = "cm", res = 120)
print(biplot(prcomp(clim_env_pca, scale = T), cex = c(0.5,0.4)))
dev.off()

all <- glm(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC3 + dif_PC4 + mean_pairdist, data = p2_alldimred)
summary(all)

#PCA of ag and nat sites, so not mean across habitats
clim <- read.table(file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/bioclim_pairstogetherAFvaper.txt", header = T)
clim$ID <- paste("pair_", clim$Pair, "_", clim$Env, sep = "")
other_env <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/amaranthus_predictors_for_libby.csv")
other_env$ID <- paste("pair_", other_env$Pair, "_", other_env$Env, sep = "")

all_metrics <- merge(clim, other_env, by = "ID")
#dropping some cols (repeated or character)
all_metrics <- all_metrics %>%
  column_to_rownames(var = "ID")
all_metrics <- all_metrics[,-c(1,2,5,25,26,27,28,30,37,41,43,45,47,49,51,53,55,57,59)]

png(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clim_afvaper/PCA_byenv.png", sep = ""), width = 30, height = 30, units = "cm", res = 120)
print(biplot(prcomp(all_metrics, scale = T), cex = c(0.5,0.4)))
dev.off()
```

#NOT THIS ONE
#do clean run of the PCA AFvaper model choosing
```{r}
#get climate at the pair level (pair mean: ag and nat)
clim <- read.table(file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/bioclim_pairstogetherAFvaper.txt", header = T)
clim_acrosshab <- clim %>%
  group_by(Pair_name) %>%
  summarise("B1_Ann_T" = mean(B1_Ann_T),"B2_Diurnal_Range" = mean(B2_Diurnal_Range),"B3_Isothermality"= mean(B3_Isothermality),"B4_T_Seasonality"= mean(B4_T_Seasonality),"B5_MaxT_Wrmst_Month"= mean(B5_MaxT_Wrmst_Month),"B6_MinT_Cldst_Month"= mean(B6_MinT_Cldst_Month),"B7_T_Ann_Range"= mean(B7_T_Ann_Range),"B8_T_Wettest_Qtr"= mean(B8_T_Wettest_Qtr),"B9_T_Driest_Qtr"= mean(B9_T_Driest_Qtr),"B10_T_Wrmst_Qtr"= mean(B10_T_Wrmst_Qtr),"B11_T_Cldst_Qtr"= mean(B11_T_Cldst_Qtr),"B12_Ann_Precip"= mean(B12_Ann_Precip),"B13_P_Wettest_Month"= mean(B13_P_Wettest_Month),"B14_P_Driest_Month"= mean(B14_P_Driest_Month),"B15_P_Seasonality"= mean(B15_P_Seasonality),"B16_P_Wettest_Qtr"= mean(B16_P_Wettest_Qtr),"B17_P_Driest_Qtr"= mean(B17_P_Driest_Qtr),"B18_P_Wrmst_Qtr"= mean(B18_P_Wrmst_Qtr),"B19_P_Cldst_Qtr"= mean(B19_P_Cldst_Qtr))

#run PCA/dimension reduction on this climate
clim_acrosshab <- clim_acrosshab %>%
  column_to_rownames(var = "Pair_name")

png(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clim_afvaper/climPCA.png", sep = ""), width = 30, height = 30, units = "cm", res = 120)
print(biplot(prcomp(clim_acrosshab, scale = T), cex = c(0.5,0.4)))
dev.off()

pca <- prcomp(clim_acrosshab, scale = T)
summary(pca)

pca_scores <- as.data.frame(pca$x)
pca_scores$pair_name <- rownames(pca_scores)

#add in count data and get dif in pc scores between hab within same pair
p2_contrasts <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Afvaper_pairsfoundtogether.txt") 
for(c in 1:17){
  #dif between pair i and pair j
  varname <- paste0("dif_PC", c)

  p2_contrasts <- p2_contrasts %>%
    rowwise() %>%
    mutate(val1 = pca_scores[[which(pca_scores$pair_name == V1), c]], val2 = pca_scores[[which(pca_scores$pair_name == V2), c]], !!varname := abs(val1 - val2)) %>%
    ungroup() %>%
    select(-val1, -val2)
}

#get the difference in the proportion of ag land within 1000 km of the site (mean ag + nat at a site, site i - j)
#and difference in shannon div index
other_env <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/amaranthus_predictors_for_libby_fixed.csv") #had to edit out the commas from Kg of herbicide
other_env$Pair_name <- paste("pair_", other_env$Pair, sep = "")
p2_contrasts <- p2_contrasts %>%
  rowwise() %>%
  mutate("dif_proag2019" = abs(mean(other_env$Ag_1km_2019[which(other_env$Pair_name == V1)]) - mean(other_env$Ag_1km_2019[which(other_env$Pair_name == V2)])), "dif_proag2016" = abs(mean(other_env$Ag_1km_2016[which(other_env$Pair_name == V1)]) - mean(other_env$Ag_1km_2016[which(other_env$Pair_name == V2)])), "dif_Shannon2019" = abs(mean(other_env$Shannon_1km_2019[which(other_env$Pair_name == V1)]) - mean(other_env$Shannon_1km_2019[which(other_env$Pair_name == V2)]))) %>%
  ungroup()

#herbicide use
#ag site herbicide use
p2_contrasts <- p2_contrasts %>%
  rowwise() %>%
  mutate("ag_gly_kgha" = abs(other_env$glyphosate_kg_per_ha[which(other_env$Pair_name == V1 & other_env$Env == "AG")] - other_env$glyphosate_kg_per_ha[which(other_env$Pair_name == V2 & other_env$Env == "AG")]), "ag_atrazine_kgha" = abs(other_env$atrazine_kg_per_ha[which(other_env$Pair_name == V1 & other_env$Env == "AG")] - other_env$atrazine_kg_per_ha[which(other_env$Pair_name == V2 & other_env$Env == "AG")]), "ag_fomesafen_kgha" = abs(other_env$fomesafen_kg_per_ha[which(other_env$Pair_name == V1 & other_env$Env == "AG")] - other_env$fomesafen_kg_per_ha[which(other_env$Pair_name == V2 & other_env$Env == "AG")]), "ag_gly_kg" = abs(other_env$glyphosate_kg[which(other_env$Pair_name == V1 & other_env$Env == "AG")] - other_env$glyphosate_kg[which(other_env$Pair_name == V2 & other_env$Env == "AG")]), "ag_atrazine_kg" = abs(other_env$atrazine_kg[which(other_env$Pair_name == V1 & other_env$Env == "AG")] - other_env$atrazine_kg[which(other_env$Pair_name == V2 & other_env$Env == "AG")]), "ag_fomesafen_kg" = abs(other_env$fomesafen_kg[which(other_env$Pair_name == V1 & other_env$Env == "AG")] - other_env$fomesafen_kg[which(other_env$Pair_name == V2 & other_env$Env == "AG")])) %>%
  ungroup()

#pairwise distance 
dist <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Pairwise_distances.txt")
dist <- dist %>%
  rowwise() %>%
  mutate("p1" = paste("pair_", str_split(pop1, "_")[[1]][2], sep = ""), "e1" = str_split(pop1, "_")[[1]][3], "p2" = paste("pair_", str_split(pop2, "_")[[1]][2], sep =""),"e2" = str_split(pop2, "_")[[1]][3], "pop_contrast" = case_when(
    e1 == e2 ~ "same_env",
    e1 != e2 ~ "ag-nat"
  )) %>%
  ungroup()
dist$p1 <- factor(dist$p1, c("pair_1","pair_2",  "pair_3",  "pair_4",  "pair_5",  "pair_6",  "pair_7",  "pair_8",  "pair_9",  "pair_10", "pair_11", "pair_12", "pair_13","pair_14", "pair_15", "pair_16", "pair_17"))
dist$p2 <- factor(dist$p2, c("pair_1","pair_2",  "pair_3",  "pair_4",  "pair_5",  "pair_6",  "pair_7",  "pair_8",  "pair_9",  "pair_10", "pair_11", "pair_12", "pair_13","pair_14", "pair_15", "pair_16", "pair_17"))
p2_contrasts$mean_pairdist <- 0

for(p in 1:length(p2_contrasts$V1)){
  ind <- which((dist$p1 == p2_contrasts$V1[p] & dist$p2 == p2_contrasts$V2[p]) | (dist$p2 == p2_contrasts$V1[p] & dist$p1 == p2_contrasts$V2[p]))
  ind_use <- filter(dist[ind,])
  #have to average over ag and nat 
  p2_contrasts$mean_pairdist[p] <- sum(ind_use$Geo_km)/length(ind_use$pop1)
}

p2_contrasts$Contrast <- paste(p2_contrasts$V1, p2_contrasts$V2, sep = ":")

#add in relatedness/Fst
fst_wc <- tibble("pop1" = character(), "pop2" = character(), "chromosome" = character(), "window_pos_1" = numeric(), "window_pos_2" = numeric(), "avg_wc_fst" = numeric(), "no_snps" = numeric(),	"wc_fst_a" = numeric(),	"wc_fst_b" = numeric(),	"wc_fst_c" = numeric())

for(p in 1:16){
  pix_file <- fread(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/pixy_fst_50Mb_Scaffold_", p, ".txt", sep = ""), sep="\t", header=T)
  fst_wc <- fst_wc %>%
    add_row("pop1" = pix_file$pop1, "pop2" = pix_file$pop2, "chromosome" = pix_file$chromosome, "window_pos_1" = pix_file$window_pos_1, "window_pos_2" = pix_file$window_pos_2, "avg_wc_fst" = pix_file$avg_wc_fst, "no_snps" = pix_file$no_snps,	"wc_fst_a" = pix_file$wc_fst_a,	"wc_fst_b" = pix_file$wc_fst_b,	"wc_fst_c" = pix_file$wc_fst_c)
}

#average over all chromosomes and across p1-ag:p2-ag,p1-nat:p2-nat habitats within a pair
fst_wc_avg <- fst_wc %>%
  rowwise() %>%
  mutate("p1" = str_split(pop1, "_")[[1]][2], "p2" = str_split(pop2, "_")[[1]][2]) %>%
  ungroup() %>%
  mutate("contrast" = paste(p1, "_", p2, sep = "")) %>%
  group_by(contrast) %>%
  summarise("wcfst_allscaf" = sum(wc_fst_a) / (sum(wc_fst_a) + sum(wc_fst_b) + sum(wc_fst_c)), "pop1" = paste("pair_", first(p1), sep = ""), "pop2" = paste("pair_", first(p2), sep = "")) %>%
  ungroup()

for(p in 1:length(p2_contrasts$V1)){
  ind <- which((fst_wc_avg$pop1 == p2_contrasts$V1[p] & fst_wc_avg$pop2 == p2_contrasts$V2[p]) | (fst_wc_avg$pop2 == p2_contrasts$V1[p] & fst_wc_avg$pop1 == p2_contrasts$V2[p]))
  p2_contrasts$mean_fst[p] <- unlist(fst_wc_avg$wcfst_allscaf[ind])
}


univariate_fst <- ggplot(p2_contrasts, aes(x = mean_fst, y = count_onlyparallel)) + geom_point() + geom_smooth(method = lm) + labs(title = "count of pair cooccurance in parallel category vs population differentiation") + theme_bw() + theme(axis.text = element_text(size = 20))
ggsave(univariate_fst, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Afvaper_parallel_pairstogether_related.png", device = "png", height = 6, width = 12, units = "in")
fit <- glm(p2_contrasts$count_onlyparallel ~ p2_contrasts$mean_fst)
summary(fit)
pseudo_r2 <- 1 - (fit$deviance / fit$null.deviance)

#ok also doing with dim reduced herbicide data
#only using the herbicide data from the ag site
other_env_ag_pca <- other_env %>%
  filter(Env == "AG") 
other_env_ag_pca <- other_env_ag_pca %>%
  column_to_rownames(var = "Pair_name") 
other_env_ag_pca <- other_env_ag_pca[,seq(18,37, by = 1)]

png(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clim_afvaper/herbPCA.png", sep = ""), width = 30, height = 30, units = "cm", res = 120)
print(biplot(prcomp(other_env_ag_pca, scale = T), cex = c(0.5,0.4)))
dev.off()

herb_pca <- prcomp(other_env_ag_pca, scale = T)
summary(herb_pca)

herb_pca_scores <- as.data.frame(herb_pca$x)
herb_pca_scores$pair_name <- rownames(herb_pca_scores)

#now adding the dim red herb to the others
p2_contrasts_wPCherb <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Afvaper_pairsfoundtogether.txt") 
for(c in 1:17){
  #dif between pair i and pair j
  #herbicide
  varname <- paste0("herbdif_PC", c)

  p2_contrasts_wPCherb <- p2_contrasts_wPCherb %>%
    rowwise() %>%
    mutate(val1 = herb_pca_scores[[which(herb_pca_scores$pair_name == V1), c]], val2 = herb_pca_scores[[which(herb_pca_scores$pair_name == V2), c]], !!varname := abs(val1 - val2)) %>%
    ungroup() %>%
    select(-val1, -val2)
  
  #climate
  varname <- paste0("dif_PC", c)

  p2_contrasts_wPCherb <- p2_contrasts_wPCherb %>%
    rowwise() %>%
    mutate(val1 = pca_scores[[which(pca_scores$pair_name == V1), c]], val2 = pca_scores[[which(pca_scores$pair_name == V2), c]], !!varname := abs(val1 - val2)) %>%
    ungroup() %>%
    select(-val1, -val2)
}
#get the difference in the proportion of ag land within 1000 km of the site (mean ag + nat at a site, site i - j)
#and difference in shannon div index
p2_contrasts_wPCherb <- p2_contrasts_wPCherb %>%
  rowwise() %>%
  mutate("dif_proag2019" = abs(mean(other_env$Ag_1km_2019[which(other_env$Pair_name == V1)]) - mean(other_env$Ag_1km_2019[which(other_env$Pair_name == V2)])), "dif_proag2016" = abs(mean(other_env$Ag_1km_2016[which(other_env$Pair_name == V1)]) - mean(other_env$Ag_1km_2016[which(other_env$Pair_name == V2)])), "dif_Shannon2019" = abs(mean(other_env$Shannon_1km_2019[which(other_env$Pair_name == V1)]) - mean(other_env$Shannon_1km_2019[which(other_env$Pair_name == V2)]))) %>%
  ungroup()

p2_contrasts_wPCherb$mean_pairdist <- 0

for(p in 1:length(p2_contrasts_wPCherb$V1)){
  ind <- which((dist$p1 == p2_contrasts_wPCherb$V1[p] & dist$p2 == p2_contrasts_wPCherb$V2[p]) | (dist$p2 == p2_contrasts_wPCherb$V1[p] & dist$p1 == p2_contrasts_wPCherb$V2[p]))
  ind_use <- filter(dist[ind,])
  #have to average over ag and nat 
  p2_contrasts_wPCherb$mean_pairdist[p] <- sum(ind_use$Geo_km)/length(ind_use$pop1)
}

p2_contrasts_wPCherb$Contrast <- paste(p2_contrasts_wPCherb$V1, p2_contrasts_wPCherb$V2, sep = ":")

for(p in 1:length(p2_contrasts_wPCherb$V1)){
  ind <- which((fst_wc_avg$pop1 == p2_contrasts_wPCherb$V1[p] & fst_wc_avg$pop2 == p2_contrasts_wPCherb$V2[p]) | (fst_wc_avg$pop2 == p2_contrasts_wPCherb$V1[p] & fst_wc_avg$pop1 == p2_contrasts_wPCherb$V2[p]))
  p2_contrasts_wPCherb$mean_fst[p] <- unlist(fst_wc_avg$wcfst_allscaf[ind])
}

#------- Now have all the variables, start comparing models---------------------



#Not changing this section to mean fst--------------------------------------

#Now making models to compare between
#always keep pairdist
#look at interaction effect of land use and herbicides, and separately PC1:PC2...

everything <- glm(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + dif_proag2019 + dif_Shannon2019 + ag_gly_kgha + ag_atrazine_kgha + ag_fomesafen_kgha + ag_gly_kg + ag_atrazine_kg + ag_fomesafen_kg + dif_proag2019:ag_gly_kgha + dif_proag2019:ag_atrazine_kgha + dif_proag2019:ag_fomesafen_kgha + mean_pairdist, data = p2_contrasts, family = poisson)
summary(everything) 
#poisson AIC: 2381.5 Sig: all EXCEPT for  ag_atrazine_kg, ag_fomesafen_kg, dif_PC1:dif_PC2, dif_proag2019:ag_gly_kgha

herb <- glm(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + ag_gly_kgha + ag_atrazine_kgha + ag_fomesafen_kgha + ag_gly_kg + ag_atrazine_kg + ag_fomesafen_kg + dif_proag2019:ag_gly_kgha + dif_proag2019:ag_atrazine_kgha + dif_proag2019:ag_fomesafen_kgha + mean_pairdist, data = p2_contrasts, family = poisson) 
summary(herb)
#poisson AIC: 2529.4 Sig: PC1+2, ag_gly_kgha, ag_atrazine_kgha, ag_fomesafen_kg, ag_atrazine_kgha:dif_proag2019, mean_pairdist   

herb_kgha <- glm(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + ag_gly_kgha + ag_atrazine_kgha + ag_fomesafen_kgha + dif_proag2019:ag_gly_kgha + dif_proag2019:ag_atrazine_kgha + dif_proag2019:ag_fomesafen_kgha + mean_pairdist, data = p2_contrasts, family = poisson) 
summary(herb_kgha) 
#poisson AIC: 2538.2 Sig: PC1+2,ag_gly_kgha,ag_atrazine_kgha, ag_atrazine_kgha:dif_proag2019, mean_pairdist

land <- glm(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + dif_proag2019 + dif_Shannon2019 + mean_pairdist, data = p2_contrasts, family = poisson)
summary(land) 
#poisson AIC:2505.5 Sig: all of them

none <- glm(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + mean_pairdist, data = p2_contrasts, family = poisson)
summary(none) 
#poisson AIC: 2604.8 Sig: all except interaction of PC1+2

no_interaction <- glm(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_proag2019 + dif_Shannon2019 + mean_pairdist, data = p2_contrasts, family = poisson)
summary(no_interaction) 
#Poisson AIC: 2506.6 Sig: all of them

no_interaction2 <- glm(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_proag2019 + mean_pairdist, data = p2_contrasts, family = poisson)
summary(no_interaction2)
#poisson AIC: 2562.3 Sig: all of them

herbkgha_land <- glm(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + ag_gly_kgha + ag_atrazine_kgha + ag_fomesafen_kgha + mean_pairdist, data = p2_contrasts, family = poisson)
summary(herbkgha_land) 
#poisson AIC: 2549.3 Sig: all except for ag_fomesafen_kgha and interaction of PC1+2

#now looking with herb dim red
everything_herbred <- glm(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + dif_proag2019 + dif_Shannon2019 + herbdif_PC1 + herbdif_PC2 + herbdif_PC1:herbdif_PC2 + mean_pairdist, data = p2_contrasts_wPCherb, family = poisson)
summary(everything_herbred) 
#poisson AIC: 2464.8 Sig: all except interactions

herb_herbred_noint <- glm(count_onlyparallel ~ dif_PC1 + dif_PC2 + herbdif_PC1 + herbdif_PC2 + mean_pairdist, data = p2_contrasts_wPCherb, family = poisson)
summary(herb_herbred_noint)
#poisson AIC: 2595.4 Sig: all except for herb dif_PC2

noint_herbred <- glm(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_proag2019 + dif_Shannon2019 + herbdif_PC1 + mean_pairdist, data = p2_contrasts_wPCherb, family = poisson)
summary(noint_herbred) 
#poisson AIC: 2484.6 Sig: all

min_herbred <- glm(count_onlyparallel ~ dif_PC1 + dif_PC2 + herbdif_PC1 + mean_pairdist, data = p2_contrasts_wPCherb, family = poisson)
summary(min_herbred) 
#poisson AIC: 2593.4 sig: all

#checking if overdispersion
library(AER) 
dispersiontest(everything) #true dispersion is greater than 1 (disp = 10.79284)
dispersiontest(everything_herbred) #true dispersion is greater than 1 (disp = 11.55281)




#-------------redoing models with negative binomial------------------------------
#decided only dif_proag (not shannon) bc highly correlated
library(MASS)
everything <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + dif_proag2019 + ag_gly_kgha + ag_atrazine_kgha + ag_fomesafen_kgha + ag_gly_kg + ag_atrazine_kg + ag_fomesafen_kg + dif_proag2019:ag_gly_kgha + dif_proag2019:ag_atrazine_kgha + dif_proag2019:ag_fomesafen_kgha + mean_pairdist, data = p2_contrasts)
summary(everything) 
#AIC: 1428 Sig: PC2, pro ag, shannon, mean pair dist

everything <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + dif_Shannon2019 + ag_gly_kgha + ag_atrazine_kgha + ag_fomesafen_kgha + ag_gly_kg + ag_atrazine_kg + ag_fomesafen_kg + dif_proag2019:ag_gly_kgha + dif_proag2019:ag_atrazine_kgha + dif_proag2019:ag_fomesafen_kgha + mean_pairdist, data = p2_contrasts)
summary(everything) 
#AIC: 1425.9

herb <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + ag_gly_kgha + ag_atrazine_kgha + ag_fomesafen_kgha + ag_gly_kg + ag_atrazine_kg + ag_fomesafen_kg + dif_proag2019:ag_gly_kgha + dif_proag2019:ag_atrazine_kgha + dif_proag2019:ag_fomesafen_kgha + mean_pairdist, data = p2_contrasts) 
summary(herb)
#AIC: 1429.6 Sig: PC1+2, mean pair dist 

herb_kgha <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + ag_gly_kgha + ag_atrazine_kgha + ag_fomesafen_kgha + dif_proag2019:ag_gly_kgha + dif_proag2019:ag_atrazine_kgha + dif_proag2019:ag_fomesafen_kgha + mean_pairdist, data = p2_contrasts) 
summary(herb_kgha) 
#AIC: 1424.3 Sig: PC1+2, mean_pairdist

land <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + dif_proag2019 + mean_pairdist, data = p2_contrasts)
summary(land) 
#AIC: 1415.9 Sig: all except interaction of pc1+2

land <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + dif_Shannon2019 + mean_pairdist, data = p2_contrasts)
summary(land) 
#AIC: 1418.8

none <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + mean_pairdist, data = p2_contrasts)
summary(none) 
#AIC: 1417.6 Sig: all except interaction of PC1+2

#no_interaction <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_proag2019 + dif_Shannon2019 + mean_pairdist, data = p2_contrasts)
#summary(no_interaction) 
#AIC: 1409.6 Sig: all of them - not using bc highly corr

no_interaction2 <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_proag2019 + mean_pairdist, data = p2_contrasts)
summary(no_interaction2)
#AIC: 1414.1 Sig: all of them

no_interaction3 <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_Shannon2019 + mean_pairdist, data = p2_contrasts)
summary(no_interaction3) 
#AIC: 1416.9

herbkgha_land <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + ag_gly_kgha + ag_atrazine_kgha + ag_fomesafen_kgha + mean_pairdist, data = p2_contrasts)
summary(herbkgha_land) 
#AIC: 1420.4 Sig: PC1+2, ag_atrazine_kgha, mean_pairdist

everything_herbred <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + dif_proag2019 + herbdif_PC1 + herbdif_PC2 + herbdif_PC1:herbdif_PC2 + mean_pairdist, data = p2_contrasts_wPCherb)
summary(everything_herbred) 
#AIC: 1421.2 Sig: PC1+2, pro ag, shannon, mean pair dist

everything_herbred <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + dif_Shannon2019 + herbdif_PC1 + herbdif_PC2 + herbdif_PC1:herbdif_PC2 + mean_pairdist, data = p2_contrasts_wPCherb)
summary(everything_herbred) 
#AIC: 1423.8

herb_herbred_noint <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + herbdif_PC1 + herbdif_PC2 + mean_pairdist, data = p2_contrasts_wPCherb)
summary(herb_herbred_noint)
#AIC: 1419.3 Sig: all except for herb dif_PC2

noint_herbred <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_proag2019 + herbdif_PC1 + mean_pairdist, data = p2_contrasts_wPCherb)
summary(noint_herbred) 
#AIC: 1415.7 Sig: all except herb PC1

noint_herbred <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_Shannon2019 + herbdif_PC1 + mean_pairdist, data = p2_contrasts_wPCherb)
summary(noint_herbred) 
#AIC: 1417.9

min_herbred <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + herbdif_PC1 + mean_pairdist, data = p2_contrasts_wPCherb)
summary(min_herbred) 
#AIC: 1417.3 Sig: PC1+2, mean pair dist



#---------------neg bin with mean fst---------------------------
library(MASS)
everything <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + dif_proag2019 + ag_gly_kgha + ag_atrazine_kgha + ag_fomesafen_kgha + ag_gly_kg + ag_atrazine_kg + ag_fomesafen_kg + dif_proag2019:ag_gly_kgha + dif_proag2019:ag_atrazine_kgha + dif_proag2019:ag_fomesafen_kgha + mean_fst, data = p2_contrasts)
summary(everything) 
#AIC: 1428.5 Sig: PC1,PC2, mean fst

everything <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + dif_Shannon2019 + ag_gly_kgha + ag_atrazine_kgha + ag_fomesafen_kgha + ag_gly_kg + ag_atrazine_kg + ag_fomesafen_kg + dif_proag2019:ag_gly_kgha + dif_proag2019:ag_atrazine_kgha + dif_proag2019:ag_fomesafen_kgha + mean_fst, data = p2_contrasts)
summary(everything) 
#AIC: 1425.7

herb <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + ag_gly_kgha + ag_atrazine_kgha + ag_fomesafen_kgha + ag_gly_kg + ag_atrazine_kg + ag_fomesafen_kg + dif_proag2019:ag_gly_kgha + dif_proag2019:ag_atrazine_kgha + dif_proag2019:ag_fomesafen_kgha + mean_fst, data = p2_contrasts) 
summary(herb)
#AIC: 1428.8 Sig: PC1+2, mean fst 

herb_kgha <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + ag_gly_kgha + ag_atrazine_kgha + ag_fomesafen_kgha + dif_proag2019:ag_gly_kgha + dif_proag2019:ag_atrazine_kgha + dif_proag2019:ag_fomesafen_kgha + mean_fst, data = p2_contrasts) 
summary(herb_kgha) 
#AIC: 1424 Sig: PC1+2, mean_fst

land <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + dif_proag2019 + mean_fst, data = p2_contrasts)
summary(land) 
#AIC: 1418 Sig: all except interaction of pc1+2 AND mean_fst

land <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + dif_Shannon2019 + mean_fst, data = p2_contrasts)
summary(land) 
#AIC: 1421.4 Sig: PC2, mean_fst

none <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + mean_fst, data = p2_contrasts)
summary(none) 
#AIC: 1419.5 Sig: all except interaction of PC1+2

no_interaction2 <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_proag2019 + mean_fst, data = p2_contrasts)
summary(no_interaction2)
#AIC: 1416.1 Sig: all except mean_fst

no_interaction3 <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_Shannon2019 + mean_fst, data = p2_contrasts)
summary(no_interaction3) 
#AIC: 1419.4 Sig: all but dif Shannon div

herbkgha_land <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + ag_gly_kgha + ag_atrazine_kgha + ag_fomesafen_kgha + mean_fst, data = p2_contrasts)
summary(herbkgha_land) 
#AIC: 1420.6 Sig: PC2, ag_gly_kgha, ag_atrazine_kgha, mean_fst

everything_herbred <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + dif_proag2019 + herbdif_PC1 + herbdif_PC2 + herbdif_PC1:herbdif_PC2 + mean_fst, data = p2_contrasts_wPCherb)
summary(everything_herbred) 
#AIC: 1422.2 Sig: PC1+2, pro ag, mean_fst

everything_herbred <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + dif_Shannon2019 + herbdif_PC1 + herbdif_PC2 + herbdif_PC1:herbdif_PC2 + mean_fst, data = p2_contrasts_wPCherb)
summary(everything_herbred) 
#AIC: 1425.6 Sig: PC1+2, mean_fst

herb_herbred_noint <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + herbdif_PC1 + herbdif_PC2 + mean_fst, data = p2_contrasts_wPCherb)
summary(herb_herbred_noint)
#AIC: 1420.4 Sig: PC1+2, mean_fst

noint_herbred <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_proag2019 + herbdif_PC1 + mean_fst, data = p2_contrasts_wPCherb)
summary(noint_herbred) 
#AIC: 1417.5 Sig: all except herb PC1

noint_herbred <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_Shannon2019 + herbdif_PC1 + mean_fst, data = p2_contrasts_wPCherb)
summary(noint_herbred) 
#AIC: 1420.3 Sig: all except herb PC1

min_herbred <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + herbdif_PC1 + mean_fst, data = p2_contrasts_wPCherb)
summary(min_herbred) 
#AIC: 1418.7 Sig: PC1+2, mean fst

png("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clim_afvaper/env_corr.png")
corrplot(cor(other_env[,-c(1,2,3,14,38)]), method = 'color')
dev.off()

#decided on no interaction 2 model
no_interaction2 <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_proag2019 + mean_pairdist, data = p2_contrasts)
summary(no_interaction2)

pred_data <- data.frame(
  dif_PC2 = seq(min(p2_contrasts$dif_PC2, na.rm = TRUE),
                max(p2_contrasts$dif_PC2, na.rm = TRUE),
                length.out = 100),
  dif_PC1 = mean(p2_contrasts$dif_PC1, na.rm = TRUE),
  dif_proag2019 = mean(p2_contrasts$dif_proag2019, na.rm = TRUE),
  mean_pairdist = mean(p2_contrasts$mean_pairdist, na.rm = TRUE)
)

# Get predictions on the response scale, with standard errors
pred <- predict(no_interaction2, newdata = pred_data, type = "link", se.fit = TRUE)
pred_data$fit <- exp(pred$fit)
pred_data$lower <- exp(pred$fit - 1.96 * pred$se.fit)
pred_data$upper <- exp(pred$fit + 1.96 * pred$se.fit)

climplot_PC2 <- ggplot() + geom_point(data = p2_contrasts, aes(x = dif_PC2, y = count_onlyparallel)) + geom_ribbon(data = pred_data, aes(x = dif_PC2, ymin = lower, ymax = upper), alpha = 0.2) + geom_line(data = pred_data, aes(x = dif_PC2, y = fit), color = "blue", linewidth = 1) + labs(title = "count of pair cooccurance in parallel category vs dif_PC2 (other variables are held at their mean value") + theme_bw() + theme(axis.text = element_text(size = 20))
ggsave(climplot, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clim_afvaper/Afvaper_parallel_pairstogether_env_final_PC2.png", device = "png", height = 6, width = 12, units = "in")

#do for all sig variables
pred2_data <- data.frame(
  dif_PC1 = seq(min(p2_contrasts$dif_PC1, na.rm = TRUE),
                max(p2_contrasts$dif_PC1, na.rm = TRUE),
                length.out = 100),
  dif_PC2 = mean(p2_contrasts$dif_PC2, na.rm = TRUE),
  dif_proag2019 = mean(p2_contrasts$dif_proag2019, na.rm = TRUE),
  mean_pairdist = mean(p2_contrasts$mean_pairdist, na.rm = TRUE)
)

# Get predictions on the response scale, with standard errors
pred2 <- predict(no_interaction2, newdata = pred2_data, type = "link", se.fit = TRUE)
pred2_data$fit <- exp(pred2$fit)
pred2_data$lower <- exp(pred2$fit - 1.96 * pred2$se.fit)
pred2_data$upper <- exp(pred2$fit + 1.96 * pred2$se.fit)

climplot_PC1 <- ggplot() + geom_point(data = p2_contrasts, aes(x = dif_PC1, y = count_onlyparallel)) + geom_ribbon(data = pred2_data, aes(x = dif_PC1, ymin = lower, ymax = upper), alpha = 0.2) + geom_line(data = pred2_data, aes(x = dif_PC1, y = fit), color = "blue", linewidth = 1) + labs(title = "Difference in PC1 for pop i,j", x = "PC1,i - PC1,j", y = "# parallel windows") + theme_bw() + theme(axis.text = element_text(size = 20))

pred3_data <- data.frame(
  dif_proag2019 = seq(min(p2_contrasts$dif_proag2019, na.rm = TRUE),
                max(p2_contrasts$dif_proag2019, na.rm = TRUE),
                length.out = 100),
  dif_PC2 = mean(p2_contrasts$dif_PC2, na.rm = TRUE),
  dif_PC1 = mean(p2_contrasts$dif_PC1, na.rm = TRUE),
  mean_pairdist = mean(p2_contrasts$mean_pairdist, na.rm = TRUE)
)

# Get predictions on the response scale, with standard errors
pred3 <- predict(no_interaction2, newdata = pred3_data, type = "link", se.fit = TRUE)
pred3_data$fit <- exp(pred3$fit)
pred3_data$lower <- exp(pred3$fit - 1.96 * pred3$se.fit)
pred3_data$upper <- exp(pred3$fit + 1.96 * pred3$se.fit)

climplot_ag <- ggplot() + geom_point(data = p2_contrasts, aes(x = dif_proag2019, y = count_onlyparallel)) + geom_ribbon(data = pred3_data, aes(x = dif_proag2019, ymin = lower, ymax = upper), alpha = 0.2) + geom_line(data = pred3_data, aes(x = dif_proag2019, y = fit), color = "blue", linewidth = 1) + labs(title = "Difference in % ag in 1km of pop i,j", x = "% ag in 1km", y = "# parallel windows") + theme_bw() + theme(axis.text = element_text(size = 20))

pred4_data <- data.frame(
  mean_pairdist = seq(min(p2_contrasts$mean_pairdist, na.rm = TRUE),
                max(p2_contrasts$mean_pairdist, na.rm = TRUE),
                length.out = 100),
  dif_PC2 = mean(p2_contrasts$dif_PC2, na.rm = TRUE),
  dif_PC1 = mean(p2_contrasts$dif_PC1, na.rm = TRUE),
  dif_proag2019 = mean(p2_contrasts$dif_proag2019, na.rm = TRUE)
)

# Get predictions on the response scale, with standard errors
pred4 <- predict(no_interaction2, newdata = pred4_data, type = "link", se.fit = TRUE)
pred4_data$fit <- exp(pred4$fit)
pred4_data$lower <- exp(pred4$fit - 1.96 * pred4$se.fit)
pred4_data$upper <- exp(pred4$fit + 1.96 * pred4$se.fit)

climplot_dist <- ggplot() + geom_point(data = p2_contrasts, aes(x = mean_pairdist, y = count_onlyparallel)) + geom_ribbon(data = pred4_data, aes(x = dif_proag2019, ymin = lower, ymax = upper), alpha = 0.2) + geom_line(data = pred4_data, aes(x = mean_pairdist, y = fit), color = "blue", linewidth = 1) + labs(title = "Distance between pop i,j", x = "Distance (km)", y = "# parallel windows") + theme_bw() + theme(axis.text = element_text(size = 20))

#facet plot
library(patchwork)
all <- climplot_dist / climplot_ag / climplot_PC1
ggsave(all, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clim_afvaper/Afvaper_parallel_pairstogether_env_final_others.png", device = "png", height = 12, width = 8, units = "in")

#get rsq of this model
pseudo_r2 <- 1 - (no_interaction2$deviance / no_interaction2$null.deviance) #.404
pval <- summary(no_interaction2)$coefficients[2, 4]





#----------- Now getting 26 year avg of herbicide use --------------
#find the fips for herbicide use and acres of cropland
#pair_coord <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/paircoord.txt")
#pair_coord <- pair_coord %>%
#  mutate(V5 = paste(V3, V4, sep = "_")) %>%
#  select(V5, V1, V2)

#fips is constant across our years so we can use 2010
# latlong2fips <- function(latitude, longitude) {
#   url <- "https://geo.fcc.gov/api/census/area?lat=%s&lon=%s&censusYear=2010&format=json"
#   url <- sprintf(url, latitude, longitude)
#   json <- RJSONIO::fromJSON(url)
#   as.character(json[["results"]][[1]][["county_fips"]])
# }
# 
# n <- nrow(pair_coord)
# pair_coord$fips <- rep(NA_character_, n)
# 
# for (i in 1:n){
#   pair_coord$fips[i] <- latlong2fips(pair_coord$V1[i], pair_coord$V2[i])
# }
# 
# pair_coord <- pair_coord %>%
#   rowwise() %>%
#   mutate("fips_st" = paste(strsplit(fips, "")[[1]][1], strsplit(fips, "")[[1]][2], sep = ""), "fips_co" = paste(strsplit(fips, "")[[1]][3], strsplit(fips, "")[[1]][4], sep = "", strsplit(fips, "")[[1]][5]))
# 
# write.table(pair_coord, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/EPest_recordspair_fis.txt", quote = F, row.names = F)
pair_coord <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/EPest_recordspair_fis.txt")


#and ACRES of crop land
# Returns cropland ACRES per county-year from USDA NASS (1997-present).

# install.packages("rnassqs")
library(rnassqs)

# Get API key at https://quickstats.nass.usda.gov/api,
get_cropland_ha_census <- function(fips, years, pause = 0.5) {
  fips        <- sprintf("%05d", as.integer(fips))
  state_fips  <- substr(fips, 1, 2)
  county_fips <- substr(fips, 3, 5)

  # Ag Census only runs in years ending in 2 or 7, and only 1997+ is in the API
  census_years <- c(1997, 2002, 2007, 2012, 2017)
  valid_years  <- intersect(years, census_years)

  do.call(rbind, lapply(valid_years, function(yr) {
    dat <- tryCatch(
      nassqs(list(
        source_desc     = "CENSUS",
        short_desc      = "AG LAND, CROPLAND - ACRES",
        domain_desc     = "TOTAL",
        agg_level_desc  = "COUNTY",
        year            = yr,
        state_fips_code = state_fips,
        county_code     = county_fips
      )),
      error = function(e) {
        message(sprintf("fips %s, year %d failed: %s", fips, yr, e$message))
        NULL
      }
    )
    Sys.sleep(pause)

    if (is.null(dat) || nrow(dat) == 0) {
      return(data.frame(fips = fips, year = yr, cropland_ha = NA_real_))
    }

    # Value comes back as a formatted string (e.g. "12,345" or "(D)" if suppressed)
    acres <- suppressWarnings(as.numeric(gsub(",", "", dat$Value[1])))
    data.frame(fips = fips, year = yr, cropland_ha = round(acres * 0.404686, 1))
  }))
}

#Get acres of cropland for every 5 years 1992-2017
nassqs_auth(key = "5C7D44EE-9771-300D-B6CC-4A0396141D11")
crop_1992_18 <- do.call(rbind, lapply(pair_coord$fips, get_cropland_ha_census, years = 1992:2018))

# then join to her multi-year EPest kg and divide:
herb %>% left_join(crop, by = c("fips","year")) %>%
  mutate(kg_per_ha = kg / cropland_ha)
```

#PCA AFvaper model choosing with herbicide averaging and genetic distance
```{r}
#26 year avg of herbicide use 
#get fips for herbicide use and acres of cropland
#pair_coord <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/paircoord.txt")
#pair_coord <- pair_coord %>%
#  mutate(V5 = paste(V3, V4, sep = "_")) %>%
#  select(V5, V1, V2)

#fips is constant across our years so we can use 2010
# latlong2fips <- function(latitude, longitude) {
#   url <- "https://geo.fcc.gov/api/census/area?lat=%s&lon=%s&censusYear=2010&format=json"
#   url <- sprintf(url, latitude, longitude)
#   json <- RJSONIO::fromJSON(url)
#   as.character(json[["results"]][[1]][["county_fips"]])
# }
# 
# n <- nrow(pair_coord)
# pair_coord$fips <- rep(NA_character_, n)
# 
# for (i in 1:n){
#   pair_coord$fips[i] <- latlong2fips(pair_coord$V1[i], pair_coord$V2[i])
# }
# 
# pair_coord <- pair_coord %>%
#   rowwise() %>%
#   mutate("fips_st" = paste(strsplit(fips, "")[[1]][1], strsplit(fips, "")[[1]][2], sep = ""), "fips_co" = paste(strsplit(fips, "")[[1]][3], strsplit(fips, "")[[1]][4], sep = "", strsplit(fips, "")[[1]][5]))
# 
# write.table(pair_coord, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/EPest_recordspair_fis.txt", quote = F, row.names = F)
#
#
# pair_coord <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/EPest_recordspair_fis.txt")
# 
# target_compounds <- c("ATRAZINE", "FOMESAFEN", "GLYPHOSATE")
# 
# #Filter each EPest file to correct compounds and counties
# read_epest_year <- function(path, fips_keep) {
#   dat <- fread(path, sep = "\t")
#   dat <- dat[COMPOUND %in% target_compounds]
#   dat[, fips := sprintf("%02d%03d", as.integer(STATE_FIPS_CODE), as.integer(COUNTY_FIPS_CODE))]
#   dat <- dat[fips %in% fips_keep]
#   dat[, .(fips, year = YEAR, COMPOUND, kg = EPEST_HIGH_KG)]
# }
# 
# site_fips <- unique(pair_coord$fips)
# 
# herb_2013_2017 <- read_epest_year("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/EPest_recordsEPest_county_estimates_2013_2017_v2.txt", site_fips)
# remaining_years <- 1992:2012
# herb_1992_2012 <- rbindlist(lapply(remaining_years, function(yr) {
#   fpath <- sprintf("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/EPest_recordsEPest.county.estimates.%d.txt", yr)
#   if (!file.exists(fpath)) {
#     message("Missing file for year ", yr, ": ", fpath)
#     return(NULL)
#   }
#   read_epest_year(fpath, site_fips)
# }), fill = TRUE)
# 
# #combine all years
# herb_all <- rbindlist(list(herb_1992_2012, herb_2013_2017))
# write.table(herb_all, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/EPest_recordsherbicideuse_92_17.txt", quote = F, row.names = F)
# #compressed original EPest records as EPest_records.zip

#add in ACRES of crop land
# Returns cropland ACRES per county-year from USDA NASS (1997-present).

# install.packages("rnassqs")
library(rnassqs)

# # Get API key at https://quickstats.nass.usda.gov/api,
# get_cropland_acres_census <- function(fips, years, pause = 0.5) {
#   fips        <- sprintf("%05d", as.integer(fips))
#   state_fips  <- substr(fips, 1, 2)
#   county_fips <- substr(fips, 3, 5)
# 
#   # Ag Census only runs in years ending in 2 or 7, and only 1997+ is in the API
#   census_years <- c(1997, 2002, 2007, 2012, 2017)
#   valid_years  <- intersect(years, census_years)
# 
#   do.call(rbind, lapply(valid_years, function(yr) {
#     dat <- tryCatch(
#       nassqs(list(
#         source_desc     = "CENSUS",
#         short_desc      = "AG LAND, CROPLAND - ACRES",
#         domain_desc     = "TOTAL",
#         agg_level_desc  = "COUNTY",
#         year            = yr,
#         state_fips_code = state_fips,
#         county_code     = county_fips
#       )),
#       error = function(e) {
#         message(sprintf("fips %s, year %d failed: %s", fips, yr, e$message))
#         NULL
#       }
#     )
#     Sys.sleep(pause)
# 
#     if (is.null(dat) || nrow(dat) == 0) {
#       return(data.frame(fips = fips, year = yr, cropland_acres = NA_real_))
#     }
# 
#     # Value comes back as a formatted string (e.g. "12,345" or "(D)" if suppressed)
#     acres <- suppressWarnings(as.numeric(gsub(",", "", dat$Value[1])))
#     data.frame(fips = fips, year = yr, cropland_acres = round(acres * 0.404686, 1))
#   }))
# }
# 
# # #Get acres of cropland for every 5 years 1992-2017
# nassqs_auth(key = "5C7D44EE-9771-300D-B6CC-4A0396141D11")
# crop_1992_18 <- do.call(rbind, lapply(pair_coord$fips, get_cropland_acres_census, years = 1992:2018))
# #double checking all sites have some observations
# avg_acres <- crop_1992_18 %>%
#   group_by(fips) %>%
#   summarise("mean_acres" = mean(cropland_acres))
# no_census_obs <- avg_acres$fips[which(is.na(avg_acres$mean_acres))]
# crop_add <- do.call(rbind, lapply(no_census_obs, get_cropland_acres_census, years = 1992:2018))
# crop_1992_18 <- distinct(crop_1992_18, fips, year, .keep_all = TRUE)
# crop_add     <- distinct(crop_add, fips, year, .keep_all = TRUE)
# 
# #merge, then coalesce the two acreage columns into one
# crop_all <- merge(crop_1992_18, crop_add, by = c("fips", "year"), all = TRUE) %>%
#   mutate(cropland_acres = coalesce(cropland_acres.x, cropland_acres.y)) %>%
#   select(fips, year, cropland_acres)
# write.table(crop_all, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/EPest_recordsacres_ofcropland.txt", quote = F, row.names = F)
# crop_1992_18 <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/EPest_recordsacres_ofcropland.txt")
# 
# #then pull out how much land in the county in total
# library(readr)
# 
# # download of county land areas (static, not year-specific)
# download.file(
#   "https://www2.census.gov/geo/docs/maps-data/data/gazetteer/2020_Gazetteer/2020_Gaz_counties_national.zip",
#   destfile = "2020_Gaz_counties_national.zip"
# )
# unzip("2020_Gaz_counties_national.zip")
# 
# gaz <- read_tsv("2020_Gaz_counties_national.txt") %>%
#   transmute(
#     fips = sprintf("%05d", as.integer(GEOID)),
#     total_land_acres = ALAND / 4046.86   # ALAND is in square meters; 1 acre = 4046.86 m^2
#   )
# gaz$fips <- as.integer(gaz$fips)
# 
# # join to cropland data and compute the proportion
# cropland_with_pct <- crop_1992_18 %>%          # fips, year, cropland_acres
#   left_join(gaz, by = "fips") %>%
#   mutate(pct_cropland = cropland_acres / total_land_acres)
# 
# write.table(cropland_with_pct, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/EPest_recordspropor_ofcropland.txt", quote = F, row.names = F)

##run from here ------
# crop_1992_18 <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/EPest_recordspropor_ofcropland.txt")
# crop_1992_18 <- unique(crop_1992_18)
# crop_1992_18[, fips := sprintf("%05d", as.integer(fips))]
# 
# #cropland_ha = round(crop * 0.404686, 1))  # acres -> ha
# 
# assign_census_chunk <- function(yr) {
#   breaks <- c(-Inf, 1999.5, 2004.5, 2009.5, 2014.5, Inf)
#   labels <- c(1997, 2002, 2007, 2012, 2017)
#   as.numeric(as.character(cut(yr, breaks = breaks, labels = labels)))
# }
# 
# herb_all <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/EPest_recordsherbicideuse_92_17.txt")
# herb_all[, census_chunk := assign_census_chunk(year)]
# herb_all[, fips := sprintf("%05d", as.integer(fips))]
# crop_1992_18[, fips := sprintf("%05d", as.integer(fips))]
# 
# herb_joined <- herb_all %>%
#   left_join(crop_1992_18, by = c("fips", "census_chunk" = "year")) %>%
#   mutate(kg_peracre = kg / cropland_acres, fips = as.integer(fips))
# 
# pair_coord <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/EPest_recordspair_fis.txt")
# 
# #attach site names (V5)
# herb_site <- herb_joined %>%
#   left_join(pair_coord %>% select(V5, fips), by = "fips", relationship = "many-to-many")
#   # many-to-many is expected/fine here: each fips can map to multiple sites
#   # (e.g. 2_Nat and 2_Ag share fips 39063) and multiple years/compounds
# 
# cropland_2017 <- crop_1992_18 %>%
#   filter(year == 2017) %>%
#   mutate(fips = as.integer(fips)) %>%
#   select(fips, procropland_2017 = pct_cropland)
# 
# site_fips_lookup <- pair_coord %>% select(V5, fips)
# 
# # per-site, per-compound average across all years (1992-2017)
# site_compound_avg <- herb_site %>%
#   group_by(V5, COMPOUND) %>%
#   summarize(
#     mean_kg        = mean(kg, na.rm = TRUE),
#     mean_kg_per_acre = mean(kg_peracre, na.rm = TRUE),
#     mean_cropland_acres = mean(cropland_acres, na.rm = TRUE),
#     mean_procropland = mean(pct_cropland, na.rm = TRUE),
#     n_years        = sum(!is.na(kg)),
#     .groups = "drop"
#   ) %>%
#   left_join(site_fips_lookup, by = "V5") %>%
#   left_join(cropland_2017, by = "fips")
# 
# #drop all 7b_Nat bc this is the same county and I don't have it in any later things
# site_compound_avg <- filter(site_compound_avg, site_compound_avg$V5 != "7b_Nat")
# colnames(site_compound_avg) <- c("population", "compound", "mean_kg", "mean_kg_per_acre", "mean_cropland_acres", "mean_procropland", "n_years", "fips", "procropland_2017")
# 
# write.table(site_compound_avg, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/herbicide_and_cropland_avg.txt", row.names = F, quote = F)

#get climate at the pair level (pair mean: ag and nat)
clim <- read.table(file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/bioclim_pairstogetherAFvaper.txt", header = T)
clim_acrosshab <- clim %>%
  group_by(Pair_name) %>%
  summarise("B1_Ann_T" = mean(B1_Ann_T),"B2_Diurnal_Range" = mean(B2_Diurnal_Range),"B3_Isothermality"= mean(B3_Isothermality),"B4_T_Seasonality"= mean(B4_T_Seasonality),"B5_MaxT_Wrmst_Month"= mean(B5_MaxT_Wrmst_Month),"B6_MinT_Cldst_Month"= mean(B6_MinT_Cldst_Month),"B7_T_Ann_Range"= mean(B7_T_Ann_Range),"B8_T_Wettest_Qtr"= mean(B8_T_Wettest_Qtr),"B9_T_Driest_Qtr"= mean(B9_T_Driest_Qtr),"B10_T_Wrmst_Qtr"= mean(B10_T_Wrmst_Qtr),"B11_T_Cldst_Qtr"= mean(B11_T_Cldst_Qtr),"B12_Ann_Precip"= mean(B12_Ann_Precip),"B13_P_Wettest_Month"= mean(B13_P_Wettest_Month),"B14_P_Driest_Month"= mean(B14_P_Driest_Month),"B15_P_Seasonality"= mean(B15_P_Seasonality),"B16_P_Wettest_Qtr"= mean(B16_P_Wettest_Qtr),"B17_P_Driest_Qtr"= mean(B17_P_Driest_Qtr),"B18_P_Wrmst_Qtr"= mean(B18_P_Wrmst_Qtr),"B19_P_Cldst_Qtr"= mean(B19_P_Cldst_Qtr))

#------ run PCA/dimension reduction on this climate------
clim_acrosshab <- clim_acrosshab %>%
  column_to_rownames(var = "Pair_name")

png(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clim_afvaper/climPCA.png", sep = ""), width = 30, height = 30, units = "cm", res = 120)
print(biplot(prcomp(clim_acrosshab, scale = T), cex = c(0.5,0.4)))
dev.off()

pca <- prcomp(clim_acrosshab, scale = T)
summary(pca)

pca_scores <- as.data.frame(pca$x)
pca_scores$pair_name <- rownames(pca_scores)

#add in count data and get dif in pc scores between hab within same pair
p2_contrasts <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Afvaper_pairsfoundtogether.txt") 
for(c in 1:17){
  #dif between pair i and pair j
  varname <- paste0("dif_PC", c)

  p2_contrasts <- p2_contrasts %>%
    rowwise() %>%
    mutate(val1 = pca_scores[[which(pca_scores$pair_name == V1), c]], val2 = pca_scores[[which(pca_scores$pair_name == V2), c]], !!varname := abs(val1 - val2)) %>%
    ungroup() %>%
    select(-val1, -val2)
}

#get the difference in the proportion of ag land, kg_acre and kg of herbicide use for 26 yr avg (1992-2017) for a county (mean ag + nat at a site, site i - j)
#ag site herbicide use
site_compound_avg <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/herbicide_and_cropland_avg.txt")
site_compound_avg <- site_compound_avg %>%
  rowwise() %>%
  mutate("pair_num" = str_split(population, "_")[[1]][1], "Env" = str_split(population, "_")[[1]][2], "Env" = case_when(
    Env == "Ag" ~ "AG",
    Env == "Nat" ~ "NAT"
  )) %>%
    ungroup()
site_compound_avg$Pair_name <- paste("pair_", site_compound_avg$pair_num, sep = "")
p2_contrasts <- p2_contrasts %>%
  rowwise() %>%
  mutate("dif_proag_26avg" = abs(mean(site_compound_avg$mean_procropland[which(site_compound_avg$Pair_name == V1)]) - mean(site_compound_avg$mean_procropland[which(site_compound_avg$Pair_name == V2)])), "dif_proag_2017" = abs(mean(site_compound_avg$procropland_2017[which(site_compound_avg$Pair_name == V1)]) - mean(site_compound_avg$procropland_2017[which(site_compound_avg$Pair_name == V2)]))) %>%
  mutate("ag_gly_kg_acre" = abs(site_compound_avg$mean_kg_per_acre[which(site_compound_avg$Pair_name == V1 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "GLYPHOSATE")] - site_compound_avg$mean_kg_per_acre[which(site_compound_avg$Pair_name == V2 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "GLYPHOSATE")]), "ag_atrazine_kg_acre" = abs(site_compound_avg$mean_kg_per_acre[which(site_compound_avg$Pair_name == V1 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "ATRAZINE")] - site_compound_avg$mean_kg_per_acre[which(site_compound_avg$Pair_name == V2 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "ATRAZINE")]), "ag_fomesafen_kg_acre" = abs(site_compound_avg$mean_kg_per_acre[which(site_compound_avg$Pair_name == V1 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "FOMESAFEN")] - site_compound_avg$mean_kg_per_acre[which(site_compound_avg$Pair_name == V2 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "FOMESAFEN")]), "ag_gly_kg" = abs(site_compound_avg$mean_kg[which(site_compound_avg$Pair_name == V1 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "GLYPHOSATE")] - site_compound_avg$mean_kg[which(site_compound_avg$Pair_name == V2 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "GLYPHOSATE")]), "ag_atrazine_kg" = abs(site_compound_avg$mean_kg[which(site_compound_avg$Pair_name == V1 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "ATRAZINE")] - site_compound_avg$mean_kg[which(site_compound_avg$Pair_name == V2 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "ATRAZINE")]), "ag_fomesafen_kg" = abs(site_compound_avg$mean_kg[which(site_compound_avg$Pair_name == V1 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "FOMESAFEN")] - site_compound_avg$mean_kg[which(site_compound_avg$Pair_name == V2 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "FOMESAFEN")])) %>%
  ungroup()

#add in pairwise distance bw i and j
dist <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Pairwise_distances.txt")
dist <- dist %>%
  rowwise() %>%
  mutate("p1" = paste("pair_", str_split(pop1, "_")[[1]][2], sep = ""), "e1" = str_split(pop1, "_")[[1]][3], "p2" = paste("pair_", str_split(pop2, "_")[[1]][2], sep =""),"e2" = str_split(pop2, "_")[[1]][3], "pop_contrast" = case_when(
    e1 == e2 ~ "same_env",
    e1 != e2 ~ "ag-nat"
  )) %>%
  ungroup()
dist$p1 <- factor(dist$p1, c("pair_1","pair_2",  "pair_3",  "pair_4",  "pair_5",  "pair_6",  "pair_7",  "pair_8",  "pair_9",  "pair_10", "pair_11", "pair_12", "pair_13","pair_14", "pair_15", "pair_16", "pair_17"))
dist$p2 <- factor(dist$p2, c("pair_1","pair_2",  "pair_3",  "pair_4",  "pair_5",  "pair_6",  "pair_7",  "pair_8",  "pair_9",  "pair_10", "pair_11", "pair_12", "pair_13","pair_14", "pair_15", "pair_16", "pair_17"))
p2_contrasts$mean_pairdist <- 0

for(p in 1:length(p2_contrasts$V1)){
  ind <- which((dist$p1 == p2_contrasts$V1[p] & dist$p2 == p2_contrasts$V2[p]) | (dist$p2 == p2_contrasts$V1[p] & dist$p1 == p2_contrasts$V2[p]))
  ind_use <- filter(dist[ind,])
  #have to average over ag and nat 
  p2_contrasts$mean_pairdist[p] <- sum(ind_use$Geo_km)/length(ind_use$pop1)
}

p2_contrasts$Contrast <- paste(p2_contrasts$V1, p2_contrasts$V2, sep = ":")

#Add in the genetic distance between pair i and j 
fst_wc <- tibble("pop1" = character(), "pop2" = character(), "chromosome" = character(), "window_pos_1" = numeric(), "window_pos_2" = numeric(), "avg_wc_fst" = numeric(), "no_snps" = numeric(),	"wc_fst_a" = numeric(),	"wc_fst_b" = numeric(),	"wc_fst_c" = numeric())

for(p in 1:16){
  pix_file <- fread(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/pixy_fst_50Mb_Scaffold_", p, ".txt", sep = ""), sep="\t", header=T)
  fst_wc <- fst_wc %>%
    add_row("pop1" = pix_file$pop1, "pop2" = pix_file$pop2, "chromosome" = pix_file$chromosome, "window_pos_1" = pix_file$window_pos_1, "window_pos_2" = pix_file$window_pos_2, "avg_wc_fst" = pix_file$avg_wc_fst, "no_snps" = pix_file$no_snps,	"wc_fst_a" = pix_file$wc_fst_a,	"wc_fst_b" = pix_file$wc_fst_b,	"wc_fst_c" = pix_file$wc_fst_c)
}

#average over all chromosomes and across p1-ag:p2-ag,p1-nat:p2-nat habitats within a pair
fst_wc_avg <- fst_wc %>%
  rowwise() %>%
  mutate("p1" = str_split(pop1, "_")[[1]][2], "p2" = str_split(pop2, "_")[[1]][2]) %>%
  ungroup() %>%
  mutate("contrast" = paste(p1, "_", p2, sep = "")) %>%
  group_by(contrast) %>%
  summarise("wcfst_allscaf" = sum(wc_fst_a) / (sum(wc_fst_a) + sum(wc_fst_b) + sum(wc_fst_c)), "pop1" = paste("pair_", first(p1), sep = ""), "pop2" = paste("pair_", first(p2), sep = "")) %>%
  ungroup()

p2_contrasts$mean_fst <- 0

for(p in 1:length(p2_contrasts$V1)){
  ind <- which((fst_wc_avg$pop1 == p2_contrasts$V1[p] & fst_wc_avg$pop2 == p2_contrasts$V2[p]) | (fst_wc_avg$pop2 == p2_contrasts$V1[p] & fst_wc_avg$pop1 == p2_contrasts$V2[p]))
  p2_contrasts$mean_fst[p] <- unlist(fst_wc_avg$wcfst_allscaf[ind])
}


#---- Run univariate model on genetic distance----------
univariate_fst <- ggplot(p2_contrasts, aes(x = mean_fst, y = count_onlyparallel)) + geom_point() + geom_smooth(method = lm) + labs(title = "count of pair cooccurance in parallel category vs population differentiation") + theme_bw() + theme(axis.text = element_text(size = 20))
ggsave(univariate_fst, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clim_afvaper/Afvaper_parallel_pairstogether_related.png", device = "png", height = 6, width = 12, units = "in")
fit <- glm(p2_contrasts$count_onlyparallel ~ p2_contrasts$mean_fst)
summary(fit)
pseudo_r2 <- 1 - (fit$deviance / fit$null.deviance)


#----- Dimenstionality reduction on herbicide use ----------
#only using the herbicide data from the ag site
herb_ag_pca <- site_compound_avg %>%
  filter(Env == "AG") 

#make only one obs per pop
herb_ag_pca_wide <- herb_ag_pca %>%
  pivot_wider(
    id_cols = c(population, pair_num, Env, Pair_name),
    names_from = compound,
    values_from = c(mean_kg, mean_kg_per_acre, mean_cropland_acres, mean_procropland, procropland_2017),
    names_glue = "{.value}_{compound}"
  )
herb_ag_pca_wide$mean_cropland_acres <- herb_ag_pca_wide$mean_cropland_acres_ATRAZINE
herb_ag_pca_wide$mean_procropland <- herb_ag_pca_wide$mean_procropland_ATRAZINE
herb_ag_pca_wide$procropland_2017 <- herb_ag_pca_wide$procropland_2017_ATRAZINE
herb_ag_pca_wide <- herb_ag_pca_wide[-seq(11,19, by = 1)]

herb_ag_pca <- herb_ag_pca_wide %>%
  column_to_rownames(var = "Pair_name") 
herb_ag_pca <- herb_ag_pca[,seq(4,12, by = 1)]

png(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clim_afvaper/herbPCA.png", sep = ""), width = 30, height = 30, units = "cm", res = 120)
print(biplot(prcomp(herb_ag_pca, scale = T), cex = c(0.5,0.4)))
dev.off()

herb_pca <- prcomp(herb_ag_pca, scale = T)
summary(herb_pca)

herb_pca_scores <- as.data.frame(herb_pca$x)
herb_pca_scores$pair_name <- rownames(herb_pca_scores)

#now adding the dim red herb to the others
p2_contrasts_wPCherb <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Afvaper_pairsfoundtogether.txt") 
for(c in 1:9){
  #dif between pair i and pair j
  #herbicide
  varname <- paste0("herbdif_PC", c)

  p2_contrasts_wPCherb <- p2_contrasts_wPCherb %>%
    rowwise() %>%
    mutate(val1 = herb_pca_scores[[which(herb_pca_scores$pair_name == V1), c]], val2 = herb_pca_scores[[which(herb_pca_scores$pair_name == V2), c]], !!varname := abs(val1 - val2)) %>%
    ungroup() %>%
    select(-val1, -val2)
  
  #climate
  varname <- paste0("dif_PC", c)

  p2_contrasts_wPCherb <- p2_contrasts_wPCherb %>%
    rowwise() %>%
    mutate(val1 = pca_scores[[which(pca_scores$pair_name == V1), c]], val2 = pca_scores[[which(pca_scores$pair_name == V2), c]], !!varname := abs(val1 - val2)) %>%
    ungroup() %>%
    select(-val1, -val2)
}

#get the difference in the proportion of ag land within 1000 km of the site (mean ag + nat at a site, site i - j)
#and difference in shannon div index
p2_contrasts_wPCherb <- p2_contrasts_wPCherb %>%
  rowwise() %>%
  mutate("dif_proag_26avg" = abs(mean(site_compound_avg$mean_procropland[which(site_compound_avg$Pair_name == V1)]) - mean(site_compound_avg$mean_procropland[which(site_compound_avg$Pair_name == V2)])), "dif_proag_2017" = abs(mean(site_compound_avg$procropland_2017[which(site_compound_avg$Pair_name == V1)]) - mean(site_compound_avg$procropland_2017[which(site_compound_avg$Pair_name == V2)]))) %>%
  mutate("ag_gly_kg_acre" = abs(site_compound_avg$mean_kg_per_acre[which(site_compound_avg$Pair_name == V1 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "GLYPHOSATE")] - site_compound_avg$mean_kg_per_acre[which(site_compound_avg$Pair_name == V2 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "GLYPHOSATE")]), "ag_atrazine_kg_acre" = abs(site_compound_avg$mean_kg_per_acre[which(site_compound_avg$Pair_name == V1 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "ATRAZINE")] - site_compound_avg$mean_kg_per_acre[which(site_compound_avg$Pair_name == V2 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "ATRAZINE")]), "ag_fomesafen_kg_acre" = abs(site_compound_avg$mean_kg_per_acre[which(site_compound_avg$Pair_name == V1 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "FOMESAFEN")] - site_compound_avg$mean_kg_per_acre[which(site_compound_avg$Pair_name == V2 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "FOMESAFEN")]), "ag_gly_kg" = abs(site_compound_avg$mean_kg[which(site_compound_avg$Pair_name == V1 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "GLYPHOSATE")] - site_compound_avg$mean_kg[which(site_compound_avg$Pair_name == V2 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "GLYPHOSATE")]), "ag_atrazine_kg" = abs(site_compound_avg$mean_kg[which(site_compound_avg$Pair_name == V1 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "ATRAZINE")] - site_compound_avg$mean_kg[which(site_compound_avg$Pair_name == V2 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "ATRAZINE")]), "ag_fomesafen_kg" = abs(site_compound_avg$mean_kg[which(site_compound_avg$Pair_name == V1 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "FOMESAFEN")] - site_compound_avg$mean_kg[which(site_compound_avg$Pair_name == V2 & site_compound_avg$Env == "AG" & site_compound_avg$compound == "FOMESAFEN")])) %>%
  ungroup()

#Add in the distance between pair i and j (averaged over the habitat types distances)
p2_contrasts_wPCherb$mean_pairdist <- 0

for(p in 1:length(p2_contrasts_wPCherb$V1)){
  ind <- which((dist$p1 == p2_contrasts_wPCherb$V1[p] & dist$p2 == p2_contrasts_wPCherb$V2[p]) | (dist$p2 == p2_contrasts_wPCherb$V1[p] & dist$p1 == p2_contrasts_wPCherb$V2[p]))
  ind_use <- filter(dist[ind,])
  #have to average over ag and nat 
  p2_contrasts_wPCherb$mean_pairdist[p] <- sum(ind_use$Geo_km)/length(ind_use$pop1)
}

p2_contrasts_wPCherb$Contrast <- paste(p2_contrasts_wPCherb$V1, p2_contrasts_wPCherb$V2, sep = ":")

#Add in the genetic distance between pair i and j (averaged over the habitat types distances)
p2_contrasts_wPCherb$mean_fst <- 0
for(p in 1:length(p2_contrasts_wPCherb$V1)){
  ind <- which((fst_wc_avg$pop1 == p2_contrasts_wPCherb$V1[p] & fst_wc_avg$pop2 == p2_contrasts_wPCherb$V2[p]) | (fst_wc_avg$pop2 == p2_contrasts_wPCherb$V1[p] & fst_wc_avg$pop1 == p2_contrasts_wPCherb$V2[p]))
  p2_contrasts_wPCherb$mean_fst[p] <- unlist(fst_wc_avg$wcfst_allscaf[ind])
}

#------- Now have all the variables, start comparing models---------------------
everything <- glm(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + dif_proag_26avg + ag_gly_kg_acre + ag_atrazine_kg_acre + ag_fomesafen_kg_acre + ag_gly_kg + ag_atrazine_kg + ag_fomesafen_kg + dif_proag_26avg:ag_gly_kg_acre + dif_proag_26avg:ag_atrazine_kg_acre + dif_proag_26avg:ag_fomesafen_kg_acre + mean_fst, data = p2_contrasts, family = poisson)
summary(everything) 

#checking if overdispersion
library(AER) 
dispersiontest(everything) #dispersion greater than 1, disp = 12.55132 


#---------------neg bin with mean fst for only bioclim dim reduction ---------------------------
library(MASS)
everything <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + dif_proag_26avg + dif_proag_2017 + ag_gly_kg_acre + ag_atrazine_kg_acre + ag_fomesafen_kg_acre + ag_gly_kg + ag_atrazine_kg + ag_fomesafen_kg + dif_proag_26avg:ag_gly_kg_acre + dif_proag_26avg:ag_atrazine_kg_acre + dif_proag_26avg:ag_fomesafen_kg_acre + mean_fst, data = p2_contrasts)
summary(everything) 
#AIC: 1437.6 Sig: PC1,PC2, mean fst

herb <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + ag_gly_kg_acre + ag_atrazine_kg_acre + ag_fomesafen_kg_acre + ag_gly_kg + ag_atrazine_kg + ag_fomesafen_kg + dif_proag_26avg:ag_gly_kg_acre + dif_proag_26avg:ag_atrazine_kg_acre + dif_proag_26avg:ag_fomesafen_kg_acre + mean_fst, data = p2_contrasts) 
summary(herb)
#AIC: 1434.3 Sig: PC1+2, mean fst 

herb_kgacre <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + ag_gly_kg_acre + ag_atrazine_kg_acre + ag_fomesafen_kg_acre + dif_proag_26avg:ag_gly_kg_acre + dif_proag_26avg:ag_atrazine_kg_acre + dif_proag_26avg:ag_fomesafen_kg_acre + mean_fst, data = p2_contrasts) 
summary(herb_kgacre) 
#AIC: 1429.6 Sig: PC1+2, mean_fst

land <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + dif_proag_26avg + dif_proag_2017 + mean_fst, data = p2_contrasts)
summary(land) 
#AIC: 1423.3 Sig: pc1+2 AND mean_fst

none <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + mean_fst, data = p2_contrasts)
summary(none) 
#AIC: 1419.5 Sig: all except interaction of PC1+2

no_interaction <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + mean_fst, data = p2_contrasts)
summary(no_interaction)
#AIC: 1417.5 Sig: All

no_interaction2 <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_proag_26avg + mean_fst, data = p2_contrasts)
summary(no_interaction2)
#AIC: 1419.5 Sig: all except proportion of agriculture

no_interaction3 <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_proag_2017 + mean_fst, data = p2_contrasts)
summary(no_interaction3)
#AIC: 1419.5 Sig: all except proportion of agriculture

herbkgacre_land <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + ag_gly_kg_acre + ag_atrazine_kg_acre + ag_fomesafen_kg_acre + mean_fst, data = p2_contrasts)
summary(herbkgacre_land) 
#AIC: 1423.8 Sig: PC1+2, mean_fst

#---------------neg bin with mean fst for bioclim AND herbicide dim reduction ---------------------------

everything_herbred <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + dif_proag_26avg + dif_proag_2017 + herbdif_PC1 + herbdif_PC2 + herbdif_PC1:herbdif_PC2 + mean_fst, data = p2_contrasts_wPCherb)
summary(everything_herbred) 
#AIC: 1427.9 Sig: PC1+2, mean_fst

herbred_noint <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + herbdif_PC1 + herbdif_PC2 + mean_fst, data = p2_contrasts_wPCherb)
summary(herbred_noint)
#AIC: 1420.4 Sig: PC1+2, mean_fst

noint_herbred <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_proag_26avg + herbdif_PC1 + herbdif_PC2 + mean_fst, data = p2_contrasts_wPCherb)
summary(noint_herbred) 
#AIC: 1422 Sig: PC1+2, mean_fst
#dif_proag_26avg is better than the model with dif_proag_2017 (rsq .385 vs .40)

min_herbred <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + herbdif_PC1 + mean_fst, data = p2_contrasts_wPCherb)
summary(min_herbred) 
#AIC: 1419.3 Sig: PC1+2, mean fst


#----- trying with mean pair distance to see if that or the average is why environment dropped out --------
everything_herbred_dist <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + dif_proag_26avg + dif_proag_2017 + herbdif_PC1 + herbdif_PC2 + herbdif_PC1:herbdif_PC2 + mean_pairdist, data = p2_contrasts_wPCherb)
summary(everything_herbred_dist) 
#AIC: 1425.7 Sig: PC1+2, mean_pairdist

everything_dist <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_PC1:dif_PC2 + dif_proag_26avg + dif_proag_2017 + ag_gly_kg_acre + ag_atrazine_kg_acre + ag_fomesafen_kg_acre + ag_gly_kg + ag_atrazine_kg + ag_fomesafen_kg + dif_proag_26avg:ag_gly_kg_acre + dif_proag_26avg:ag_atrazine_kg_acre + dif_proag_26avg:ag_fomesafen_kg_acre + mean_pairdist, data = p2_contrasts)
summary(everything_dist) 
#AIC: 1436.3 Sig: PC1+2, mean_pairdist


#----- Looking into best model -----------------

#decided on no interaction herbicide reduction model
noint_herbred <- glm.nb(count_onlyparallel ~ dif_PC1 + dif_PC2 + dif_proag_26avg + herbdif_PC1 + herbdif_PC2 + mean_fst, data = p2_contrasts_wPCherb)
summary(noint_herbred) 

pred_data <- data.frame(
  dif_PC2 = seq(min(p2_contrasts_wPCherb$dif_PC2, na.rm = TRUE),
                max(p2_contrasts_wPCherb$dif_PC2, na.rm = TRUE),
                length.out = 100),
  dif_PC1 = mean(p2_contrasts_wPCherb$dif_PC1, na.rm = TRUE),
  dif_proag_26avg = mean(p2_contrasts_wPCherb$dif_proag_26avg, na.rm = TRUE),
  herbdif_PC1 = mean(p2_contrasts_wPCherb$herbdif_PC1, na.rm = TRUE),
  herbdif_PC2 = mean(p2_contrasts_wPCherb$herbdif_PC2, na.rm = TRUE),
  mean_fst = mean(p2_contrasts_wPCherb$mean_fst, na.rm = TRUE)
)

# Get predictions on the response scale, with standard errors
pred <- predict(noint_herbred, newdata = pred_data, type = "link", se.fit = TRUE)
pred_data$fit <- exp(pred$fit)
pred_data$lower <- exp(pred$fit - 1.96 * pred$se.fit)
pred_data$upper <- exp(pred$fit + 1.96 * pred$se.fit)

climplot_PC2 <- ggplot() + geom_point(data = p2_contrasts_wPCherb, aes(x = dif_PC2, y = count_onlyparallel)) + geom_ribbon(data = pred_data, aes(x = dif_PC2, ymin = lower, ymax = upper), alpha = 0.2) + geom_line(data = pred_data, aes(x = dif_PC2, y = fit), color = "blue", linewidth = 1) + labs(title = "count of pair cooccurance in parallel category vs dif_PC2 (other variables are held at their mean value", x = "PC2,i - PC2,j", y = "# parallel windows") + theme_bw() + theme(axis.text = element_text(size = 20))
ggsave(climplot_PC2, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clim_afvaper/Afvaper_parallel_pairstogether_env_final_PC2.png", device = "png", height = 6, width = 12, units = "in")

#do for all sig variables
pred2_data <- data.frame(
  dif_PC1 = seq(min(p2_contrasts_wPCherb$dif_PC1, na.rm = TRUE),
                max(p2_contrasts_wPCherb$dif_PC1, na.rm = TRUE),
                length.out = 100),
  dif_PC2 = mean(p2_contrasts_wPCherb$dif_PC2, na.rm = TRUE),
  dif_proag_26avg = mean(p2_contrasts_wPCherb$dif_proag_26avg, na.rm = TRUE),
  herbdif_PC1 = mean(p2_contrasts_wPCherb$herbdif_PC1, na.rm = TRUE),
  herbdif_PC2 = mean(p2_contrasts_wPCherb$herbdif_PC2, na.rm = TRUE),
  mean_fst = mean(p2_contrasts_wPCherb$mean_fst, na.rm = TRUE)
)

# Get predictions on the response scale, with standard errors
pred2 <- predict(noint_herbred, newdata = pred2_data, type = "link", se.fit = TRUE)
pred2_data$fit <- exp(pred2$fit)
pred2_data$lower <- exp(pred2$fit - 1.96 * pred2$se.fit)
pred2_data$upper <- exp(pred2$fit + 1.96 * pred2$se.fit)

climplot_PC1 <- ggplot() + geom_point(data = p2_contrasts_wPCherb, aes(x = dif_PC1, y = count_onlyparallel)) + geom_ribbon(data = pred2_data, aes(x = dif_PC1, ymin = lower, ymax = upper), alpha = 0.2) + geom_line(data = pred2_data, aes(x = dif_PC1, y = fit), color = "blue", linewidth = 1) + labs(title = "Difference in PC1 for pop i,j", x = "PC1,i - PC1,j", y = "# parallel windows") + theme_bw() + theme(axis.text = element_text(size = 20))

pred3_data <- data.frame(
  dif_PC1 = mean(p2_contrasts_wPCherb$dif_PC1, na.rm = TRUE),
  dif_PC2 = mean(p2_contrasts_wPCherb$dif_PC2, na.rm = TRUE),
  dif_proag_26avg = seq(min(p2_contrasts_wPCherb$dif_proag_26avg, na.rm = TRUE),
                max(p2_contrasts_wPCherb$dif_proag_26avg, na.rm = TRUE),
                length.out = 100),
  herbdif_PC1 = mean(p2_contrasts_wPCherb$herbdif_PC1, na.rm = TRUE),
  herbdif_PC2 = mean(p2_contrasts_wPCherb$herbdif_PC2, na.rm = TRUE),
  mean_fst = mean(p2_contrasts_wPCherb$mean_fst, na.rm = TRUE)
)

# Get predictions on the response scale, with standard errors
pred3 <- predict(noint_herbred, newdata = pred3_data, type = "link", se.fit = TRUE)
pred3_data$fit <- exp(pred3$fit)
pred3_data$lower <- exp(pred3$fit - 1.96 * pred3$se.fit)
pred3_data$upper <- exp(pred3$fit + 1.96 * pred3$se.fit)

climplot_ag <- ggplot() + geom_point(data = p2_contrasts_wPCherb, aes(x = dif_proag_26avg, y = count_onlyparallel)) + geom_ribbon(data = pred3_data, aes(x = dif_proag_26avg, ymin = lower, ymax = upper), alpha = 0.2) + geom_line(data = pred3_data, aes(x = dif_proag_26avg, y = fit), color = "blue", linewidth = 1) + labs(title = "Difference in % ag in county for pop i,j", x = "% ag in county", y = "# parallel windows") + theme_bw() + theme(axis.text = element_text(size = 20))

pred4_data <- data.frame(
  herbdif_PC1 = seq(min(p2_contrasts_wPCherb$herbdif_PC1, na.rm = TRUE),
                max(p2_contrasts_wPCherb$herbdif_PC1, na.rm = TRUE),
                length.out = 100),
  dif_PC2 = mean(p2_contrasts_wPCherb$dif_PC2, na.rm = TRUE),
  dif_proag_26avg = mean(p2_contrasts_wPCherb$dif_proag_26avg, na.rm = TRUE),
  dif_PC1 = mean(p2_contrasts_wPCherb$dif_PC1, na.rm = TRUE),
  herbdif_PC2 = mean(p2_contrasts_wPCherb$herbdif_PC2, na.rm = TRUE),
  mean_fst = mean(p2_contrasts_wPCherb$mean_fst, na.rm = TRUE)
)

# Get predictions on the response scale, with standard errors
pred4 <- predict(noint_herbred, newdata = pred4_data, type = "link", se.fit = TRUE)
pred4_data$fit <- exp(pred4$fit)
pred4_data$lower <- exp(pred4$fit - 1.96 * pred4$se.fit)
pred4_data$upper <- exp(pred4$fit + 1.96 * pred4$se.fit)

climplot_herbPC1 <- ggplot() + geom_point(data = p2_contrasts_wPCherb, aes(x = herbdif_PC1, y = count_onlyparallel)) + geom_ribbon(data = pred4_data, aes(x = herbdif_PC1, ymin = lower, ymax = upper), alpha = 0.2) + geom_line(data = pred4_data, aes(x = herbdif_PC1, y = fit), color = "blue", linewidth = 1) + labs(title = "Difference in herbPC1 for pop i,j", x = "herb_PC1,i - herb_PC1,j", y = "# parallel windows") + theme_bw() + theme(axis.text = element_text(size = 20))

pred5_data <- data.frame(
  herbdif_PC2 = seq(min(p2_contrasts_wPCherb$herbdif_PC2, na.rm = TRUE),
                max(p2_contrasts_wPCherb$herbdif_PC2, na.rm = TRUE),
                length.out = 100),
  dif_PC2 = mean(p2_contrasts_wPCherb$dif_PC2, na.rm = TRUE),
  dif_proag_26avg = mean(p2_contrasts_wPCherb$dif_proag_26avg, na.rm = TRUE),
  dif_PC1 = mean(p2_contrasts_wPCherb$dif_PC1, na.rm = TRUE),
  herbdif_PC1 = mean(p2_contrasts_wPCherb$herbdif_PC1, na.rm = TRUE),
  mean_fst = mean(p2_contrasts_wPCherb$mean_fst, na.rm = TRUE)
)

# Get predictions on the response scale, with standard errors
pred5 <- predict(noint_herbred, newdata = pred5_data, type = "link", se.fit = TRUE)
pred5_data$fit <- exp(pred5$fit)
pred5_data$lower <- exp(pred5$fit - 1.96 * pred5$se.fit)
pred5_data$upper <- exp(pred5$fit + 1.96 * pred5$se.fit)

climplot_herbPC2 <- ggplot() + geom_point(data = p2_contrasts_wPCherb, aes(x = herbdif_PC2, y = count_onlyparallel)) + geom_ribbon(data = pred5_data, aes(x = herbdif_PC2, ymin = lower, ymax = upper), alpha = 0.2) + geom_line(data = pred5_data, aes(x = herbdif_PC2, y = fit), color = "blue", linewidth = 1) + labs(title = "Difference in herbPC2 for pop i,j", x = "herb_PC1,i - herb_PC1,j", y = "# parallel windows") + theme_bw() + theme(axis.text = element_text(size = 20))

pred6_data <- data.frame(
  mean_fst = seq(min(p2_contrasts_wPCherb$mean_fst, na.rm = TRUE),
                max(p2_contrasts_wPCherb$mean_fst, na.rm = TRUE),
                length.out = 100),
  dif_PC1 = mean(p2_contrasts_wPCherb$dif_PC1, na.rm = TRUE),
  dif_PC2 = mean(p2_contrasts_wPCherb$dif_PC2, na.rm = TRUE),
  dif_proag_26avg = mean(p2_contrasts_wPCherb$dif_proag_26avg, na.rm = TRUE),
  herbdif_PC1 = mean(p2_contrasts_wPCherb$herbdif_PC1, na.rm = TRUE),
  herbdif_PC2 = mean(p2_contrasts_wPCherb$herbdif_PC2, na.rm = TRUE)
)

# Get predictions on the response scale, with standard errors
pred6 <- predict(noint_herbred, newdata = pred6_data, type = "link", se.fit = TRUE)
pred6_data$fit <- exp(pred6$fit)
pred6_data$lower <- exp(pred6$fit - 1.96 * pred6$se.fit)
pred6_data$upper <- exp(pred6$fit + 1.96 * pred6$se.fit)

climplot_fst <- ggplot() + geom_point(data = p2_contrasts_wPCherb, aes(x = mean_fst, y = count_onlyparallel)) + geom_ribbon(data = pred6_data, aes(x = mean_fst, ymin = lower, ymax = upper), alpha = 0.2) + geom_line(data = pred6_data, aes(x = mean_fst, y = fit), color = "blue", linewidth = 1) + labs(title = "Difference in fst for pop i,j", x = "Fsti - Fstj", y = "# parallel windows") + theme_bw() + theme(axis.text = element_text(size = 20))

#facet plot
library(patchwork)
all <- (climplot_PC1 + climplot_PC2) / (climplot_ag + climplot_fst) / (climplot_herbPC1 + climplot_herbPC2)
ggsave(all, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clim_afvaper/Afvaper_parallel_pairstogether_env_final_others.png", device = "png", height = 12, width = 12, units = "in")

#get rsq of this model
pseudo_r2 <- 1 - (noint_herbred$deviance / noint_herbred$null.deviance) #.404
pval <- summary(noint_herbred)$coefficients[2, 4]
```

#compare with cmh and xpehh hits
```{r}
sigwin <- readRDS("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/afvapr/afvapr_sigwindows.txt")
alleigen_sig_windows <- as.vector(unlist(sigwin, use.names = FALSE))
alleigen_sig_windows <- as.data.table(alleigen_sig_windows)
colnames(alleigen_sig_windows) <- "window_id"
alleigen_sig_windows <- alleigen_sig_windows %>%
  rowwise() %>%
  mutate("sep" = str_split(window_id, ":"), "Chr" = as.numeric(str_split(sep[[1]][1], "_")[[1]][2]), "pos" = str_split(window_id, ":")[[1]][2], "start" = as.numeric(str_split(pos, "-")[[1]][1]), "end" = as.numeric(str_split(pos, "-")[[1]][2]), "length" = end - start, "center_snp" = round((end-start)/2)) %>%
  select(window_id, Chr, start, end, length, center_snp) %>%
  ungroup()

#par = eig 1, multipar = sig on eig 2-4
parallel_sigwin <- alleigen_sig_windows[1:1049,]
multiparallel_sigwin <- alleigen_sig_windows[1050:2066,]
all_parallel <- alleigen_sig_windows[1:2066,]
all_parallel <- all_parallel %>%
  select(Chr, start, end)

#output table for 10kb window overlap check
write.table(all_parallel, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/afvaperhits_dcgm.bed", quote = F, col.names = F, row.names = F)

#check for overlaps with cmh and xpehh
xpehh_clumpsites <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/xpehh_sigsites.txt")
cmh_clumpsites <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clumphits_dcgm.txt")

overlap_sites <- tibble("par_sites" = length(parallel_sigwin$window_id), "multipar_sites" = length(multiparallel_sigwin$window_id), "par_xpehh_count" = 0, "par_cmh_count" = 0, "multipar_xpehh_count" = 0, "multipar_cmh_count" = 0, "par_xpehh" = c("sites"), "par_cmh" = c("sites"), "multipar_xpehh" = c("sites"), "multipar_cmh" = c("sites"))

for(i in 1:length(xpehh_clumpsites$CHROM)){
  par_site <- which(parallel_sigwin$Chr == xpehh_clumpsites$CHROM[i] & parallel_sigwin$start <= xpehh_clumpsites$POS[i] & parallel_sigwin$end >= xpehh_clumpsites$POS[i])
  multipar_site <- which(multiparallel_sigwin$Chr == xpehh_clumpsites$CHROM[i] & multiparallel_sigwin$start <= xpehh_clumpsites$POS[i] & multiparallel_sigwin$end >= xpehh_clumpsites$POS[i])
  if(length(par_site) > 0){
    overlap_sites$par_xpehh_count[1] <- overlap_sites$par_xpehh_count[1] + 1
    overlap_sites$par_xpehh <- as.list(overlap_sites$par_xpehh)
    overlap_sites$par_xpehh[[1]] <- c(overlap_sites$par_xpehh[[1]], as.character(xpehh_clumpsites$locus_id[i]))
    if(length(multipar_site) > 0){
      overlap_sites$multipar_xpehh_count[1] <- overlap_sites$multipar_xpehh_count[1] + 1
      overlap_sites$multipar_xpehh <- as.list(overlap_sites$multipar_xpehh)
      overlap_sites$multipar_xpehh[[1]] <- c(overlap_sites$multipar_xpehh[[1]], as.character(xpehh_clumpsites$locus_id[i]))
    }
  }
}

for(i in 1:length(cmh_clumpsites$CHROM)){
  par_site <- which(parallel_sigwin$Chr == cmh_clumpsites$CHROM[i] & parallel_sigwin$start <= cmh_clumpsites$BP[i] & parallel_sigwin$end >= cmh_clumpsites$BP[i])
  multipar_site <- which(multiparallel_sigwin$Chr == cmh_clumpsites$CHROM[i] & multiparallel_sigwin$start <= cmh_clumpsites$BP[i] & multiparallel_sigwin$end >= cmh_clumpsites$BP[i])
  if(length(par_site) > 0){
    overlap_sites$par_cmh_count[1] <- overlap_sites$par_cmh_count[1] + 1
    overlap_sites$par_cmh <- as.list(overlap_sites$par_cmh)
    overlap_sites$par_cmh[[1]] <- c(overlap_sites$par_cmh[[1]], as.character(cmh_clumpsites$locus_id[i]))
    if(length(multipar_site) > 0){
      overlap_sites$multipar_cmh_count[1] <- overlap_sites$multipar_cmh_count[1] + 1
      overlap_sites$multipar_cmh <- as.list(overlap_sites$multipar_cmh)
      overlap_sites$multipar_cmh[[1]] <- c(overlap_sites$multipar_cmh[[1]], as.character(cmh_clumpsites$locus_id[i]))
    }
  }
}

which(overlap_sites$par_xpehh[[1]] %in% overlap_sites$par_cmh[[1]]) #no overlap

write_rds(overlap_sites, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/Afvaper_overlaps.txt")


#----------- comparing xpehh, cmh and afvaper hits in 10 kb windows
xpehh <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_xpehh_all.txt")
colnames(xpehh) <- c("CHROM", "POS1", "POS2", "NUM_XPEHH_CRIT")
cmh <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_cmh_all.txt")
colnames(cmh) <- c("CHROM", "POS1", "POS2", "NUM_CMH_CRIT")
afvaper <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_afvaper_all.txt")
colnames(afvaper) <- c("CHROM", "POS1", "POS2", "NUM_AFVAPER_CRIT")

hits <- merge(cmh, xpehh)
hits <- merge(hits, afvaper)

hits <- hits %>%
  rowwise() %>%
  mutate("sep" = str_split(CHROM, ":"), "chr" = as.numeric(str_split(sep[[1]][1], "_")[[1]][2])) %>%
  select(-sep) %>%
  ungroup() %>%
  filter(NUM_CMH_CRIT > 0 | NUM_XPEHH_CRIT > 0 | NUM_AFVAPER_CRIT > 0) #3281

overlaps <- filter(hits, NUM_CMH_CRIT > 0 & NUM_XPEHH_CRIT > 0 & NUM_AFVAPER_CRIT > 0) #21
cmhandxpehh <- filter(hits, NUM_CMH_CRIT > 0 & NUM_XPEHH_CRIT > 0) #95
cmhandafvaper <- filter(hits, NUM_CMH_CRIT > 0 & NUM_AFVAPER_CRIT > 0) #133
afvaperandxpehh <- filter(hits, NUM_AFVAPER_CRIT > 0 & NUM_XPEHH_CRIT > 0) #57
cmhonly <- filter(hits, NUM_CMH_CRIT > 0 & NUM_XPEHH_CRIT == 0 & NUM_AFVAPER_CRIT == 0) #931
xpehhonly <- filter(hits, NUM_CMH_CRIT == 0 & NUM_XPEHH_CRIT > 0 & NUM_AFVAPER_CRIT == 0) #258
afvaperonly <- filter(hits, NUM_CMH_CRIT == 0 & NUM_XPEHH_CRIT == 0 & NUM_AFVAPER_CRIT > 0) #1849
xpehhall <- filter(hits, NUM_XPEHH_CRIT > 0) #389
afvaperall <- filter(hits, NUM_AFVAPER_CRIT > 0) #2018
cmhall <- filter(hits, NUM_CMH_CRIT > 0) #1138
cmh_or_xpehh <- filter(hits, NUM_CMH_CRIT > 0 | NUM_XPEHH_CRIT > 0)

hits_cum <- hits %>%
  filter(chr == 1) %>%
  mutate("cum_bp" = POS2)

for (s in 2:16){
  max <- max(hits_cum$cum_bp)
  hits_s <- hits %>%
    filter(chr == s) %>%
    mutate("cum_bp" = POS2 + max)
  hits_cum <- add_row(hits_cum, hits_s)
  
  hit_chr_plot <- ggplot() + geom_point(data = filter(hits_s, hits_s$NUM_CMH_CRIT > 0), aes(x = POS1, y = 1, color = "CMH")) + geom_point(data = filter(hits_s, hits_s$NUM_XPEHH_CRIT > 0), aes(x = POS1, y = 2, color = "XPEHH")) + geom_point(data = filter(hits_s, hits_s$NUM_AFVAPER_CRIT > 0), aes(x = POS1, y = 3, color = "AFVAPER")) + scale_color_manual(values = c("#382A54FF", "#3497A9FF", "#60CEACFF")) + theme_bw() + theme(axis.text = element_text(size = 15))
  ggsave(plot = hit_chr_plot, filename = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/overlaps/overlaps_allmethod_chr", s, ".png", sep = ""), device = "png", dpi = 300, height = 6, width = 10, units = "in")
}

hit_plot <- ggplot() + geom_point(data = filter(hits_cum, hits_cum$NUM_CMH_CRIT > 0), aes(x = cum_bp, y = 1, color = "CMH")) + geom_point(data = filter(hits_cum, hits_cum$NUM_XPEHH_CRIT > 0), aes(x = cum_bp, y = 2, color = "XPEHH")) + geom_point(data = filter(hits_cum, hits_cum$NUM_AFVAPER_CRIT > 0), aes(x = cum_bp, y = 3, color = "AFVAPER")) + scale_color_manual(values = c("#382A54FF", "#3497A9FF", "#60CEACFF")) + theme_bw() + theme(axis.text = element_text(size = 15))
  ggsave(plot = hit_plot, filename = paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/overlaps/overlaps_allmethod.png", sep = ""), device = "png", dpi = 300, height = 6, width = 30, units = "in")
```

#H12 with all sites (not recombination rate filtered)
```{r}
#first- rerun visualizations!!!
# for sc in {1..16}; do Rscript /Users/libbypolston/Desktop/h12/H12_viz.R /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/dcgm_hap_stats_Scaffold_${sc}.txt /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/dcgm_hap_peaks_Scaffold_${sc}.txt /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/H12/dcgm_Scaffold_${sc}_h12scan.pdf 10; done
# 
# for sc in {1..16}; do Rscript /Users/libbypolston/Desktop/h12/hapSpectrum_viz.R /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/dcgm_hap_peaks_Scaffold_${sc}.txt /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/H12/dcgm_Scaffold_${sc}_spectrum.pdf 10 440; done


#load in ag, nat and pooled run
peaks_ag <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/dcgm_hap_peaks_Scaffold_", 1,"_ag.txt", sep = ""))
colnames(peaks_ag) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123","smallest_edge_coordofpeak", "largest_edge_coordofpeak")
peaks_nat <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/dcgm_hap_peaks_Scaffold_", 1,"_nat.txt", sep = ""))
colnames(peaks_nat) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123","smallest_edge_coordofpeak", "largest_edge_coordofpeak")
peaks_all <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/dcgm_hap_peaks_Scaffold_", 1,".txt", sep = ""))
colnames(peaks_all) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123","smallest_edge_coordofpeak", "largest_edge_coordofpeak")

#have to repeat for the hapstats (not just hap peaks)
H12Scan_all <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/dcgm_hap_stats_Scaffold_", 1,".txt", sep = ""))
colnames(H12Scan_all) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123")
H12Scan_all$chr <- rep(1, length(H12Scan_all$win_center))
H12Scan_all$bpcum <- H12Scan_all$win_center
H12Scan_ag <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/dcgm_hap_stats_Scaffold_", 1,"_ag.txt", sep = ""))
colnames(H12Scan_ag) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123")
H12Scan_ag$chr <- rep(1, length(H12Scan_ag$win_center))
H12Scan_ag$bpcum <- H12Scan_ag$win_center
H12Scan_nat <- read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/dcgm_hap_stats_Scaffold_", 1,"_nat.txt", sep = ""))
colnames(H12Scan_nat) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123")
H12Scan_nat$chr <- rep(1, length(H12Scan_nat$win_center))
H12Scan_nat$bpcum <- H12Scan_nat$win_center

#get the top 10 peaks for each chrom
peaks_ag$top_peak <- c(rep(T, 1))
peaks_ag$chr <- rep(1, length(peaks_ag$win_center))
peaks_ag$bpcum <- peaks_ag$win_center
peaks_nat$top_peak <- c(rep(T, 10), rep(F, length(peaks_nat$win_center) - 10))
peaks_nat$chr <- rep(1, length(peaks_nat$win_center))
peaks_nat$bpcum <- peaks_nat$win_center
peaks_all$top_peak <- c(rep(T, 1))
peaks_all$chr <- rep(1, length(peaks_all$win_center))
peaks_all$bpcum <- peaks_all$win_center

H12Scan_all$chr <- rep(1, length(H12Scan_all$win_center))
H12Scan_all$bpcum <- H12Scan_all$win_center
H12Scan_ag$chr <- rep(1, length(H12Scan_ag$win_center))
H12Scan_ag$bpcum <- H12Scan_ag$win_center
H12Scan_nat$chr <- rep(1, length(H12Scan_nat$win_center))
H12Scan_nat$bpcum <- H12Scan_nat$win_center

for(s in 2:16){
  peaks_ag_in = read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/dcgm_hap_peaks_Scaffold_", s,"_ag.txt", sep = ""))
  colnames(peaks_ag_in) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123","smallest_edge_coordofpeak", "largest_edge_coordofpeak")

  H12Scan_ag_in = read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/dcgm_hap_stats_Scaffold_", s,"_ag.txt", sep = ""))
  colnames(H12Scan_ag_in) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123")

  H12Scan_ag_in$chr <- rep(s, length(H12Scan_ag_in$win_center))
  peaks_ag_in$top_peak <- c(rep(T, 10), rep(F, length(peaks_ag_in$win_center) - 10))
  peaks_ag_in$chr <- rep(s, length(peaks_ag_in$win_center))

  #make bpcum
  prelength <- max(H12Scan_ag$bpcum)
  peaks_ag_in$bpcum <- peaks_ag_in$win_center + prelength
  H12Scan_ag_in$bpcum <- H12Scan_ag_in$win_center + prelength

  peaks_ag <- rbind(peaks_ag, peaks_ag_in)
  H12Scan_ag <- rbind(H12Scan_ag, H12Scan_ag_in)
  
  #nat run
  peaks_nat_in = read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/dcgm_hap_peaks_Scaffold_", s,"_nat.txt", sep = ""))
  colnames(peaks_nat_in) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123","smallest_edge_coordofpeak", "largest_edge_coordofpeak")
  H12Scan_nat_in = read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/dcgm_hap_stats_Scaffold_", s,"_nat.txt", sep = ""))
  colnames(H12Scan_nat_in) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123")

  H12Scan_nat_in$chr <- rep(s, length(H12Scan_nat_in$win_center))

  peaks_nat_in$top_peak <- c(rep(T, 10), rep(F, length(peaks_nat_in$win_center) - 10))
  peaks_nat_in$chr <- rep(s, length(peaks_nat_in$win_center))

  #make bpcum
  prelength <- max(H12Scan_nat$bpcum)
  peaks_nat_in$bpcum <- peaks_nat_in$win_center + prelength
  peaks_nat <- rbind(peaks_nat, peaks_nat_in)
  H12Scan_nat_in$bpcum <- H12Scan_nat_in$win_center + prelength
  H12Scan_nat <- rbind(H12Scan_nat, H12Scan_nat_in)
  
  #all run
  peaks_all_in = read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/dcgm_hap_peaks_Scaffold_", s,".txt", sep = ""))
  colnames(peaks_all_in) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123","smallest_edge_coordofpeak", "largest_edge_coordofpeak")
  H12Scan_all_in = read.table(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/dcgm_hap_stats_Scaffold_", s,".txt", sep = ""))
  colnames(H12Scan_all_in) <- c("win_center", "win_pos_1", "win_pos_2", "num_unique_hap", "hap_freq_spectrum", "num_hap_in_each_freqbin", "H1", "H2", "H12", "H2H1", "H123")

  H12Scan_all_in$chr <- rep(s, length(H12Scan_all_in$win_center))

  peaks_all_in$top_peak <- c(rep(T, 10), rep(F, length(peaks_all_in$win_center) - 10))
  peaks_all_in$chr <- rep(s, length(peaks_all_in$win_center))

  #make bpcum
  prelength <- max(H12Scan_all$bpcum)
  peaks_all_in$bpcum <- peaks_all_in$win_center + prelength
  H12Scan_all_in$bpcum <- H12Scan_all_in$win_center + prelength

  peaks_all <- rbind(peaks_all, peaks_all_in)
  H12Scan_all <- rbind(H12Scan_all, H12Scan_all_in)
}

#now compare top 50 peaks - this is based on top H12 value
peaks_all <- peaks_all %>%
  arrange(desc(H12))

peaks_all$top50 <- c(seq(1, 50, by = 1), rep(NA, length(peaks_all$win_center)-50))
peaks_all <- peaks_all %>%
  mutate("locusID" = paste(chr, ":", win_center, sep = ""))

peaks_ag <- peaks_ag %>%
  arrange(desc(H12))

peaks_ag$top50 <- c(seq(1, 50, by = 1), rep(NA, length(peaks_ag$win_center)-50))
peaks_ag <- peaks_ag %>%
  mutate("locusID" = paste(chr, ":", win_center, sep = ""))

peaks_nat <- peaks_nat %>%
  arrange(desc(H12))

peaks_nat$top50 <- c(seq(1, 50, by = 1), rep(NA, length(peaks_nat$win_center)-50))
peaks_nat <- peaks_nat %>%
  mutate("locusID" = paste(chr, ":", win_center, sep = ""))

H12Scan_all$top50 <- c(seq(1, 50, by = 1), rep(NA, length(H12Scan_all$win_center)-50))
H12Scan_all <- H12Scan_all %>%
  mutate("locusID" = paste(chr, ":", win_center, sep = ""))

H12Scan_ag$top50 <- c(seq(1, 50, by = 1), rep(NA, length(H12Scan_ag$win_center)-50))
H12Scan_ag <- H12Scan_ag %>%
  mutate("locusID" = paste(chr, ":", win_center, sep = ""))

H12Scan_nat$top50 <- c(seq(1, 50, by = 1), rep(NA, length(H12Scan_nat$win_center)-50))
H12Scan_nat <- H12Scan_nat %>%
  mutate("locusID" = paste(chr, ":", win_center, sep = ""))

H12Scan_all <- H12Scan_all %>%
  arrange(desc(H12)) %>% 
  distinct(locusID, .keep_all = TRUE)

H12Scan_ag <- H12Scan_ag %>%
  arrange(desc(H12)) %>% 
  distinct(locusID, .keep_all = TRUE)

H12Scan_nat <- H12Scan_nat %>%
  arrange(desc(H12)) %>% 
  distinct(locusID, .keep_all = TRUE)

write.table(peaks_ag, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/Aggregated_peaks_ag", quote = F, row.names = F)
write.table(peaks_nat, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/Aggregated_peaks_nat", quote = F, row.names = F)
write.table(peaks_all, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/Aggregated_peaks_all", quote = F, row.names = F)

write.table(H12Scan_all, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/Aggregated_H12scan_all", quote = F, row.names = F)
write.table(H12Scan_ag, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/Aggregated_H12scan_ag", quote = F, row.names = F)
write.table(H12Scan_nat, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/Aggregated_H12scan_nat", quote = F, row.names = F)



#--------Run from here---------------------


peaks_ag <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/Aggregated_peaks_ag")
peaks_nat <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/Aggregated_peaks_nat")
peaks_all <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/Aggregated_peaks_all")

H12Scan_all <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/Aggregated_H12scan_all")
H12Scan_ag <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/Aggregated_H12scan_ag")
H12Scan_nat <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/H12/Aggregated_H12scan_nat")

#filter to top 50 h12 vals genome wide, then pull out the H2H1 vals in ag and nat and then plot these
toph12all_agnatvals <- H12Scan_all %>%
  filter(is.na(top50) == F) %>%
  rowwise() %>%
  mutate(H12_ag   = H12Scan_ag$H12[match(locusID, H12Scan_ag$locusID)],
         H12_nat  = H12Scan_nat$H12[match(locusID, H12Scan_nat$locusID)],
         H2H1_ag  = H12Scan_ag$H2H1[match(locusID, H12Scan_ag$locusID)],
         H2H1_nat = H12Scan_nat$H2H1[match(locusID, H12Scan_nat$locusID)])

#now look at H2/H1 (so H12 is a sweep, H2/H1 is the softness of it)
h2h1plot <- ggplot() + geom_point(data = toph12all_agnatvals, aes(x = H2H1_ag, y = H2H1_nat))  + ggtitle("Values of H2H1 in ag and nat\nfor top H12 values in pooled run") + geom_abline(intercept = 0, slope = 1, linetype = 2) + theme_bw() + theme(axis.text = element_text(size = 15)) + scale_color_manual(values = c("#60CEACFF","#395D9CFF","#382A54FF","gray80")) + theme_bw() + theme(axis.text = element_text(size = 15))
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/h2h1_agvnat_fortop50.png", plot = h2h1plot, device = "png", dpi = 300, height = 6, width = 10, units = "in")


#Maybe above plot isn’t the right thing to ask- distribution of top 50 h2h1 vals in pooled, distribution of h2h1 vals for cmh/xpehh sites
top50pool <- filter(H12Scan_all, is.na(top50) == F)

cmh_clumpsites <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clumphits_dcgm.txt")
cmh_clumpsites$snp <- paste(cmh_clumpsites$CHROM, ":", cmh_clumpsites$BP, sep ="")

xpehh_clumpsites <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/xpehh_sigsites.txt")

#get sites that the locus ID matches (the cmh site is the window center)
cmh_ag <- H12Scan_ag %>%
  filter(locusID %in% cmh_clumpsites$snp)

cmh_nat <- H12Scan_nat %>%
  filter(locusID %in% cmh_clumpsites$snp)

cmh_pool <- H12Scan_all %>%
  filter(locusID %in% cmh_clumpsites$snp)

xpehh_ag <- H12Scan_ag %>%
  filter(locusID %in% xpehh_clumpsites$locus_id)

xpehh_nat <- H12Scan_nat %>%
  filter(locusID %in% xpehh_clumpsites$locus_id)

xpehh_pool <- H12Scan_all %>%
  filter(locusID %in% xpehh_clumpsites$locus_id)

#for the remaining sites, take top h12 value from the windows that contain the snp
for(site in 1:length(cmh_clumpsites$snp)){
  if(cmh_clumpsites$snp[site] %in% cmh_ag$locusID) {
  }else{
    row_ag <- H12Scan_ag[which(H12Scan_ag$chr == cmh_clumpsites$CHROM[site] & cmh_clumpsites$BP[site] >= H12Scan_ag$win_pos_1 & cmh_clumpsites$BP[site] <= H12Scan_ag$win_pos_2),]
    row_ag <- arrange(row_ag, desc(H12))
    cmh_ag <- add_row(cmh_ag, row_ag[1,])
  }
}

for(site in 1:length(cmh_clumpsites$snp)){
  if(cmh_clumpsites$snp[site] %in% cmh_nat$locusID) {
  }else{
    row_nat <- H12Scan_nat[which(H12Scan_nat$chr == cmh_clumpsites$CHROM[site] & cmh_clumpsites$BP[site] >= H12Scan_nat$win_pos_1 & cmh_clumpsites$BP[site] <= H12Scan_nat$win_pos_2),]
    row_nat <- arrange(row_nat, desc(H12))
    cmh_nat <- add_row(cmh_nat, row_nat[1,])
  }
}

for(site in 1:length(cmh_clumpsites$snp)){
  if(cmh_clumpsites$snp[site] %in% cmh_pool$locusID) {
  }else{  
    row_pool <- H12Scan_all[which(H12Scan_all$chr == cmh_clumpsites$CHROM[site] & cmh_clumpsites$BP[site] >= H12Scan_all$win_pos_1 & cmh_clumpsites$BP[site] <= H12Scan_all$win_pos_2),]
    row_pool <- arrange(row_pool, desc(H12))
    cmh_pool <- add_row(cmh_pool, row_pool[1,])
  }
}

for(site in 1:length(xpehh_clumpsites$locus_id)){
  if(xpehh_clumpsites$locus_id[site] %in% xpehh_pool$locusID) {
  }else{
    row_pool <- H12Scan_all[which(H12Scan_all$chr == xpehh_clumpsites$CHROM[site] & xpehh_clumpsites$POS[site] >= H12Scan_all$win_pos_1 & xpehh_clumpsites$POS[site] <= H12Scan_all$win_pos_2),] 
    row_pool <- arrange(row_pool, desc(H12))
    xpehh_pool <- add_row(xpehh_pool, row_pool[1,])
  }
}

for(site in 1:length(xpehh_clumpsites$locus_id)){
  if(xpehh_clumpsites$locus_id[site] %in% xpehh_ag$locusID) {
  }else{
    row_ag <- H12Scan_ag[which(H12Scan_ag$chr == xpehh_clumpsites$CHROM[site] & xpehh_clumpsites$POS[site] >= H12Scan_ag$win_pos_1 & xpehh_clumpsites$POS[site] <= H12Scan_ag$win_pos_2),]
    row_ag <- arrange(row_ag, desc(H12))
    xpehh_ag <- add_row(xpehh_ag, row_ag[1,])
  }
}

for(site in 1:length(xpehh_clumpsites$locus_id)){
  if(xpehh_clumpsites$locus_id[site] %in% xpehh_nat$locusID) {
  }else{
    row_nat <- H12Scan_nat[which(H12Scan_nat$chr == xpehh_clumpsites$CHROM[site] & xpehh_clumpsites$POS[site] >= H12Scan_nat$win_pos_1 & xpehh_clumpsites$POS[site] <= H12Scan_nat$win_pos_2),]
    row_nat <- arrange(row_nat, desc(H12))
    xpehh_nat <- add_row(xpehh_nat, row_nat[1,])

  }
}

#get num of sites shown - all sites!
length(xpehh_ag$win_center) - sum(is.na(xpehh_ag$win_center)) #1147
length(xpehh_nat$win_center) - sum(is.na(xpehh_nat$win_center)) #1147

length(cmh_ag$win_center) - sum(is.na(cmh_ag$win_center)) #865
length(cmh_nat$win_center) - sum(is.na(cmh_nat$win_center)) #865

length(cmh_pool$win_center) - sum(is.na(cmh_pool$win_center)) #865
length(xpehh_pool$win_center) - sum(is.na(xpehh_pool$win_center)) #1147

H2h1_dist <- ggplot() + geom_density(data = top50pool, aes(x = H2H1, fill = "top50H12"), alpha = 0.5, position = "identity") + geom_density(data = cmh_ag, aes(x = H2H1, fill = "CMH_ag"), alpha = 0.5, position = "identity") + geom_density(data = cmh_nat, aes(x = H2H1, fill = "CMH_nat"), alpha = 0.5, position = "identity") + geom_density(data = xpehh_ag, aes(x = H2H1, fill = "XPEHH_ag"), alpha = 0.5, position = "identity") + geom_density(data = xpehh_nat, aes(x = H2H1, fill = "XPEHH_nat"), alpha = 0.5, position = "identity") + geom_density(data = cmh_pool, aes(x = H2H1, fill = "CMH_pool"), alpha = 0.5, position = "identity") + geom_density(data = xpehh_pool, aes(x = H2H1, fill = "XPEHH_pool"), alpha = 0.5, position = "identity") + labs(title = "H2H1 for: top 50 h12 vals in pooled set, cmh/xpehh loci in ag and nat") + theme_bw() + theme(axis.text = element_text(size = 15))

ggsave(plot = H2h1_dist, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/h2h1_dist.png", device = "png", dpi = 300, height = 6, width = 10, units = "in")

H2h1_dist_cmh <- ggplot() + geom_histogram(data = cmh_pool, aes(y = H2H1), alpha = 0.5, position = "identity") + labs(title = "H2H1 for: cmh loci in pooled run") + theme_bw() + theme(axis.text = element_text(size = 20))
ggsave(plot = H2h1_dist_cmh, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/h2h1_dist_cmh.png", device = "png", dpi = 300, height = 6, width = 10, units = "in")

H2h1_dist_xpehh <- ggplot() + geom_histogram(data = xpehh_pool, aes(y = H2H1), alpha = 0.5, position = "identity") + labs(title = "H2H1 for: xpehh loci in pooled run") + theme_bw() + theme(axis.text = element_text(size = 20))
ggsave(plot = H2h1_dist_xpehh, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/h2h1_dist_xpehh.png", device = "png", dpi = 300, height = 6, width = 10, units = "in")

h12_dist <- ggplot() + geom_histogram(data = top50pool, aes(x = H12, fill = "top50H12"), alpha = 0.5, position = "identity") + geom_histogram(data = cmh_ag, aes(x = H12, fill = "CMH_ag"), alpha = 0.5, position = "identity") + geom_histogram(data = cmh_nat, aes(x = H12, fill = "CMH_nat"), alpha = 0.5, position = "identity") + geom_histogram(data = xpehh_ag, aes(x = H12, fill = "XPEHH_ag"), alpha = 0.5, position = "identity") + geom_histogram(data = xpehh_nat, aes(x = H12, fill = "XPEHH_nat"), alpha = 0.5, position = "identity") + geom_histogram(data = cmh_pool, aes(x = H12, fill = "CMH_pool"), alpha = 0.5, position = "identity") + geom_histogram(data = xpehh_pool, aes(x = H12, fill = "XPEHH_pool"), alpha = 0.5, position = "identity") + labs(title = "H12 for: top 50 h12 vals in pooled set, cmh/xpehh loci in ag and nat") + theme_bw() + theme(axis.text = element_text(size = 15))

ggsave(plot = h12_dist, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/h12_dist.png", device = "png", dpi = 300, height = 6, width = 10, units = "in")

#plot chromosome wide h12 values for pooled, ag and nat runs
chr_colors <- ifelse(H12Scan_all$chr %% 2 == 0, "black", "gray70")

png("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/H12/Aggregated_H12scan_all.png",width=18,height=6, units="in", res=300)
plot(H12Scan_all$bpcum, H12Scan_all$H12, pch=20, col=chr_colors,
     xlab='Position', ylab='H12', main='Pooled H12 scan all scaffolds')
points(peaks_all$bpcum[is.na(peaks_all$top50) == F], peaks_all$H12[is.na(peaks_all$top50) == F], col='red', pch=20, cex=1)
dev.off()

#plotting H2h1 for top 50
png("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/H12/Aggregated_H2H1scan_all.png",width=18,height=6, units="in", res=300)
plot(H12Scan_all$bpcum, H12Scan_all$H2H1, pch=20, col=chr_colors,
     xlab='Position', ylab='H2H1', main='Pooled H2H1 scan all scaffolds', ylim = c(1,0))
points(peaks_all$bpcum[is.na(peaks_all$top50) == F], peaks_all$H2H1[is.na(peaks_all$top50) == F], col='red', pch=20, cex=1)
dev.off()

#plotting cmh or xpehh sites in red (instead of top50 h12)
chr_colors <- ifelse(H12Scan_all$chr %% 2 == 0, "#3497A9FF", "#A0DFB9FF")
png("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/H12/Aggregated_H2H1_putsitesflip.png",width=18,height=6, units="in", res=300)
plot(H12Scan_all$bpcum, H12Scan_all$H2H1, pch=20, col=chr_colors,
     xlab='Position', ylab='H2H1', main='Pooled H2H1 scan all scaffolds', ylim = c(0,1), auto.key = TRUE)
points(cmh_pool$bpcum, cmh_pool$H2H1, col=adjustcolor("black", alpha.f = 0.3), pch=20, cex=1)
points(xpehh_pool$bpcum, xpehh_pool$H2H1, col=adjustcolor("black", alpha.f = 0.3), pch=20, cex=1)
dev.off()

chr_colors <- ifelse(H12Scan_all$chr %% 2 == 0, "#3497A9FF", "#A0DFB9FF")
png("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/H12/Aggregated_H2H1_putsites.png",width=18,height=6, units="in", res=300)
plot(H12Scan_all$bpcum, H12Scan_all$H2H1, pch=20, col=chr_colors,
     xlab='Position', ylab='H2H1', main='Pooled H2H1 scan all scaffolds', ylim = c(1,0), auto.key = TRUE)
points(cmh_pool$bpcum, cmh_pool$H2H1, col=adjustcolor("black", alpha.f = 0.3), pch=20, cex=1)
points(xpehh_pool$bpcum, xpehh_pool$H2H1, col=adjustcolor("black", alpha.f = 0.3), pch=20, cex=1)
dev.off()

chr_colors <- ifelse(H12Scan_ag$chr %% 2 == 0, "black", "gray70")

png("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/H12/Aggregated_H12scan_ag.png",width=18,height=6, units="in", res=300)
plot(H12Scan_ag$bpcum, H12Scan_ag$H12, pch=20, col=chr_colors,
     xlab='Position', ylab='H12', main='Ag H12 scan all scaffolds')
points(peaks_ag$bpcum[is.na(peaks_ag$top50) == F], peaks_ag$H12[is.na(peaks_ag$top50) == F], col='red', pch=20, cex=1)
dev.off()

chr_colors <- ifelse(H12Scan_nat$chr %% 2 == 0, "black", "gray70")

png("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/H12/Aggregated_H12scan_nat.png",width=18,height=6, units="in", res=300)
plot(H12Scan_nat$bpcum, H12Scan_nat$H12, pch=20, col=chr_colors,
     xlab='Position', ylab='H12', main='Nat H12 scan all scaffolds')
points(peaks_nat$bpcum[is.na(peaks_nat$top50) == F], peaks_nat$H12[is.na(peaks_nat$top50) == F], col='red', pch=20, cex=1)
dev.off()
```
#venn diagram for cmh and xpehh overlap in 10kb windows
```{r}
#put all the sites into 10kb windows (from above)
xpehh <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_xpehh_all.txt")
colnames(xpehh) <- c("CHROM", "POS1", "POS2", "NUM_XPEHH_CRIT")
cmh <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_cmh_all.txt")
colnames(cmh) <- c("CHROM", "POS1", "POS2", "NUM_CMH_CRIT")
afvaper <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_afvaper_all.txt")
colnames(afvaper) <- c("CHROM", "POS1", "POS2", "NUM_AFVAPER_CRIT")

hits <- merge(cmh, xpehh)
hits <- merge(hits, afvaper)

hits <- hits %>%
  rowwise() %>%
  mutate("sep" = str_split(CHROM, ":"), "chr" = as.numeric(str_split(sep[[1]][1], "_")[[1]][2])) %>%
  select(-sep) %>%
  ungroup() %>%
  filter(NUM_CMH_CRIT > 0 | NUM_XPEHH_CRIT > 0 | NUM_AFVAPER_CRIT > 0)

library(sf)
library(polylabelr)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggforce)
library(eulerr)

# ------------------------------------------------------------------
# Helper 1: turn an eulerr fit into sf circle polygons
# ------------------------------------------------------------------
circles_from_fit <- function(fit) {
  geom <- as.data.frame(fit$ellipses)
  geom$set <- rownames(geom)

  circles <- lapply(seq_len(nrow(geom)), function(i) {
    sf::st_sfc(
      sf::st_buffer(sf::st_point(c(geom$h[i], geom$k[i])),
                     dist = geom$a[i], nQuadSegs = 120)
    )
  })
  names(circles) <- geom$set
  list(circles = circles, geom = geom)
}

# ------------------------------------------------------------------
# Helper 2: for every exclusive region (single, pairwise, triple...),
# compute the polygon and its "pole of inaccessibility" label point
# ------------------------------------------------------------------
region_labels <- function(circles) {
  sets <- names(circles)
  n <- length(sets)
  subsets <- unlist(lapply(seq_len(n), function(k) combn(sets, k, simplify = FALSE)),
                     recursive = FALSE)

  purrr::map_dfr(subsets, function(s) {
    inside  <- circles[s]
    outside <- circles[setdiff(sets, s)]

    region <- Reduce(sf::st_intersection, inside)
    if (length(outside) > 0) region <- Reduce(sf::st_difference, outside, region)

    if (length(region) == 0 || sf::st_is_empty(region) ||
        as.numeric(sf::st_area(region)) < 1e-6) {
      return(NULL)
    }

    # if a sliver leaves multiple polygon parts, keep only the largest
    region <- sf::st_cast(sf::st_sfc(sf::st_union(region)), "POLYGON")
    if (length(region) > 1) {
      region <- region[which.max(sf::st_area(region))]
    }

    coords <- sf::st_coordinates(region)[, c("X", "Y")]
    pos <- polylabelr::poi(coords, precision = 0.001)

    data.frame(set = paste(s, collapse = "&"), x = pos$x, y = pos$y)
  })
}

# ------------------------------------------------------------------
# Helper 3: build + save a ggplot Euler diagram for any number of sets
# ------------------------------------------------------------------
plot_euler <- function(counts, colors, title, outfile,
                        label_col = "white", set_label_col = "black") {
  fit <- eulerr::euler(counts)
  built <- circles_from_fit(fit)
  circ_geom <- built$geom
  circ_geom$set <- rownames(circ_geom)

  labels <- region_labels(built$circles)

  circles_df <- data.frame(
    x0 = circ_geom$h, y0 = circ_geom$k, r = circ_geom$a, set = circ_geom$set
  )

  p <- ggplot() +
    geom_circle(
      data = circles_df,
      aes(x0 = x0, y0 = y0, r = r, fill = set),
      alpha = 0.75, color = "white", linewidth = 1.2
    ) +
    scale_fill_manual(name = NULL, values = colors) +
    geom_text(
      data = labels,
      aes(x = x, y = y, label = ifelse(grepl("&", set),
                                        as.character(round(exp(1))), "")), # placeholder, replaced below
      size = 0
    )

  # attach the actual counts to each label row, then draw them
  counts_lookup <- setNames(as.numeric(counts), names(counts))
  labels$n <- counts_lookup[labels$set]

  p <- ggplot() +
    geom_circle(
      data = circles_df,
      aes(x0 = x0, y0 = y0, r = r, fill = set),
      alpha = 0.75, color = "white", linewidth = 1.2
    ) +
    scale_fill_manual(name = NULL, values = colors) +
    geom_text(
      data = labels,
      aes(x = x, y = y, label = n),
      size = 6, fontface = "bold", color = label_col
    ) +
    geom_text(
      data = circles_df,
      aes(x = x0, y = y0 + r + 0.4, label = set),
      size = 5.5, fontface = "bold", color = set_label_col
    ) +
    coord_fixed() +
    labs(title = title) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
      legend.position = "none"
    )

  print(p)
  ggsave(outfile, plot = p, width = 7, height = 7, dpi = 200)
  p
}

# ------------------------------------------------------------------
# 2-set diagram: XPEHH vs CMH
# ------------------------------------------------------------------
cmh_only   <- filter(hits, NUM_CMH_CRIT > 0 & NUM_XPEHH_CRIT == 0)
xpehh_only <- filter(hits, NUM_CMH_CRIT == 0 & NUM_XPEHH_CRIT > 0)
both       <- filter(hits, NUM_CMH_CRIT > 0 & NUM_XPEHH_CRIT > 0)

plot_euler(
  counts = c("XPEHH" = nrow(xpehh_only), "CMH" = nrow(cmh_only),
             "XPEHH&CMH" = nrow(both)),
  colors = c("XPEHH" = "#8AD9B1FF", "CMH" = "#395D9CFF"),
  title  = "XPEHH vs CMH hits in 10kb windows Overlap",
  outfile = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/venn_xpehh_cmh.png"
)

# ------------------------------------------------------------------
# 3-set diagram: CMH vs XPEHH vs AFVAPER
# ------------------------------------------------------------------
cmhonly           <- filter(hits, NUM_CMH_CRIT > 0 & NUM_XPEHH_CRIT == 0 & NUM_AFVAPER_CRIT == 0)
xpehhonly         <- filter(hits, NUM_CMH_CRIT == 0 & NUM_XPEHH_CRIT > 0 & NUM_AFVAPER_CRIT == 0)
afvaperonly       <- filter(hits, NUM_CMH_CRIT == 0 & NUM_XPEHH_CRIT == 0 & NUM_AFVAPER_CRIT > 0)
cmhxpehh_only     <- filter(hits, NUM_CMH_CRIT > 0 & NUM_XPEHH_CRIT > 0 & NUM_AFVAPER_CRIT == 0)
cmhafvaper_only   <- filter(hits, NUM_CMH_CRIT > 0 & NUM_AFVAPER_CRIT > 0 & NUM_XPEHH_CRIT == 0)
xpehhafvaper_only <- filter(hits, NUM_XPEHH_CRIT > 0 & NUM_AFVAPER_CRIT > 0 & NUM_CMH_CRIT == 0)
overlaps          <- filter(hits, NUM_CMH_CRIT > 0 & NUM_XPEHH_CRIT > 0 & NUM_AFVAPER_CRIT > 0)

stopifnot(
  nrow(cmhonly) + nrow(xpehhonly) + nrow(afvaperonly) +
  nrow(cmhxpehh_only) + nrow(cmhafvaper_only) + nrow(xpehhafvaper_only) +
  nrow(overlaps) == nrow(hits)
)

plot_euler(
  counts = c(
    "CMH" = nrow(cmhonly), "XPEHH" = nrow(xpehhonly), "AFVAPER" = nrow(afvaperonly),
    "CMH&XPEHH" = nrow(cmhxpehh_only), "CMH&AFVAPER" = nrow(cmhafvaper_only),
    "XPEHH&AFVAPER" = nrow(xpehhafvaper_only),
    "CMH&XPEHH&AFVAPER" = nrow(overlaps)
  ),
  colors = c("CMH" = "#395D9CFF", "XPEHH" = "#8AD9B1FF", "AFVAPER" = "#382A54FF"),
  title  = "CMH vs XPEHH vs AFVAPER hits in 10kb windows Overlap",
  outfile = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/venn_cmh_xpehh_afvaper.png"
)
```

#looking for GO enrichment
```{r}
# xpehh <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_xpehh_all.txt")
# colnames(xpehh) <- c("CHROM", "POS1", "POS2", "NUM_XPEHH_CRIT")
# xpehh <- xpehh %>% filter(NUM_XPEHH_CRIT > 0)
# 
# cmh <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_cmh_all.txt")
# colnames(cmh) <- c("CHROM", "POS1", "POS2", "NUM_CMH_CRIT")
# cmh <- cmh %>% filter(NUM_CMH_CRIT > 0)
# 
# afvaper <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/window_afvaper_all.txt")
# colnames(afvaper) <- c("CHROM", "POS1", "POS2", "NUM_AFVAPER_CRIT")
# afvaper <- afvaper %>% filter(NUM_AFVAPER_CRIT > 0)
# 
# write.table(xpehh, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/xpehhonly_GO.bed", quote = F, col.names = F, row.names = F, sep = "\t")
# write.table(cmh, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/cmhonly_GO.bed", quote = F, col.names = F, row.names = F, sep = "\t")
# write.table(afvaper, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/afvaperonly_GO.bed", quote = F, col.names = F, row.names = F, sep = "\t")

#get atub gene from gff
#bedtools intersect -a /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/gffgene.bed -b /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/cmhonly_GO.bed -wo > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/cmh_bedtools.txt
#bedtools intersect -a /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/gffgene.bed -b /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/xpehhonly_GO.bed -wo > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/xpehh_bedtools.txt
#bedtools intersect -a /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/gffgene.bed -b /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/afvaperonly_GO.bed -wo > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/afvaper_bedtools.txt

cmh_genes <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/cmh_bedtools.txt")
colnames(cmh_genes) <- c("CHROM_gff", "POS1_gff", "POS2_gff", "INFO", "CHROM_methodhit", "POS1_methodhit", "POS2_methodhit", "NUM_CMH_CRIT", "NUM_BP_OVERLAP")
xpehh_genes <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/xpehh_bedtools.txt")
colnames(xpehh_genes) <- c("CHROM_gff", "POS1_gff", "POS2_gff", "INFO", "CHROM_methodhit", "POS1_methodhit", "POS2_methodhit", "NUM_XPEHH_CRIT", "NUM_BP_OVERLAP")
afvaper_genes <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/afvaper_bedtools.txt")
colnames(afvaper_genes) <- c("CHROM_gff", "POS1_gff", "POS2_gff", "INFO", "CHROM_methodhit", "POS1_methodhit", "POS2_methodhit", "NUM_AFVAPER_CRIT", "NUM_BP_OVERLAP")

#get similar to name for each gene that overlaps
cmh_unique_genes <- cmh_genes %>%
  separate(INFO, into = c("ID", "INFO"), sep = ";") %>%
  separate(ID, into = c("INFO2", "ID"), sep = "ID=") %>%
  group_by(ID) %>%
  summarise(ID_Atub = unique(ID), MAX_NUM_HIT = max(NUM_CMH_CRIT), N_WIN_SHAREGENE = n()) %>%
  ungroup()
xpehh_unique_genes <- xpehh_genes %>%
  separate(INFO, into = c("ID", "INFO"), sep = ";") %>%
  separate(ID, into = c("INFO2", "ID"), sep = "ID=") %>%
  group_by(ID) %>%
  summarise(ID_Atub = unique(ID), MAX_NUM_HIT = max(NUM_XPEHH_CRIT), N_WIN_SHAREGENE = n()) %>%
  ungroup()
afvaper_unique_genes <- afvaper_genes %>%
  separate(INFO, into = c("ID", "INFO"), sep = ";") %>%
  separate(ID, into = c("INFO2", "ID"), sep = "ID=") %>%
  group_by(ID) %>%
  summarise(ID_Atub = unique(ID), MAX_NUM_HIT = max(NUM_AFVAPER_CRIT), N_WIN_SHAREGENE = n()) %>%
  ungroup()

#atub genes in arabidopsis
genes <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/Arabidopsis_ATgenes_tuberculatsmatches.txt")
genes_clean <- genes %>%
  separate_rows(`Atub_193_hap2.all.maker.proteins`, sep = ",\\s*") %>%
  mutate(ID_Atub_clean = str_remove(`Atub_193_hap2.all.maker.proteins`, "-RA$") %>% str_trim())

matches_cmh <- cmh_unique_genes %>%
  inner_join(genes_clean, by = c("ID_Atub" = "ID_Atub_clean"))
matches_xpehh <- xpehh_unique_genes %>%
  inner_join(genes_clean, by = c("ID_Atub" = "ID_Atub_clean"))
matches_afvaper <- afvaper_unique_genes %>%
  inner_join(genes_clean, by = c("ID_Atub" = "ID_Atub_clean"))

#get list of arabidopsis proteins for GO
matches_arab_cmh <- separate_rows(matches_cmh, Arabdopsis_protiens, sep = ",\\s*")
cmh_toanalyze <- unique(matches_arab_cmh$Arabdopsis_protiens)
matches_arab_afvaper <- separate_rows(matches_afvaper, Arabdopsis_protiens, sep = ",\\s*")
afvaper_toanalyze <- unique(matches_arab_afvaper$Arabdopsis_protiens)
matches_arab_xpehh <- separate_rows(matches_xpehh, Arabdopsis_protiens, sep = ",\\s*")
xpehh_toanalyze <- unique(matches_arab_xpehh$Arabdopsis_protiens)

hit_inany_method_toanalyze <- c(cmh_toanalyze, afvaper_toanalyze)
hit_inany_method_toanalyze <- c(hit_inany_method_toanalyze, xpehh_toanalyze)
hit_inany_method_toanalyze <- unique(hit_inany_method_toanalyze)

#get list of all possible arabidopsis proteins for go
ref_arab <- separate_rows(genes, Arabdopsis_protiens, sep = ",\\s*")
ref <- ref_arab$Arabdopsis_protiens

write.table(cmh_toanalyze, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/cmh_atub_forGO.txt", col.names = F, row.names = F, quote = F)
write.table(xpehh_toanalyze, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/xpehh_atub_forGO.txt", col.names = F, row.names = F, quote = F)
write.table(afvaper_toanalyze, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/afvaper_atub_forGO.txt", col.names = F, row.names = F, quote = F)
write.table(hit_inany_method_toanalyze, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/hit_inany_method_toanalyze_forGO.txt", col.names = F, row.names = F, quote = F)
write.table(ref, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ref_forGO.txt", col.names = F, row.names = F, quote = F)

#repeating for the set of 2-3 methods overlap in that window
hits <- merge(cmh, xpehh)
hits <- merge(hits, afvaper)

hits <- hits %>%
  rowwise() %>%
  mutate("sep" = str_split(CHROM, ":"), "chr" = as.numeric(str_split(sep[[1]][1], "_")[[1]][2])) %>%
  select(-sep) %>%
  ungroup() %>%
  filter(NUM_CMH_CRIT > 0 | NUM_XPEHH_CRIT > 0 | NUM_AFVAPER_CRIT > 0)

cmhonly           <- filter(hits, NUM_CMH_CRIT > 0 & NUM_XPEHH_CRIT == 0 & NUM_AFVAPER_CRIT == 0)
xpehhonly         <- filter(hits, NUM_CMH_CRIT == 0 & NUM_XPEHH_CRIT > 0 & NUM_AFVAPER_CRIT == 0)
afvaperonly       <- filter(hits, NUM_CMH_CRIT == 0 & NUM_XPEHH_CRIT == 0 & NUM_AFVAPER_CRIT > 0)
cmhxpehh_only     <- filter(hits, NUM_CMH_CRIT > 0 & NUM_XPEHH_CRIT > 0 & NUM_AFVAPER_CRIT == 0)
cmhafvaper_only   <- filter(hits, NUM_CMH_CRIT > 0 & NUM_AFVAPER_CRIT > 0 & NUM_XPEHH_CRIT == 0)
xpehhafvaper_only <- filter(hits, NUM_XPEHH_CRIT > 0 & NUM_AFVAPER_CRIT > 0 & NUM_CMH_CRIT == 0)
overlaps          <- filter(hits, NUM_CMH_CRIT > 0 & NUM_XPEHH_CRIT > 0 & NUM_AFVAPER_CRIT > 0)


#making list of regions with 2 or 3 overlaps
shared_hits <- merge(cmhafvaper_only, cmhxpehh_only, all = T)
shared_hits <- merge(shared_hits, xpehhafvaper_only, all = T)
shared_hits <- merge(shared_hits, overlaps, all = T)

#confirm all unique
shared_hits <- shared_hits %>%
  mutate("snp_id" = paste(CHROM, "_", POS1, ":", POS2, sep = ""))
length(unique(shared_hits$snp_id))
shared_hits <- shared_hits %>%
  select(CHROM, POS1, POS2, NUM_CMH_CRIT, NUM_XPEHH_CRIT, NUM_AFVAPER_CRIT)

write.table(shared_hits, file = "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/dcgm_methodoverlap.bed", sep = "\t", row.names = F, col.names = F, quote = F)

#bedtools intersect -a /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/gffgene.bed -b /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/dcgm_methodoverlap.bed -wo > /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/dcgm_methodoverlap_bedtools.txt

shared_hits_genes <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/dcgm_methodoverlap_bedtools.txt")
colnames(shared_hits_genes) <- c("CHROM_gff", "POS1_gff", "POS2_gff", "INFO", "CHROM_methodhit", "POS1_methodhit", "POS2_methodhit", "NUM_CMH_CRIT", "NUM_XPEHH_CRIT", "NUM_AFVAPER_CRIT", "NUM_BP_OVERLAP")

#get similar to name for each gene that overlaps
shared_unique_genes <- shared_hits_genes %>%
  separate(INFO, into = c("ID", "INFO"), sep = ";") %>%
  separate(ID, into = c("INFO2", "ID"), sep = "ID=") %>%
  group_by(ID) %>%
  summarise(ID_Atub = unique(ID), MAX_CMH = max(NUM_CMH_CRIT), MAX_XPEHH = max(NUM_XPEHH_CRIT), MAX_AFVAPER = max(NUM_AFVAPER_CRIT), N_WIN_SHAREGENE = n()) %>%
  ungroup()

#atub genes in arabidopsis
genes <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/Arabidopsis_ATgenes_tuberculatsmatches.txt")
genes_clean <- genes %>%
  separate_rows(`Atub_193_hap2.all.maker.proteins`, sep = ",\\s*") %>%
  mutate(ID_Atub_clean = str_remove(`Atub_193_hap2.all.maker.proteins`, "-RA$") %>% str_trim())

matches <- shared_unique_genes %>%
  inner_join(genes_clean, by = c("ID_Atub" = "ID_Atub_clean"))

#get list of arabidopsis proteins for GO
matches_arab <- separate_rows(matches, Arabdopsis_protiens, sep = ",\\s*")
matches_toanalyze <- unique(matches_arab$Arabdopsis_protiens)

#get list of all possible arabidopsis proteins for go
ref_arab <- separate_rows(genes, Arabdopsis_protiens, sep = ",\\s*")
ref <- ref_arab$Arabdopsis_protiens

#this table gives the unique gene names for all of our windows that are shared by two or more methods (243 originally), then gives the number of windows that had that gene name (in case 2 windows overlapped w single gene), and the value of the hits for each method that were the greatest for that gene
write.table(matches, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/genes_sharedacrossmethods.txt", row.names = F, quote = F)

write.table(matches_toanalyze, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/matches_atub_forGO.txt", col.names = F, row.names = F, quote = F)
write.table(ref, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ref_forGO.txt", col.names = F, row.names = F, quote = F)
```

#Quick info

```{r}
xpehh_clumpsites <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/xpehh_sigsites.txt")
cmh_clumpsites <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/clumphits_dcgm.txt")
dcgm_.05 <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/P.05FDR_dcgm")

colors <- c("#0B0405FF", "#382A54FF", "#395D9CFF", "#3497A9FF", "#60CEACFF", "#DEF5E5FF")
#add to every plot: 
theme_bw() + theme(axis.text = element_text(size = 15))
```

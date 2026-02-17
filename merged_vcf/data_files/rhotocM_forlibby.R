# Recombination Rate to cM conversion

# Load required libraries
library(data.table)   # Fast data reading and manipulation
library(tidyverse)    # Data manipulation and visualization toolkit


# Read Y chromosome recombination data (males only)
recomb_LDhat <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/allchroms_193reference_ldmapresults.txt2",sep="\t")  # tabs or multiple spaces

# Standardize column names
names(recomb_LDhat) <- c("Scaf", "Loci", "Mean_rho", "Median_rho", "L95", "U95")
recomb_LDhat$pos<-recomb_LDhat$Loci * 1000 #to use bp positions - easier
unique(recomb_LDhat$Scaf)
recomb_LDhat$chr_num <- as.numeric(gsub("Scaffold_([0-9]+).*", "\\1", recomb_LDhat$Scaf))

# Quick histogram of variation in mean rho
hist(recomb_LDhat$Mean_rho)

# Create windowed analysis by chromosome ()
winrecomb <-  recomb_LDhat %>% 
  mutate(win = floor(pos/100000)) %>%
  group_by(chr_num,win) %>%
  dplyr::summarise(mean_winrho = mean(Mean_rho))

# Plot chromosomal recombination landscape with region boundaries
recomb_landscape <- ggplot(data = winrecomb, aes(win/10, mean_winrho)) +
  geom_point() +
  geom_smooth(span = .1) +
  labs(y = "Mean 4NeR\n in 100kb Windows",
       x = "Genomic Position (Mb)") +
  theme_bw() +
  facet_grid(~chr_num, scales = "free_x",space = "free_x")
  
recomb_landscape

# Calculate cumulative recombination rate across each chromosome
recomb_LDhat <- recomb_LDhat %>%
  group_by(chr_num) %>%
  mutate(cumulativerho = cumsum(Mean_rho)) %>%
  ungroup()

winrecombcum <-  recomb_LDhat %>% 
  mutate(win = floor(pos/100000)) %>%
  group_by(chr_num,win) %>%
  dplyr::summarise(mean_winrhocum = mean(cumulativerho))

recombcum_landscape <- ggplot(data = winrecombcum, aes(win/10, mean_winrhocum)) +
  geom_point() +
  geom_smooth(span = .1) +
  labs(y = "Cumulative 4NeR\n in 100kb Windows",
       x = "Genomic Position (Mb)") +
  theme_bw() +
  facet_grid(~chr_num, scales = "free_x",space = "free_x")

recombcum_landscape #note that some chromosomes seem to have longer maps despite being shorter! (i.e. 2 vs 1)


#########

#Ne from Kreiner et al. PNAS 2018
ne=5411804

#convert to cM
recomb_LDhat$r<-recomb_LDhat$cumulativerho/(4*ne)
recomb_LDhat$cM<-recomb_LDhat$r*100


map_forshapeit<-recomb_LDhat %>% select(pos, Scaf, cM)
write.table(map_forshapeit,"/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/LDbased_cMestimates_bychrom.txt",quote=F,col.names = F, row.names = F)

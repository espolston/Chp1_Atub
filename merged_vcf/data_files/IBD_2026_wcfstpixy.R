#!/usr/bin/env Rscript
# adapted from IBD_2026.R from Julia
# Tests IBD (isolation by distance) and IBE (isolation by environment)
# was using pairwise wc Fst, geographic distance, and environment (Ag vs Nat)
#updated 6/2/26 for wc fst on each chrom

#pop1_1 is pair 1, ag pop1_2 is pair 1, nat

library(tidyverse)
library(geosphere)   # haversine distances
library(vegan)       # mantel test
library(ggplot2)
library(data.table)

# ── 0. Load pixy data to convert from wc fst to wc 50 mb fst ─────────────
fst_wc <- tibble("pop1" = character(), "pop2" = character(), "chromosome" = character(), "window_pos_1" = numeric(), "window_pos_2" = numeric(), "avg_wc_fst" = numeric(), "no_snps" = numeric(),	"wc_fst_a" = numeric(),	"wc_fst_b" = numeric(),	"wc_fst_c" = numeric())

for(p in 1:16){
  pix_file <- fread(paste("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pixy/pixy_fst_50Mb_Scaffold_", p, ".txt", sep = ""), sep="\t", header=T)
  fst_wc <- fst_wc %>%
    add_row("pop1" = pix_file$pop1, "pop2" = pix_file$pop2, "chromosome" = pix_file$chromosome, "window_pos_1" = pix_file$window_pos_1, "window_pos_2" = pix_file$window_pos_2, "avg_wc_fst" = pix_file$avg_wc_fst, "no_snps" = pix_file$no_snps,	"wc_fst_a" = pix_file$wc_fst_a,	"wc_fst_b" = pix_file$wc_fst_b,	"wc_fst_c" = pix_file$wc_fst_c)
}

#average over all chromosomes 
fst_wc_avg <- fst_wc %>%
  mutate("contrast" = paste(pop1, "_", pop2, sep = "")) %>%
  group_by(contrast) %>%
  summarise("wcfst_allscaf" = sum(wc_fst_a) / (sum(wc_fst_a) + sum(wc_fst_b) + sum(wc_fst_c)), "pop1" = first(pop1), "pop2" = first(pop2)) %>%
  ungroup()

#I took out this coords because i think the coordinates are wrong and then numbers -> ag or nat I think are flipped???
#coords <- read_tsv("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/pop_coords.tsv")  # Pair_Env, Lat, Long, Environment
#coords <- coords %>%
#  rowwise() %>%
#  mutate("Pair" = strsplit(Pair_Env, "_")[[1]][1], "Pair_num" = strsplit(Pair, "pop")[[1]][2], "Env" = case_when(
#    Environment == 1 ~ "AG",
#    Environment == 2 ~ "NAT"
#  ), "Pair_Env_final" = paste("pair_", Pair_num, "_", Env, sep = "")) %>%
#  ungroup()
#I've update POP# -> pop#, Pair_Env -> Pair_Env_final, and now HUDSON_FST -> wcfst_allscaf

#using coords from https://docs.google.com/spreadsheets/d/1XN9sKfZL9NgFWdRhXZCtK0cyJifjuPaeXkcjWryU5K4/edit?gid=0#gid=0
coords <- read_tsv("/Users/libbypolston/Desktop/Popcoordinates.tsv")
coords <- coords %>%
  rowwise() %>%
  mutate("Pair_Env_final" = paste("pair_", Pair, "_", Env, sep = "")) %>%
  ungroup()

# ── 1. Load data ─────────────────────────────────────────────────────────────
#setwd("/Users/libbypolston/Desktop/")
#fst <- read_tsv("fst_pairwise.fst.summary", comment = "#",
#                col_names = c("POP1","POP2","HUDSON_FST"))

#coords <- read_tsv("pop_coords.tsv")  # Pair_Env, Lat, Long, Environment

# ── 2. Attach coordinates and environment to each population pair ─────────────
fst <- fst_wc_avg %>%
  left_join(coords, by = c("pop1" = "Pair_Env_final")) %>%
  #rename(Lat1 = Lat, Long1 = Long, Env1 = Environment) %>%
  rename(Lat1 = Lat, Long1 = Long, Env1 = Env, Pair1 = Pair) %>%
  left_join(coords, by = c("pop2" = "Pair_Env_final")) %>%
  rename(Lat2 = Lat, Long2 = Long, Env2 = Env, Pair2 = Pair)
  #rename(Lat2 = Lat, Long2 = Long, Env2 = Environment)

# ── 3. Compute geographic distance (haversine, km) ───────────────────────────
fst <- fst %>%
  rowwise() %>%
  mutate(
    Geo_km   = distHaversine(c(Long1, Lat1), c(Long2, Lat2)) / 1000,
    log_Geo  = log(Geo_km),
    Fst_lin  = pmax(wcfst_allscaf, 0) / (1 - pmax(wcfst_allscaf, 0)),
    Env_diff = as.integer(Env1 != Env2),
    comp_type = case_when(
      Env1 == Env2 & Env1 == "AG" ~ "Ag vs Ag",
      Env1 == Env2 & Env1 == "NAT" ~ "Nat vs Nat",
      TRUE                      ~ "Ag vs Nat"
    ),
    Pair1 = str_extract(pop1, "pair_\\d+"),
    Pair2 = str_extract(pop2, "pair_\\d+")
  ) %>%
  ungroup() %>%
  filter(Geo_km > 0, is.finite(Fst_lin))

cat(sprintf("Total pairwise comparisons: %d\n", nrow(fst)))
print(table(fst$comp_type))

# ── 4. Split into IBD and IBE subsets ────────────────────────────────────────
fst_same   <- fst %>% filter(Env_diff == 0)
fst_within <- fst %>% filter(Pair1 == Pair2, Env_diff == 1)

cat(sprintf("\nIBD same-env pairs: %d\n", nrow(fst_same)))
cat(sprintf("IBE within-pair comparisons: %d\n", nrow(fst_within)))

#------- 4a. adding in fst by distance---------------
#plotting like in /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/FST_genomewideavg_bypair.png
#filter to just the within pair contrasts
fst_within %>%
  ggplot(aes(x = wcfst_allscaf, fill = T)) +
  geom_density(alpha = 0.25, linewidth = 0.6) + scale_fill_manual(values = c("gray")) +
  labs(x = "Mean FST") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 15))

fst_difhab <- fst %>% filter(Env_diff == 1) #both pair 1 ag- pair 1 nat and pair 1 ag - pair 2 nat ...

fst_agag <- fst_same %>% filter(Env1 == "AG")
fst_natnat <- fst_same %>% filter(Env1 == "NAT")

#get between all slope, within env, between env slopes
anyenv <- lm(wcfst_allscaf ~ Geo_km, fst)
bwenv <- lm(wcfst_allscaf ~ Geo_km, fst_difhab)
withinenv <- lm(wcfst_allscaf ~ Geo_km, fst_same)
withinenv_ag <- lm(wcfst_allscaf ~ Geo_km, fst_agag)
withinenv_nat <- lm(wcfst_allscaf ~ Geo_km, fst_natnat)

#compare fit with geom_smooth(method="loess") geom_smooth(method="gam") and lm
anyenv_loess <- loess(wcfst_allscaf ~ Geo_km, fst)
bwenv_loess <- loess(wcfst_allscaf ~ Geo_km, fst_difhab)
withinenv_loess <- loess(wcfst_allscaf ~ Geo_km, fst_same)
withinenv_ag_loess <- loess(wcfst_allscaf ~ Geo_km, fst_agag)
withinenv_nat_loess <- loess(wcfst_allscaf ~ Geo_km, fst_natnat)

library(mgcv)
anyenv_gam <- gam(wcfst_allscaf ~ Geo_km, data = fst)
bwenv_gam <- gam(wcfst_allscaf ~ Geo_km, data = fst_difhab)
withinenv_gam <- gam(wcfst_allscaf ~ Geo_km, data = fst_same)
withinenv_ag_gam <- gam(wcfst_allscaf ~ Geo_km, data = fst_agag)
withinenv_nat_gam <- gam(wcfst_allscaf ~ Geo_km, data = fst_natnat)

loessRsq <- function(model_obj, dataset){
  y_obs <- dataset$wcfst_allscaf
  y_mean <- mean(y_obs)
  
  tss <- sum((y_obs - y_mean)^2)
  
  rss_loess <- sum(residuals(model_obj)^2)
  r2_loess <- 1 - (rss_loess / tss)
  
  return(r2_loess)
  
  #to check correct method
  #r2_lm    <- 1 - (rss_lm / tss)
  #rss_lm <- sum(residuals(anyenv)^2) 
}

summary(anyenv) #rsq = 0.3725
loessRsq(anyenv_loess, fst) #.386
summary(anyenv_gam) #rsq = .371

summary(bwenv) #rsq = 0.3752
loessRsq(bwenv_loess, fst_difhab) #.389
summary(bwenv_gam) #rsq = 0.373

summary(withinenv) #rsq = 0.372
loessRsq(withinenv_loess, fst_same) #.383
summary(withinenv_gam) #rsq = 0.37

summary(withinenv_ag) #rsq = 0.3484
loessRsq(withinenv_ag_loess, fst_agag) #.358
summary(withinenv_ag_gam) #rsq = 0.344

summary(withinenv_nat) #rsq = 0.4096
loessRsq(withinenv_nat_loess, fst_natnat) #.436
summary(withinenv_nat_gam) #rsq = 0.405

#best model is loess

#now plot fst vs distance
#dist_withinpair <- ggplot(fst_within, aes(x = Geo_km, y = wcfst_allscaf)) + geom_point() + labs(x = "Mean FST", y = "Distance (km)") 
#ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ibd_output/Fst_withinpair_distance.png", dist_withinpair, width = 9, height = 8)

#because above looked so bad go to between all ag-nat, ag-ag, nat-nat pairs
#dist <- ggplot(fst, aes(x = Geo_km, y = wcfst_allscaf)) + geom_point() + labs(y = "Mean FST", x = "Distance (km)", title = "Fst ~ Distance for all ag-nat, ag-ag, nat-nat pairs") + geom_abline(aes(intercept = coef(anyenv)[1], slope = coef(anyenv)[2], color = "All pairs")) + geom_abline(aes(intercept = coef(bwenv)[1], slope = coef(bwenv)[2], color = "Ag-nat pairs")) + geom_abline(aes(intercept = coef(withinenv)[1], slope = coef(withinenv)[2], color = "Same env pairs"))

#ggplot(fst, aes(x = Geo_km, y = wcfst_allscaf, color=comp_type)) + geom_point() + geom_smooth(method="gam") + labs(y = "Mean FST", x = "Distance (km)", title = "Fst ~ Distance for all ag-nat, ag-ag, nat-nat pairs") 
  
#ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ibd_output/Fst_allpair_ibd_linear.png", dist, width = 9, height = 8)

dist <- ggplot(fst, aes(x = Geo_km, y = wcfst_allscaf, color = comp_type)) + geom_point() + labs(y = "Mean FST", x = "Distance (km)", title = "Fst ~ Distance for all ag-nat, ag-ag, nat-nat pairs") + geom_smooth(method="loess")
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ibd_output/Fst_allpair_ibd_linear.png", dist, width = 9, height = 8)

#compare fst of ag-ag and ag-nat of similar distances to each other 
dist_bwpair <- ggplot(fst_same, aes(x = Geo_km, y = wcfst_allscaf, color = comp_type)) + geom_point() + labs(y = "Mean FST", x = "Distance (km)", title = "Fst ~ Distance for ag-ag or nat-nat \nDashed line is any within env comparison") + geom_smooth(method="loess") + geom_smooth(data = fst_same, aes(x = Geo_km, y = wcfst_allscaf),method = "loess", inherit.aes = FALSE, color = "black", linetype = "dashed")
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ibd_output/Fst_samenv_ibd.png", dist_bwpair, width = 9, height = 8)

#ag-nat of any pair
dist_bwenv <- ggplot(fst_difhab, aes(x = Geo_km, y = wcfst_allscaf, color = comp_type)) + geom_point() + labs(y = "Mean FST", x = "Distance (km)", title = "Fst ~ Distance for ag-nat of any pair") + geom_smooth(method="loess")
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ibd_output/Fst_difenv_ibd.png", dist_bwenv, width = 9, height = 8)

# ── 5. Mantel tests ───────────────────────────────────────────────────────────
pops_same <- sort(unique(c(fst_same$pop1, fst_same$pop2)))

make_mat <- function(df, val_col, pops) {
  n <- length(pops)
  m <- matrix(0, n, n, dimnames = list(pops, pops))
  for (i in seq_len(nrow(df))) {
    p1 <- df$pop1[i]; p2 <- df$pop2[i]
    if (p1 %in% pops & p2 %in% pops) {
      v <- df[[val_col]][i]
      m[p1, p2] <- v; m[p2, p1] <- v
    }
  }
  m
}

fst_mat <- make_mat(fst_same, "Fst_lin", pops_same)
geo_mat <- make_mat(fst_same, "log_Geo", pops_same)
env_mat <- make_mat(fst_same, "Env_diff", pops_same)

cat("\n── Partial Mantel: IBD controlling for environment ──\n")
pops_all <- sort(unique(c(fst$pop1, fst$pop2)))
fst_mat_all <- make_mat(fst, "Fst_lin", pops_all)
geo_mat_all <- make_mat(fst, "log_Geo", pops_all)
env_mat_all <- make_mat(fst, "Env_diff", pops_all)

cat("\n── Partial Mantel: IBD controlling for environment (all pairs) ──\n")
mantel_ibd_partial <- mantel.partial(as.dist(fst_mat_all), as.dist(geo_mat_all),
                                     as.dist(env_mat_all),
                                     method = "pearson", permutations = 999)
print(mantel_ibd_partial)

cat("\n── Partial Mantel: IBE controlling for distance (all pairs) ──\n")
mantel_ibe_partial <- mantel.partial(as.dist(fst_mat_all), as.dist(env_mat_all),
                                     as.dist(geo_mat_all),
                                     method = "pearson", permutations = 999)
print(mantel_ibe_partial)

# ── MMRR implementation (Wang 2013) ──────────────────────────────────────────
unfold <- function(X) {
  x <- vector()
  for (i in 2:nrow(X)) x <- c(x, X[i, 1:(i-1)])
  x
}


unfold <- function(X) {
  x <- vector()
  for (i in 2:nrow(X)) x <- c(x, X[i, 1:(i-1)])
  x
}
MMRR <- function(Y, X, nperm = 999) {
  nrowsY <- nrow(Y)
  y      <- unfold(Y)
  if (is.null(names(X))) names(X) <- paste0("X", seq_along(X))
  Xmats  <- sapply(X, unfold)
  fit    <- lm(y ~ Xmats)
  coeffs <- fit$coefficients
  summ   <- summary(fit)
  r2     <- summ$r.squared
  tstat  <- summ$coefficients[, "t value"]
  Fstat  <- summ$fstatistic[1]
  tprob  <- rep(1, length(tstat))
  Fprob  <- 1
  
  for (i in seq_len(nperm)) {
    rand  <- sample(seq_len(nrowsY))
    yperm <- unfold(Y[rand, rand])
    fit   <- lm(yperm ~ Xmats)
    summ  <- summary(fit)
    Fprob <- Fprob + as.numeric(summ$fstatistic[1] >= Fstat)
    tprob <- tprob + as.numeric(abs(summ$coefficients[, "t value"]) >= abs(tstat))
  }
  
  tp <- tprob / (nperm + 1)
  Fp <- Fprob / (nperm + 1)
  names(coeffs) <- c("Intercept", names(X))
  names(tstat)  <- paste0(c("Intercept", names(X)), "(t)")
  names(tp)     <- paste0(c("Intercept", names(X)), "(p)")
  
  list(r.squared    = r2,
       coefficients = coeffs,
       tstatistic   = tstat,
       tpvalue      = tp,
       Fstatistic   = Fstat,
       Fpvalue      = Fp)
}

# Build interaction matrix: element-wise product of geo and env matrices
# (analogous to log_Geo * Env_diff interaction term in lm)
pops_all    <- sort(unique(c(fst$pop1, fst$pop2)))
fst_mat_all <- make_mat(fst, "Fst_lin",  pops_all)
geo_mat_all <- make_mat(fst, "log_Geo",  pops_all)
env_mat_all <- make_mat(fst, "Env_diff", pops_all)

# Interaction matrix = elementwise product
geo_x_env <- geo_mat_all * env_mat_all

# Run MMRR
result <- MMRR(fst_mat_all,
               list(geography = geo_mat_all,
                    environment = env_mat_all,
                    geo_x_env = geo_x_env),
               nperm = 999)
print(result)
# ── MMRR: IBD + IBE + interaction ────────────────────────────────────────────
pops_all    <- sort(unique(c(fst$pop1, fst$pop2)))
fst_mat_all <- make_mat(fst, "Fst_lin",  pops_all)
geo_mat_all <- make_mat(fst, "log_Geo",  pops_all)
env_mat_all <- make_mat(fst, "Env_diff", pops_all)
geo_x_env   <- geo_mat_all * env_mat_all  # interaction matrix

mmrr_res <- MMRR(fst_mat_all,
                 list(geography   = geo_mat_all,
                      environment = env_mat_all,
                      geo_x_env   = geo_x_env),
                 nperm = 999)

cat("\n── MMRR: IBD + IBE + interaction ──\n")
cat(sprintf("R² = %.3f\n", mmrr_res$r.squared))
print(data.frame(
  coefficient = mmrr_res$coefficients,
  t_stat      = mmrr_res$tstatistic,
  p_value     = mmrr_res$tpvalue
))
cat(sprintf("F = %.3f, p = %.3f\n", mmrr_res$Fstatistic, mmrr_res$Fpvalue))

#test seperately
# Ag vs Ag pairs only
pops_ag  <- sort(unique(c(filter(fst_same, comp_type == "Ag vs Ag")$pop1,
                          filter(fst_same, comp_type == "Ag vs Ag")$pop2)))
fst_ag   <- make_mat(filter(fst_same, comp_type == "Ag vs Ag"),  "Fst_lin", pops_ag)
geo_ag   <- make_mat(filter(fst_same, comp_type == "Ag vs Ag"),  "log_Geo", pops_ag)

# Nat vs Nat pairs only
pops_nat <- sort(unique(c(filter(fst_same, comp_type == "Nat vs Nat")$pop1,
                          filter(fst_same, comp_type == "Nat vs Nat")$pop2)))
fst_nat  <- make_mat(filter(fst_same, comp_type == "Nat vs Nat"), "Fst_lin", pops_nat)
geo_nat  <- make_mat(filter(fst_same, comp_type == "Nat vs Nat"), "log_Geo", pops_nat)

cat("\n── MMRR: IBD in Ag populations only ──\n")
res_ag  <- MMRR(fst_ag,  list(geography = geo_ag),  nperm = 999)
cat(sprintf("Ag slope  = %.4f, p = %.3f, R² = %.3f\n",
            res_ag$coefficients["geography"],
            res_ag$tpvalue["geography(p)"],
            res_ag$r.squared))

cat("\n── MMRR: IBD in Nat populations only ──\n")
res_nat <- MMRR(fst_nat, list(geography = geo_nat), nperm = 999)
cat(sprintf("Nat slope = %.4f, p = %.3f, R² = %.3f\n",
            res_nat$coefficients["geography"],
            res_nat$tpvalue["geography(p)"],
            res_nat$r.squared))





# ── 7. Plots ─────────────────────────────────────────────────────────────────
cols <- c("Ag vs Ag" = "#60CEACFF", "Nat vs Nat" = "#382A54FF")

# ── Plot 1: IBD faceted by comparison type ───────────────────────────────────
p_ibd <- ggplot(fst_same, aes(Geo_km, Fst_lin, colour = comp_type)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", formula = y ~ log(x), se = TRUE, linewidth = 0.8) +
  scale_x_continuous(
    trans  = "log",
    breaks = c(1, 10, 50, 100, 500, 1000),
    labels = c("1", "10", "50", "100", "500", "1000"),
    limits = c(min(fst_same$Geo_km) * 0.9, max(fst_same$Geo_km) * 1.1)
  ) +
  coord_cartesian(ylim = c(0, NA)) +   # floor at 0, don't clip data
  scale_colour_manual(values = cols, name = NULL) +
  facet_wrap(~comp_type) +
  labs(x = "Geographic distance (km)",
       y = expression(F[ST] / (1 - F[ST])),
       title = "Isolation by Distance — same environment pairs") +
  theme_classic(base_size = 12) +
  theme(legend.position  = "none",
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold"), 
        axis.text = element_text(size = 15))
p_ibd
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ibd_output/IBD_plot.pdf", p_ibd, width = 8, height = 4)
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ibd_output/IBD_plot.png", p_ibd, width = 8, height = 4, dpi = 300)

# ── Plot 2: IBD both env types overlaid ──────────────────────────────────────
p_ibd_overlay <- ggplot(fst_same, aes(Geo_km, Fst_lin, colour = comp_type)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", formula = y ~ log(x), se = TRUE, linewidth = 0.9) +
  scale_x_continuous(
    trans  = "log",
    breaks = c(1, 10, 50, 100, 500, 1000),
    labels = c("1", "10", "50", "100", "500", "1000"),
    limits = c(2, 1100)   # start at ~min data, no extrapolation left
  ) +
  coord_cartesian(ylim = c(0, NA)) +  # floor y at 0
  scale_colour_manual(values = cols, name = NULL) +
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("Mantel r = %.3f, p = %.3f",
                           mantel_ibd_partial$statistic, mantel_ibd_partial$signif),
           hjust = 1.05, vjust = -0.5, size = 3.5, colour = "grey20") +
  labs(x = "Geographic distance (km)",
       y = expression(F[ST] / (1 - F[ST])),
       title = "Isolation by Distance") +
  theme_bw() +
  theme(legend.position   = c(0.15, 0.85),
        legend.background = element_blank(), 
        axis.text = element_text(size = 15)) +
  scale_x_continuous(
  trans  = "log",
  breaks = c(1, 10, 50, 100, 500, 1000),
  labels = c("1", "10", "50", "100", "500", "1000"),
  limits = c(2, 1100)   # start at ~min data, no extrapolation left
) +
coord_cartesian(ylim = c(0, NA))  # floor y at 0
p_ibd_overlay
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ibd_output/IBD_overlay.pdf", p_ibd_overlay, width = 6, height = 5)
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ibd_output/IBD_overlay.png", p_ibd_overlay, width = 6, height = 5, dpi = 300)

# ── Plot 3: IBE — within-pair Ag vs Nat Fst ──────────────────────────────────
p_ibe <- ggplot(fst_within, aes(x = 1, y = Fst_lin)) +
  geom_boxplot(width = 0.3, fill = "#4DAF4A", alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.05, alpha = 0.7, size = 2.5, colour = "#4DAF4A") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 1, y = max(fst_within$Fst_lin),
           label = sprintf("Mean Fst = %.4f\nn = %d sites",
                           mean(fst_within$wcfst_allscaf), nrow(fst_within)),
           vjust = -0.5, size = 3.5, colour = "grey20") +
  labs(x = NULL,
       y = expression(F[ST] / (1 - F[ST])),
       title = "Isolation by Environment",
       subtitle = "Within-pair Ag vs Nat (geography controlled)") +
  theme_classic(base_size = 13) +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text = element_text(size = 15))
p_ibe
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ibd_output/IBE_plot.pdf", p_ibe, width = 4, height = 5)
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ibd_output/IBE_plot.png", p_ibe, width = 4, height = 5, dpi = 300)

# ── Plot 4: Fst heatmap ordered west to east by longitude ────────────────────
pop_order <- coords %>%
  arrange(Long) %>%
  pull(Pair_Env_final)

# Make symmetric for heatmap
fst_heatmap <- fst %>%
  dplyr::select(pop1, pop2, wcfst_allscaf) %>%
  bind_rows(tibble(pop1 = fst$pop2, pop2 = fst$pop1, wcfst_allscaf = fst$wcfst_allscaf)) %>%
  # add diagonal (self = 0)
  bind_rows(tibble(pop1 = pop_order, pop2 = pop_order, wcfst_allscaf = 0)) %>%
  mutate(pop1 = factor(pop1, levels = pop_order),
         pop2 = factor(pop2, levels = pop_order))

#a couple that end up w neg fst and thats not showing up so set them to 0
negrow <- which(fst_heatmap$wcfst_allscaf < 0)
fst_heatmap$wcfst_allscaf[negrow] <- 0

p_heatmap <- ggplot(fst_heatmap, aes(pop1, pop2, fill = wcfst_allscaf)) +
  geom_tile() +
  scale_fill_viridis_c(option = "mako", name = expression(F[ST]),
                       limits = c(0, NA)) +
  labs(x = NULL, y = NULL, title = "Pairwise wc Fst",
       subtitle = "Populations ordered west to east") +
  theme_classic(base_size = 10) +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y     = element_text(size = 10),
        legend.position = "right",
        panel.grid      = element_blank())
fst_heatmap
p_heatmap
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ibd_output/Fst_heatmap.pdf", p_heatmap, width = 9, height = 8)
ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ibd_output/Fst_heatmap.png", p_heatmap, width = 9, height = 8, dpi = 300)

# ── 8. Save results ───────────────────────────────────────────────────────────
write_tsv(fst,        "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ibd_output/ibd_ibe_pairwise.tsv")
write_tsv(fst_same,   "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ibd_output/ibd_same_env.tsv")
write_tsv(fst_within, "/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ibd_output/ibe_within_pair.tsv")

cat("\nDone. Outputs in /Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Results/ibd_output/:\n")
cat("  ibd_ibe_pairwise.tsv\n")
cat("  IBD_plot.pdf/png       — faceted by comparison type\n")
cat("  IBD_overlay.pdf/png    — both env types overlaid\n")
cat("  IBE_plot.pdf/png       — within-pair Ag vs Nat\n")
cat("  Fst_heatmap.pdf/png    — all pairwise Fst, ordered west to east\n")


#adding distance to fst plot

# Amaranthus tuberculatus — variable maps over coordinates (34 pops, 17 Ag/Nat pairs)
# tidyverse (readr/dplyr/ggplot2) + patchwork. Faint lines connect each Ag/Nat pair.
# Run: Rscript make_maps.R   (expects ../data/master_all_variables.csv)

#julia made and ran this but i have the output (keeping so i know how it was made)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(patchwork); library(maps)
})

#i think this might be amaranthus_predictors_for_libby.csv
df <- read_csv("../data/master_all_variables.csv", show_col_types = FALSE) |>
  mutate(Site = if_else(env == "Ag", "Agricultural", "Natural"))

# fixed jitter (mostly vertical) so pair-connecting lines meet the plotted points
set.seed(42)
df <- df |>
  mutate(xj = decimalLongitude + runif(n(), -0.04, 0.04),
         yj = decimalLatitude  + runif(n(), -0.16, 0.16))

states <- map_data("state")
xlim <- c(-96.2, -83.0); ylim <- c(37.8, 42.2)
shape_vals <- c("Agricultural" = 21, "Natural" = 24)

panel <- function(fill_col, title, palette, direction = 1, guide_title = "") {
  ggplot() +
    geom_polygon(data = states, aes(long, lat, group = group),
                 fill = NA, colour = "grey60", linewidth = 0.3) +
    geom_line(data = df, aes(xj, yj, group = pair),
              colour = "grey45", linewidth = 0.3, alpha = 0.7) +           # pair links
    geom_point(data = df, aes(xj, yj, fill = .data[[fill_col]], shape = Site),
               colour = "black", stroke = 0.35, size = 5, alpha = 0.9) +
    scale_shape_manual(values = shape_vals, name = NULL) +
    scale_fill_distiller(palette = palette, direction = direction,
                         name = guide_title, guide = guide_colourbar(order = 1)) +
    guides(shape = guide_legend(order = 2, override.aes = list(fill = "grey50"))) +
    coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          panel.grid = element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey70"),
          legend.position = "right",
          legend.key.height = unit(1, "null"),   # colourbar = panel height
          legend.key.width  = unit(0.4, "cm"),
          legend.title = element_text(size = 9))
}

p1 <- panel("proportion_agricultural_1000m_2019", "Proportion agricultural (1 km, NLCD 2019)", "YlGn", 1, "prop.")
p2 <- panel("glyphosate_kg_per_ha_crop", "Glyphosate intensity (kg / ha cropland, 2018)", "PuRd", 1, "kg/ha")
p3 <- panel("bio1_annual_mean_temp_C", "Annual mean temperature (BIO1, °C)", "RdYlBu", -1, "°C")
p4 <- panel("bio12_annual_precip_mm", "Annual precipitation (BIO12, mm)", "YlGnBu", 1, "mm")

combined <- (p1 | p2) / (p3 | p4) +
  plot_annotation(title = "Amaranthus tuberculatus — 34 populations (17 Ag/Nat pairs): land use, herbicide & climate (2018)",
    theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5)))

ggsave("../maps/variable_maps.pdf", combined, width = 13, height = 6.5, device = cairo_pdf)
message("wrote ../maps/variable_maps.pdf")

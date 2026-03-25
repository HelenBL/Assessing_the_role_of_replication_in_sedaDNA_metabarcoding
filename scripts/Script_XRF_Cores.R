### XRF Data Analysis

library(readxl)
library(tidyverse)
library(readr)      


archivo <- "STO_CCO_CLR_XRF.xlsx"  
metales <- c("As", "Cr", "Pb", "Cu", "Zn", "Fe")

sitio_labels <- c(
  CCO0 = "Cicero (CCO)",
  CLRm = "Colindres (CLR)",
  STOm = "Santoña (STO)"
)

clean_lod <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  
  x[x %in% c("<LOD", "< LOD", "LOD", "ld", "LD")] <- NA
  
  out <- parse_number(x)
  return(out)
}

sitios <- excel_sheets(archivo)

df_list <- lapply(sitios, function(s) {
  read_excel(archivo, sheet = s) %>%
    select(any_of(c("Depth","Year" , metales))) %>%  
    mutate(
      Depth = as.numeric(Depth),
      across(any_of(metales), clean_lod),
      Sitio = s
    )
})

df_all <- bind_rows(df_list)

df_all <- df_all %>%
  mutate(
    Depth = as.numeric(Depth),
    across(all_of(metales), clean_lod),
    Sitio = recode(Sitio, !!!sitio_labels)
  )

metales_presentes <- intersect(metales, names(df_all))

# FIGURE Raw concentrations


df_long <- df_all %>%
  pivot_longer(
    cols = all_of(metales_presentes),
    names_to = "Metal",
    values_to = "Concentration"
  )

p_conc <- ggplot(df_long, aes(x = Concentration, y = Depth, color = Metal)) +
  geom_path(linewidth = 0.9) +
  scale_y_reverse() +
  facet_wrap(~ Sitio, scales = "free") +
  scale_color_viridis_d(option = "D", end = 0.9) +  
  theme_bw(base_size = 13) +
  labs(
    title = "Heavy metals (XRF) – Raw concentrations",
    x = "Concentration (ppm)",
    y = "Depth (cm)",
    color = "Metal"
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave("Fig1_concentrations.png", p_conc,
       width = 12, height = 11, dpi = 300)


# FIGURE Scaled metals (Z-score)


df_z <- df_long %>%
  group_by(Sitio, Metal) %>%
  mutate(Z = as.numeric(scale(Concentration))) %>%
  ungroup()

p_z <- ggplot(df_z, aes(x = Z, y = Depth, color = Metal)) +
  geom_path(linewidth = 0.9) +
  scale_y_reverse() +
  facet_wrap(~ Sitio, scales = "free") +
  scale_color_viridis_d(option = "D", end = 0.9) +
  theme_bw(base_size = 13) +
  labs(
    title = "Heavy metals (XRF) – Scaled values (Z-score)",
    x = "Z-score",
    y = "Depth (cm)",
    color = "Metal"
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave("Fig2_Zscore_STO.png", p_z,
       width = 12, height = 11, dpi = 300)


# FIGURE Integrated Pollution Index (IPI)

df_icp <- df_all %>%
  group_by(Sitio) %>%
  mutate(
    across(all_of(metales_presentes),
           ~ as.numeric(scale(.)),
           .names = "z_{.col}"),
    IPI = rowMeans(across(starts_with("z_")), na.rm = TRUE)
  ) %>%
  ungroup()

p_ipi <- ggplot(df_icp, aes(x = IPI, y = Depth)) +
  geom_path(linewidth = 1, color = "black") +
  scale_y_reverse() +
  facet_wrap(~ Sitio, scales = "free") +
  theme_bw(base_size = 13) +
  labs(
    title = "Integrated heavy metal pollution index",
    x = "Pollution index (mean Z-score of As, Cr, Pb, Cu, Zn)",
    y = "Depth (cm)"
  ) +
  theme(
    panel.grid.minor = element_blank()
  )

ggsave("Fig3_IPI_STO.png", p_ipi,
       width = 12, height = 11, dpi = 300)


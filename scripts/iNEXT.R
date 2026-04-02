### script done with HPC 

library(iNEXT)
library(ggplot2)
library(openxlsx)
library(svglite)

# ==================== CONFIG ====================

outdir <- "iNEXT_results/COI"
dir.create(outdir, showWarnings = FALSE)

input_file <- "Clean_Data_COI_Reps_AbTot.xlsx"
metadata_file <- "metadata_repliques_f.xlsx"

colors_site <- c(
  "CCO" = "#8B3A3A",
  "CLR" = "cyan4",
  "STO" = "#CDAD00"
)

nboot_value <- 150

# ==================== LEER TABLA 18S ====================

tab18 <- read.xlsx(input_file)

rownames(tab18) <- tab18[[1]]
tab18[[1]] <- NULL

tab18[] <- lapply(tab18, as.numeric)

tab18t <- t(tab18)


inext_input18 <- apply(tab18t, 1, function(x) x[x > 0])

inext_input18 <- inext_input18[sapply(inext_input18, sum) > 0]

# ==================== EJECUTAR iNEXT ====================

out18 <- iNEXT(
  inext_input18,
  q = 0,
  datatype = "abundance",
  nboot = nboot_value
)

# ==================== GUARDAR OBJETO COMPLETO ====================

saveRDS(out18, file = file.path(outdir, "COI_iNEXT_output.rds"))

# ==================== GUARDAR TABLAS ====================

write.csv(
  out18$iNextEst$size_based,
  file = file.path(outdir, "COI_iNextEst_size_based.csv"),
  row.names = FALSE
)

write.csv(
  out18$iNextEst$coverage_based,
  file = file.path(outdir, "COI_iNextEst_coverage_based.csv"),
  row.names = FALSE
)

write.csv(
  out18$AsyEst,
  file = file.path(outdir, "COI_AsyEst.csv"),
  row.names = FALSE
)

write.csv(
  out18$DataInfo,
  file = file.path(outdir, "COI_DataInfo.csv"),
  row.names = FALSE
)

# ==================== GRAFICOS BASE DE iNEXT ====================

p_18s_size <- ggiNEXT(out18, type = 1) + theme_bw()
p_18s_cov  <- ggiNEXT(out18, type = 3) + theme_bw()

ggsave(
  file.path(outdir, "COI_curve_sample_size.png"),
  plot = p_18s_size,
  width = 10, height = 7, dpi = 300
)

ggsave(
  file.path(outdir, "COI_curve_sample_size.svg"),
  plot = p_18s_size,
  width = 10, height = 7,
  device= svglite::svglite
)

ggsave(
  file.path(outdir, "COI_curve_coverage_default.png"),
  plot = p_18s_cov,
  width = 10, height = 7, dpi = 300
)

ggsave(
  file.path(outdir, "COI_curve_coverage_default.svg"),
  plot = p_18s_cov,
  width = 10, height = 7, device= svglite::svglite
)

# ==================== LEER METADATA ====================

metadata <- read.xlsx(metadata_file)

if ("X1" %in% colnames(metadata)) {
  colnames(metadata)[colnames(metadata) == "X1"] <- "Assemblage"
}

# si ya existe una columna Sample, la renombramos a Assemblage
if ("Sample" %in% colnames(metadata)) {
  colnames(metadata)[colnames(metadata) == "Sample"] <- "Assemblage"
}

if (!"Assemblage" %in% colnames(metadata)) {
  stop("No encuentro la columna de nombres de muestra en metadata. Debe llamarse 'Assemblage', 'Sample' o 'X1'.")
}

if (!"Site" %in% colnames(metadata)) {
  stop("No encuentro la columna 'Site' en metadata_repliques_f.xlsx")
}

metadata$Assemblage <- as.character(metadata$Assemblage)
metadata$Site <- as.character(metadata$Site)

# ==================== UNIR iNEXT + METADATA ====================

cov_df <- out18$iNextEst$coverage_based
cov_df$Assemblage <- as.character(cov_df$Assemblage)

cov_meta <- merge(cov_df, metadata, by = "Assemblage", all.x = TRUE)

write.csv(
  cov_meta,
  file = file.path(outdir, "COI_iNextEst_coverage_based_with_metadata.csv"),
  row.names = FALSE
)

# ==================== GRAFICO COBERTURA COLOREADO POR SITE ====================

p_18s_cov_site <- ggplot(
  cov_meta,
  aes(x = SC, y = qD, group = Assemblage, color = Site)
) +
  geom_line(alpha = 0.5, linewidth = 0.7) +
  scale_color_manual(values = colors_site, drop = FALSE) +
  labs(
    x = "Sample coverage",
    y = "Species richness (q = 0)",
    color = "Site",
    title = "iNEXT coverage-based rarefaction/extrapolation by Site"
  ) +
  theme_bw()

ggsave(
  file.path(outdir, "COI_curve_coverage_by_Site.png"),
  plot = p_18s_cov_site,
  width = 10, height = 7, dpi = 300
)

ggsave(
  file.path(outdir, "COI_curve_coverage_by_Site.svg"),
  plot = p_18s_cov_site,
  width = 10, height = 7, device= svglite::svglite
)

# ==================== GRAFICO RESUMIDO POR SITE ====================

site_summary <- aggregate(
  qD ~ Site + SC,
  data = cov_meta,
  FUN = mean,
  na.rm = TRUE
)

write.csv(
  site_summary,
  file = file.path(outdir, "COI_coverage_summary_by_Site.csv"),
  row.names = FALSE
)

p_18s_cov_site_mean <- ggplot(
  site_summary,
  aes(x = SC, y = qD, color = Site)
) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = colors_site, drop = FALSE) +
  labs(
    x = "Sample coverage",
    y = "Mean species richness (q = 0)",
    color = "Site",
    title = "Mean coverage-based rarefaction/extrapolation by Site"
  ) +
  theme_bw()

ggsave(
  file.path(outdir, "COI_curve_coverage_mean_by_Site.png"),
  plot = p_18s_cov_site_mean,
  width = 10, height = 7, dpi = 300
)
ggsave(
  file.path(outdir, "COI_curve_coverage_mean_by_Site.svg"),
  plot = p_18s_cov_site_mean,
  width = 10, height = 7, device= svglite::svglite
)


cat("Análisis iNEXT COI completado correctamente.\n")
cat("Resultados guardados en:", outdir, "\n")




##### 18S

library(iNEXT)
library(ggplot2)
library(openxlsx)
library(svglite)

# ==================== CONFIG ====================

outdir <- "iNEXT_results"
dir.create(outdir, showWarnings = FALSE)

input_file <- "Clean_Data_18S_Reps_AbTot.xlsx"
metadata_file <- "metadata_repliques_f.xlsx"

colors_site <- c(
  "CCO" = "#8B3A3A",
  "CLR" = "cyan4",
  "STO" = "#CDAD00"
)

nboot_value <- 150

# ==================== LEER TABLA 18S ====================

tab18 <- read.xlsx(input_file)

# asumimos que la primera columna contiene el nombre del ASV
rownames(tab18) <- tab18[[1]]
tab18[[1]] <- NULL

# asegurar que todo es numérico
tab18[] <- lapply(tab18, as.numeric)

# transponer: filas = muestras, columnas = ASVs
tab18t <- t(tab18)

# iNEXT con datatype = "abundance" espera lista
# cada elemento = una muestra, con abundancias > 0
inext_input18 <- apply(tab18t, 1, function(x) x[x > 0])

# quitar muestras vacías por seguridad
inext_input18 <- inext_input18[sapply(inext_input18, sum) > 0]

# ==================== EJECUTAR iNEXT ====================

out18 <- iNEXT(
  inext_input18,
  q = 0,
  datatype = "abundance",
  nboot = nboot_value
)

# ==================== GUARDAR OBJETO COMPLETO ====================

saveRDS(out18, file = file.path(outdir, "18S_iNEXT_output.rds"))

# ==================== GUARDAR TABLAS ====================

# iNextEst es una lista con size_based y coverage_based
write.csv(
  out18$iNextEst$size_based,
  file = file.path(outdir, "18S_iNextEst_size_based.csv"),
  row.names = FALSE
)

write.csv(
  out18$iNextEst$coverage_based,
  file = file.path(outdir, "18S_iNextEst_coverage_based.csv"),
  row.names = FALSE
)

write.csv(
  out18$AsyEst,
  file = file.path(outdir, "18S_AsyEst.csv"),
  row.names = FALSE
)

write.csv(
  out18$DataInfo,
  file = file.path(outdir, "18S_DataInfo.csv"),
  row.names = FALSE
)

# ==================== GRAFICOS BASE DE iNEXT ====================

p_18s_size <- ggiNEXT(out18, type = 1) + theme_bw()
p_18s_cov  <- ggiNEXT(out18, type = 3) + theme_bw()

ggsave(
  file.path(outdir, "18S_curve_sample_size.png"),
  plot = p_18s_size,
  width = 10, height = 7, dpi = 300
)

ggsave(
  file.path(outdir, "18S_curve_sample_size.svg"),
  plot = p_18s_size,
  width = 10, height = 7,
  device= svglite::svglite
)

ggsave(
  file.path(outdir, "18S_curve_coverage_default.png"),
  plot = p_18s_cov,
  width = 10, height = 7, dpi = 300
)

ggsave(
  file.path(outdir, "18S_curve_coverage_default.svg"),
  plot = p_18s_cov,
  width = 10, height = 7, device= svglite::svglite
)

# ==================== LEER METADATA ====================

metadata <- read.xlsx(metadata_file)

if ("X1" %in% colnames(metadata)) {
  colnames(metadata)[colnames(metadata) == "X1"] <- "Assemblage"
}

if ("Sample" %in% colnames(metadata)) {
  colnames(metadata)[colnames(metadata) == "Sample"] <- "Assemblage"
}

if (!"Assemblage" %in% colnames(metadata)) {
  stop("No encuentro la columna de nombres de muestra en metadata. Debe llamarse 'Assemblage', 'Sample' o 'X1'.")
}

if (!"Site" %in% colnames(metadata)) {
  stop("No encuentro la columna 'Site' en metadata_repliques_f.xlsx")
}

metadata$Assemblage <- as.character(metadata$Assemblage)
metadata$Site <- as.character(metadata$Site)

# ==================== UNIR iNEXT + METADATA ====================

cov_df <- out18$iNextEst$coverage_based
cov_df$Assemblage <- as.character(cov_df$Assemblage)

cov_meta <- merge(cov_df, metadata, by = "Assemblage", all.x = TRUE)

write.csv(
  cov_meta,
  file = file.path(outdir, "18S_iNextEst_coverage_based_with_metadata.csv"),
  row.names = FALSE
)

# ==================== GRAFICO COBERTURA COLOREADO POR SITE ====================

p_18s_cov_site <- ggplot(
  cov_meta,
  aes(x = SC, y = qD, group = Assemblage, color = Site)
) +
  geom_line(alpha = 0.5, linewidth = 0.7) +
  scale_color_manual(values = colors_site, drop = FALSE) +
  labs(
    x = "Sample coverage",
    y = "Species richness (q = 0)",
    color = "Site",
    title = "iNEXT coverage-based rarefaction/extrapolation by Site"
  ) +
  theme_bw()

ggsave(
  file.path(outdir, "18S_curve_coverage_by_Site.png"),
  plot = p_18s_cov_site,
  width = 10, height = 7, dpi = 300
)

ggsave(
  file.path(outdir, "18S_curve_coverage_by_Site.svg"),
  plot = p_18s_cov_site,
  width = 10, height = 7, device= svglite::svglite
)

# ==================== GRAFICO RESUMIDO POR SITE ====================

site_summary <- aggregate(
  qD ~ Site + SC,
  data = cov_meta,
  FUN = mean,
  na.rm = TRUE
)

write.csv(
  site_summary,
  file = file.path(outdir, "18S_coverage_summary_by_Site.csv"),
  row.names = FALSE
)

p_18s_cov_site_mean <- ggplot(
  site_summary,
  aes(x = SC, y = qD, color = Site)
) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = colors_site, drop = FALSE) +
  labs(
    x = "Sample coverage",
    y = "Mean species richness (q = 0)",
    color = "Site",
    title = "Mean coverage-based rarefaction/extrapolation by Site"
  ) +
  theme_bw()

ggsave(
  file.path(outdir, "18S_curve_coverage_mean_by_Site.png"),
  plot = p_18s_cov_site_mean,
  width = 10, height = 7, dpi = 300
)
ggsave(
  file.path(outdir, "18S_curve_coverage_mean_by_Site.svg"),
  plot = p_18s_cov_site_mean,
  width = 10, height = 7, device= svglite::svglite
)


cat("Análisis iNEXT 18S completado correctamente.\n")
cat("Resultados guardados en:", outdir, "\n")
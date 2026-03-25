
#############
#PLOT KINGDOM
library(tidyverse)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(openxlsx)
library(wesanderson)
library(patchwork)
library(ggh4x)

### COOOI
metadata<- read.xlsx("metadata_f.xlsx")
arel_mitjCOI_kp <- read.xlsx("Clean_Data_COI_AbRel.xlsx")

rownames(arel_mitjCOI_kp) <- arel_mitjCOI_kp$X1
arel_mitjCOI_kp <- arel_mitjCOI_kp[,c(2:21)]

colSums(arel_mitjCOI_kp[, c("CCO.C1.45", "CCO.C1.5", "CCO.C2.45", "CCO.C2.5", "STO.C2.50", "STO.C3.50", "STO.C3.5" )])  #comprovar si sumen 1:)


arel_mitjCOI_kp_long <- pivot_longer(arel_mitjCOI_kp,
                         cols = matches("^CCO|^CLR|^STO"),  # Exemple si vols totes les mostres que comencin així
                         names_to = "Sample",
                         values_to = "ReadsCOI")



abdCOI_meta <- left_join(arel_mitjCOI_kp_long, metadata, by = "Sample")
abdCOI_meta$Core_Rep <- as.factor(abdCOI_meta$Core_Rep)
abdCOI_meta$Depth <- as.factor(abdCOI_meta$Depth)
abdCOI_meta$Site <- as.factor(abdCOI_meta$Site)



### 18S
arel_mitj18S_kp <- read.xlsx("Clean_Data_18S_AbRel.xlsx")

rownames(arel_mitj18S_kp) <- arel_mitj18S_kp$X1
arel_mitj18S_kp <- arel_mitj18S_kp[,c(2:21)]

colSums(arel_mitj18S_kp[, c("CCO.C1.45", "CCO.C1.5", "CCO.C2.45", "CCO.C2.5", "STO.C2.50", "STO.C3.50", "STO.C3.5" )])  #comprovar si sumen 1:)


arel_mitj18S_kp_long <- pivot_longer(arel_mitj18S_kp,
                                     cols = matches("^CCO|^CLR|^STO"),  # Exemple si vols totes les mostres que comencin així
                                     names_to = "Sample",
                                     values_to = "Reads18S")



abd18S_meta <- left_join(arel_mitj18S_kp_long, metadata, by = "Sample")
abd18S_meta$Core_Rep <- as.factor(abd18S_meta$Core_Rep)
abd18S_meta$Depth <- as.factor(abd18S_meta$Depth)
abd18S_meta$Site <- as.factor(abd18S_meta$Site)


#PLOT KINGDOM BY SITE AND MARKER

# STO COI
plotSTO_COI<- abdCOI_meta %>%
  filter(Site == "STO") %>%  
  ggplot(aes(x = Core_Rep, y = ReadsCOI, fill = kingdom)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Depth) +  
  scale_fill_manual(values = c("darkolivegreen3","cadetblue", "darkgoldenrod1", "lightblue","darkgoldenrod3", "seashell2" ),
                    na.value = "grey95") +
  labs(
    title = "COI",  
    x = "",
    y = "Relative Abundance",
    fill = "Kingdom"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    panel.spacing.x = unit(1.5, "lines")
  )

#STO 18S
plotSTO_18S<-abd18S_meta %>%
  filter(Site == "STO") %>%  
  ggplot(aes(x = Core_Rep, y = Reads18S, fill = kingdom)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Depth) +  
  scale_fill_manual(values = c("darkolivegreen3","cadetblue", "darkgoldenrod1", "lightblue","darkgoldenrod3", "seashell2" ),
                    na.value = "grey95") +
  labs(
    title = "18S",  
    x = "Core",
    y = "Relative Abundance",
    fill = "Kingdom"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    panel.spacing.x = unit(1.5, "lines")
  )

plotSTO_COI / plotSTO_18S

combined_STO <- (plotSTO_COI / plotSTO_18S) + 
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Relative Abundance by Kingdom - Locality STO (COI & 18S)",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      legend.position = "right",  # OR "bottom", your choice
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 12)
    )
  )

combined_STO


#CCO COI
plotCCO_COI<- abdCOI_meta %>%
  filter(Site == "CCO") %>%  
  ggplot(aes(x = Core_Rep, y = ReadsCOI, fill = kingdom)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Depth) +  
  scale_fill_manual(values = c("darkolivegreen3","cadetblue", "darkgoldenrod1", "lightblue","darkgoldenrod3", "seashell2" ),
                    na.value = "grey95") +
  labs(
    title = "COI", 
    x = "",
    y = "Relative Abundance",
    fill = "Kingdom"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    panel.spacing.x = unit(1.5, "lines")
  )

#CCO 18S
plotCCO_18S<-abd18S_meta %>%
  filter(Site == "CCO") %>%  
  ggplot(aes(x = Core_Rep, y = Reads18S, fill = kingdom)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Depth) +  
  scale_fill_manual(values = c("darkolivegreen3","cadetblue", "darkgoldenrod1", "lightblue","darkgoldenrod3", "seashell2" ),
                    na.value = "grey95") +
  labs(
    title = "18S",  
    x = "Core",
    y = "Relative Abundance",
    fill = "Kingdom"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    panel.spacing.x = unit(1.5, "lines")
  )

combined_CCO <- (plotCCO_COI / plotCCO_18S) + 
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Relative Abundance by Kingdom - Locality CCO (COI & 18S)",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      legend.position = "right",  # OR "bottom", your choice
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 12)
    )
  )

combined_CCO



#CLR COI
plotCLR_COI<- abdCOI_meta %>%
  filter(Site == "CLR") %>%  
  ggplot(aes(x = Core_Rep, y = ReadsCOI, fill = kingdom)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Depth) +  
  scale_fill_manual(values = c("darkolivegreen3","cadetblue", "darkgoldenrod1", "lightblue","darkgoldenrod3", "seashell2" ),
                    na.value = "grey95") +
  labs(
    title = "COI",  
    x = "",
    y = "Relative Abundance",
    fill = "Kingdom"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    panel.spacing.x = unit(1.5, "lines")
  )

#CLR 18S
plotCLR_18S<-abd18S_meta %>%
  filter(Site == "CLR") %>%  
  ggplot(aes(x = Core_Rep, y = Reads18S, fill = kingdom)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Depth) +  
  scale_fill_manual(values = c("darkolivegreen3","cadetblue", "darkgoldenrod1", "lightblue","darkgoldenrod3", "seashell2" ),
                    na.value = "grey95") +
  labs(
    title = "18S",  
    x = "Core",
    y = "Relative Abundance",
    fill = "Kingdom"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    panel.spacing.x = unit(1.5, "lines")
  )


combined_CLR <- (plotCLR_COI / plotCLR_18S) + 
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Relative Abundance by Kingdom - Locality CLR (COI & 18S)",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      legend.position = "right",  # OR "bottom", your choice
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 12)
    )
  )

combined_CLR






##### All Together

df_all <- bind_rows(
  abdCOI_meta %>%
    transmute(
      Site      = Site,
      Core_Rep  = as.factor(Core_Rep),
      Year     = Year,
      kingdom   = kingdom,
      Marker    = "COI",
      RelAbund  = ReadsCOI
    ),
  abd18S_meta %>%
    transmute(
      Site      = Site,
      Core_Rep  = as.factor(Core_Rep),
      Year     = Year,
      kingdom   = kingdom,
      Marker    = "18S",
      RelAbund  = Reads18S
    )
)


df_all <- df_all %>%
  mutate(
    Site = factor(Site, levels = c("STO", "CCO", "CLR")),
    FacetCol = paste0(Marker, " – ", Year)
  )


kingdom_order <- df_all %>%
  group_by(kingdom) %>%
  summarise(Total = sum(RelAbund, na.rm = TRUE)) %>%
  arrange(desc(Total)) %>%
  pull(kingdom)

df_all <- df_all %>%
  mutate(kingdom = factor(kingdom, levels = kingdom_order))


pal <- c(
  "seashell2", "cadetblue", "darkgoldenrod1",
  "lightblue", "darkgoldenrod3", "darkolivegreen3", "grey80"
)
if (length(levels(df_all$kingdom)) > length(pal)) {
  pal <- scales::hue_pal()(length(levels(df_all$kingdom)))
}


df_all <- df_all |>
  group_by(Site, FacetCol) |>
  filter(any(RelAbund > 0)) |>
  ungroup()

p_all <- ggplot(
  df_all,
  aes(x = Core_Rep, y = RelAbund, fill = kingdom)
) +
  geom_col(width = 0.85) +
  facet_grid(
    rows   = vars(Site),
    cols   = vars(FacetCol),
    switch = "y"
  ) +
  scale_fill_manual(values = pal, na.value = "grey95", name = "Kingdom") +
    labs(
    x = "Core",
    y = "Relative abundance"
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text.x = element_text(size = 15, face = "bold"),
    strip.text.y = element_text(size = 15, face = "bold"),
    panel.spacing.x = unit(1.2, "lines"),
    panel.spacing.y = unit(0.8, "lines")
  )

p_all



df_all_clean <- df_all %>% 
  group_by(Site, FacetCol, Core_Rep) %>% 
  filter(sum(RelAbund, na.rm = TRUE) > 0) %>% 
  ungroup()

df_all_clean <- df_all_clean %>%
  mutate(
    FacetCol_clean = FacetCol,
    FacetCol_clean = gsub("[–—]", "-", FacetCol_clean),  # normaliza guiones
    FacetCol_clean = gsub("\\s+", " ", FacetCol_clean),  # normaliza espacios
    FacetCol_clean = trimws(FacetCol_clean)
  )

levels_facet <- sort(unique(df_all_clean$FacetCol_clean))
df_all_clean$FacetCol_clean <- factor(df_all_clean$FacetCol_clean,
                                            levels = levels_facet)
year_cols <- c(
  "2014" = "#DCA0A0",
  "1939" = "#8B3A3A",
  "2017" = "#66C2C2",
  "1952" = "cyan4",
  "1974" = "#E6D46A",
  "1524" = "#CDAD00"
)

colors_depth_site2 <- setNames(
  year_cols[sub(".*-\\s*", "", levels_facet)],  # extrae año tras el "- "
  levels_facet
)

colors_depth_site2[is.na(colors_depth_site2)]



plot_site <- function(site_id, add_x_label = FALSE, add_y_label = FALSE) {
  
  df_site <- df_all_clean %>%
    filter(Site == site_id) %>%
    mutate(FacetCol_clean = droplevels(FacetCol_clean))  # <- CLAVE
  
  levs_site <- levels(df_site$FacetCol_clean)
  
  fill_list_site <- as.list(colors_depth_site2[levs_site])
  
  if (any(is.na(unlist(fill_list_site)))) {
    print(levs_site)
    print(colors_depth_site2[levs_site])
    stop("Hay algún FacetCol_clean sin color asignado (NA).")
  }
  
  ggplot(df_site, aes(x = Core_Rep, y = RelAbund, fill = kingdom)) +
    geom_col(width = 0.85) +
    facet_grid2(
      rows   = vars(Site),
      cols   = vars(FacetCol_clean),
      switch = "y",
      scales = "free_x",
      space  = "free_x",
      strip  = strip_themed(
        background_x = elem_list_rect(
          fill   = fill_list_site,  
          colour = "black"
        ),
        text_x = elem_list_text(
          colour = "black",
          face   = "bold",
          size   = 15
        )
      )
    ) +
    scale_fill_manual(values = pal, na.value = "grey95", name = "Phylum") +
    labs(
      x = if (add_x_label) "Biological Replicate" else NULL,
      y = if (add_y_label) "Relative read abundance" else NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      panel.grid.major.x = element_blank(),
      legend.position = "right",
      axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 15),
      strip.text.y = element_text(size = 15, face = "bold"),
      panel.spacing.x = unit(1.2, "lines"),
      panel.spacing.y = unit(0.8, "lines"),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20)
    )
}

p_sto <- plot_site("STO")
p_cco <- plot_site("CCO", add_y_label = TRUE)
p_clr <- plot_site("CLR", add_x_label = TRUE)

p_all <- p_sto / p_cco / p_clr +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18),   
    legend.text  = element_text(size = 16),   
    legend.key.size = unit(1.2, "lines")      
  )
p_all



##### Separated genes --> FIGURE OK

plot_site_marker <- function(site_id, marker_id,
                             add_x_label = FALSE, add_y_label = FALSE) {
  
  df_site <- df_all_clean %>%
    filter(Site == site_id, Marker == marker_id) %>%
    mutate(
      Year = factor(Year, levels = names(year_cols)),
      Year = droplevels(Year)
    )
  
  levs_year <- levels(df_site$Year)
  fill_list_year <- as.list(year_cols[levs_year])
  
  if (any(is.na(unlist(fill_list_year)))) {
    print(levs_year)
    print(year_cols[levs_year])
    stop("Hay algún Year sin color asignado (NA).")
  }
  
  ggplot(df_site, aes(x = Core_Rep, y = RelAbund, fill = kingdom)) +
    geom_col(width = 0.85) +
    facet_grid2(
      rows   = vars(Site),
      cols   = vars(Year),
      switch = "y",
      scales = "free_x",
      space  = "free_x",
      strip  = strip_themed(
        background_x = elem_list_rect(
          fill   = fill_list_year,   # <- colores por Year
          colour = "black"
        ),
        text_x = elem_list_text(
          colour = "black",
          face   = "bold",
          size   = 15
        )
      )
    ) +
    scale_fill_manual(values = pal, na.value = "grey95", name = "Kingdom") +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 15),
      strip.text.y = element_text(size = 15, face = "bold"),
      panel.spacing.x = unit(1.2, "lines"),
      panel.spacing.y = unit(0.8, "lines"),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20)
    )
}


col_coi_body <- (
  plot_site_marker("STO", "COI") /
    plot_site_marker("CCO", "COI") /
    plot_site_marker("CLR", "COI")
)

col_18s_body <- (
  plot_site_marker("STO", "18S") /
    plot_site_marker("CCO", "18S") /
    plot_site_marker("CLR", "18S")
)

col_coi <- wrap_elements(
  full = col_coi_body +
    plot_annotation(title = "COI") &
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.58))
)

col_18s <- wrap_elements(
  full = col_18s_body +
    plot_annotation(title = "18S") &
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.58))
)

p_leg <- plot_site_marker("STO", "COI") +
  theme(legend.position = "right")  # aquí SÍ la queremos

leg <- cowplot::get_legend(p_leg)

panel <- (col_coi | col_18s) +
  plot_annotation(caption = "Biological Replicates") &
  theme(
    plot.caption = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.caption.position = "plot"
  )

p_all <- cowplot::ggdraw() +
  cowplot::draw_label(
    "Relative read abundance",
    x = 0.02, y = 0.5, angle = 90,
    size = 20, fontface = "bold"
  ) +
  cowplot::draw_plot(panel, x = 0.06, y = 0.02, width = 0.80, height = 0.96) +
  cowplot::draw_plot(leg,   x = 0.87, y = 0.20, width = 0.12, height = 0.60)

p_all



 #############
#PLOT PHYLUM METAZOA
# COI
arel_mitjCOI_m<-arel_mitjCOI_kp%>%
  filter(kingdom=="Metazoa")

arel_mitjCOI_m_long<- pivot_longer(arel_mitjCOI_m,
                          cols = matches("^CCO|^CLR|^STO"),  
                          names_to = "Sample",
                          values_to = "ReadsCOI")

abdCOI_m_meta <- left_join(arel_mitjCOI_m_long, metadata, by = "Sample")
abdCOI_m_meta$Depth <- factor(abdCOI_m_meta$Depth)
abdCOI_m_meta$Site <- factor(abdCOI_m_meta$Site)
abdCOI_m_meta$Core_Rep <- factor(abdCOI_m_meta$Core_Rep)

# 18S
arel_mitj18S_m<-arel_mitj18S_kp%>%
  filter(kingdom=="Metazoa")

arel_mitj18S_m_long<- pivot_longer(arel_mitj18S_m,
                                   cols = matches("^CCO|^CLR|^STO"),  
                                   names_to = "Sample",
                                   values_to = "Reads18S")

abd18S_m_meta <- left_join(arel_mitj18S_m_long, metadata, by = "Sample")
abd18S_m_meta$Depth <- factor(abd18S_m_meta$Depth)
abd18S_m_meta$Site <- factor(abd18S_m_meta$Site)
abd18S_m_meta$Core_Rep <- factor(abd18S_m_meta$Core_Rep)



library(RColorBrewer)
palette1 <- brewer.pal(12, "Set3")
palette2 <- brewer.pal(8, "Set1")  

combined_palette <- c(palette1, palette2)


#STO COI
metaplotSTO_COI<- abdCOI_m_meta %>%
  filter(Site == "STO") %>%  
  ggplot(aes(x = Core_Rep, y = ReadsCOI, fill = phylum)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Depth) +  
  scale_fill_manual(values = combined_palette) +
  labs(
    title = "COI",  
    x = "",
    y = "Relative Abundance",
    fill = "Phylum"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    panel.spacing.x = unit(1.5, "lines")
  )

#STO 18S
metaplotSTO_18S<- abd18S_m_meta %>%
  filter(Site == "STO") %>%  
  ggplot(aes(x = Core_Rep, y = Reads18S, fill = phylum)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Depth) +  
  scale_fill_manual(values = combined_palette) +
  labs(
    title = "18S", 
    x = "Core",
    y = "Relative Abundance",
    fill = "Phylum"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    panel.spacing.x = unit(1.5, "lines")
  )
  
metacombined_STO <- (metaplotSTO_COI / metaplotSTO_18S) + 
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Relative Abundance by Phylum (Metazoa) - Locality STO (COI & 18S)",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      legend.position = "right",  # OR "bottom", your choice
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 12)
    )
  )

metacombined_STO
  

#CCO COI
metaplotCCO_COI<- abdCOI_m_meta %>%
  filter(Site == "CCO") %>%  
  ggplot(aes(x = Core_Rep, y = ReadsCOI, fill = phylum)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Depth) + 
  scale_fill_manual(values = combined_palette) +
  labs(
    title = "COI",  
    x = "",
    y = "Relative Abundance",
    fill = "Phylum"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    panel.spacing.x = unit(1.5, "lines")
  )

#CCO 18S
metaplotCCO_18S<- abd18S_m_meta %>%
  filter(ID == "CCO") %>%  
  ggplot(aes(x = Core_Rep, y = Reads18S, fill = phylum)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Depth) +  
  scale_fill_manual(values = combined_palette) +
  labs(
    title = "18S", 
    x = "Core",
    y = "Relative Abundance",
    fill = "Phylum"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    panel.spacing.x = unit(1.5, "lines")
  )

metacombined_CCO <- (metaplotCCO_COI / metaplotCCO_18S) + 
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Relative Abundance by Phylum (Metazoa) - Locality CCO",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      legend.position = "right",  # OR "bottom", your choice
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 12)
    )
  )

metacombined_CCO



#CLR COI
metaplotCLR_COI<- abdCOI_m_meta %>%
  filter(Site == "CLR") %>%  
  ggplot(aes(x = Core_Rep, y = ReadsCOI, fill = phylum)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Depth) +  
  scale_fill_manual(values = combined_palette) +
  labs(
    title = "COI",  
    x = "",
    y = "Relative Abundance",
    fill = "Phylum"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    panel.spacing.x = unit(1.5, "lines")
  )

#CLR 18S
metaplotCLR_18S<- abd18S_m_meta %>%
  filter(ID == "CLR") %>%  
  ggplot(aes(x = Core_Rep, y = Reads18S, fill = phylum)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Depth) +  
  scale_fill_manual(values = combined_palette) +
  labs(
    title = "18S",  
    x = "Core",
    y = "Relative Abundance",
    fill = "Phylum"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    panel.spacing.x = unit(1.5, "lines")
  )

metacombined_CLR <- (metaplotCLR_COI / metaplotCLR_18S) + 
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Relative Abundance by Phylum (Metazoa) - Locality CLR",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      legend.position = "right",  # OR "bottom", your choice
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 12)
    )
  )

metacombined_CLR




##### All Together METAZOA

df_meta <- bind_rows(
  abdCOI_meta %>%
    filter(kingdom == "Metazoa") %>%   
    transmute(
      Site      = Site,
      Core_Rep  = as.factor(Core_Rep),
      Depth     = Depth,
      Phylum    = phylum,              
      Marker    = "COI",
      RelAbund  = ReadsCOI,
      Year = Year
    ),
  abd18S_meta %>%
    filter(kingdom == "Metazoa") %>%
    transmute(
      Site      = Site,
      Core_Rep  = as.factor(Core_Rep),
      Depth     = Depth,
      Phylum    = phylum,
      Marker    = "18S",
      RelAbund  = Reads18S,
      Year = Year
    )
)

df_meta <- df_meta %>%
  mutate(
    Site     = factor(Site, levels = c("STO", "CCO", "CLR")),
    FacetCol = paste0(Marker, " – ", Year)
  )



phylum_order <- df_meta %>%
  group_by(Phylum) %>%
  summarise(Total = sum(RelAbund, na.rm = TRUE)) %>%
  arrange(desc(Total)) %>%
  pull(Phylum)

df_meta <- df_meta %>%
  mutate(Phylum = factor(Phylum, levels = phylum_order))


n_phy <- length(levels(df_meta$Phylum))

pal_phylum <- scales::hue_pal()(n_phy)


p_metazoa <- ggplot(
  df_meta,
  aes(x = Core_Rep, y = RelAbund, fill = Phylum)
) +
  geom_col(width = 0.85) +
  facet_grid(
    rows = vars(Site),
    cols = vars(FacetCol),
    switch = "y"
  ) +
  scale_fill_manual(values = pal_phylum, name = "Phylum") +
  labs(
    title = "Metazoan community composition by phylum across sites, markers and depths",
    x = "Core",
    y = "Relative abundance within Metazoa"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text.x = element_text(size = 12, face = "bold"),
    strip.text.y = element_text(size = 12, face = "bold"),
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.spacing.x = unit(1.2, "lines"),
    panel.spacing.y = unit(0.8, "lines")
  )

p_metazoa



df_all_clean_meta <- df_meta %>% 
  group_by(Site, FacetCol, Core_Rep) %>% 
  filter(sum(RelAbund, na.rm = TRUE) > 0) %>% 
  ungroup()


pal_phylum <- c(
  "#C6A300", 
  "#8F5E2D", 
  "#5B7F2B", 
  "#3F6E91", 
  "#B07C28", 
  "#C44E52", 
  "#8172B2", 
  "#55A868", 
  "#D18F59", 
  "#8C8C3E", 
  "#6B4C9A", 
  "#4C8A76", 
  "#AF7AA1", 
  "#7F7F7F", 
  "#A56A3A", 
  "#D0B240" 
)


colors_depth_site <- c(
  "COI - 2014"  = "#DCA0A0",
  "COI – 1939" = "#8B3A3A",
  "COI - 2017"  = "#66C2C2",
  "COI - 1952" = "cyan4",
  "COI - 1974"  = "#E6D46A",
  "COI - 1524" = "#CDAD00",
  "18S - 2014"  = "#DCA0A0",
  "18S - 1939" = "#8B3A3A",
  "18S - 2017"  = "#66C2C2",
  "18S - 1952" = "cyan4",
  "18S - 1974"  = "#E6D46A",
  "18S - 1524" = "#CDAD00"
)

df_all_clean_meta2 <- df_all_clean_meta %>%
  mutate(
    FacetCol_clean = gsub("[–—]", "-", FacetCol),  
    FacetCol_clean = trimws(FacetCol_clean),       
    FacetCol_clean = factor(FacetCol_clean, levels = names(colors_depth_site))
  )


setdiff(unique(df_all_clean_meta2$FacetCol_clean), names(colors_depth_site))
setdiff(levels(df_all_clean_meta2$FacetCol_clean), names(colors_depth_site))
setdiff(names(colors_depth_site), levels(df_all_clean_meta2$FacetCol_clean))


levels_facet <- levels(df_all_clean_meta2$FacetCol_clean)

assign_check <- data.frame(
  Facet = levels_facet,
  Year  = sub(".*-\\s*", "", levels_facet),
  stringsAsFactors = FALSE
)

assign_check$Color <- year_cols[assign_check$Year]

assign_check

assign_check$Year <- trimws(assign_check$Year)
assign_check$Year <- gsub("[^0-9]", "", assign_check$Year)  
assign_check$Color <- year_cols[assign_check$Year]

assign_check


colors_depth_site2 <- setNames(assign_check$Color, assign_check$Facet)

colors_depth_site2[grepl("1974", names(colors_depth_site2))]


plot_site <- function(site_id, add_x_label = FALSE, add_y_label = FALSE) {
  df_all_clean_meta2 %>% 
    filter(Site == site_id) %>% 
    ggplot(aes(x = Core_Rep, y = RelAbund, fill = Phylum)) +
    geom_col(width = 0.85) +
    facet_grid2(
      rows   = vars(Site),
      cols   = vars(FacetCol_clean),
      switch = "y",
      scales = "free_x",
      space  = "free_x",
      strip  = strip_themed(
        background_x = elem_list_rect(fill = colors_depth_site2, colour = "black"),
        text_x       = elem_list_text(colour = "black", face = "bold", size = 15)
      )
    ) +
    scale_fill_manual(values = pal_phylum, na.value = "grey95", name = "Phylum") +
    labs(
      x = if (add_x_label) "Core" else NULL,
      y = if (add_y_label) "Relative read abundance" else NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      panel.grid.major.x = element_blank(),
      legend.position = "right",
      axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 15),
      strip.text.x = element_text(size = 15, face = "bold"),
      strip.text.y = element_text(size = 15, face = "bold"),
      panel.spacing.x = unit(1.2, "lines"),
      panel.spacing.y = unit(0.8, "lines"),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20)
    )
}



p_sto_meta <- plot_site("STO")
p_cco_meta <- plot_site("CCO", add_y_label = TRUE)
p_clr_meta <- plot_site("CLR", add_x_label = TRUE)

p_all_meta <- p_sto_meta / p_cco_meta / p_clr_meta +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18),   
    legend.text  = element_text(size = 16),   
    legend.key.size = unit(1.2, "lines")      
  )
p_all_meta
p_all



library(dplyr)

df_all_clean_meta2 <- df_all_clean_meta %>%
  mutate(
    FacetCol_clean = FacetCol,
    FacetCol_clean = gsub("[–—]", "-", FacetCol_clean),  
    FacetCol_clean = gsub("\\s+", " ", FacetCol_clean),  
    FacetCol_clean = trimws(FacetCol_clean)
  )

levels_facet <- sort(unique(df_all_clean_meta2$FacetCol_clean))
df_all_clean_meta2$FacetCol_clean <- factor(df_all_clean_meta2$FacetCol_clean,
                                            levels = levels_facet)
year_cols <- c(
  "2014" = "#DCA0A0",
  "1939" = "#8B3A3A",
  "2017" = "#66C2C2",
  "1952" = "cyan4",
  "1974" = "#E6D46A",
  "1524" = "#CDAD00"
)

colors_depth_site2 <- setNames(
  year_cols[sub(".*-\\s*", "", levels_facet)],  
  levels_facet
)

colors_depth_site2[is.na(colors_depth_site2)]



df_all_clean_meta2 <- df_all_clean_meta2 %>%
  mutate(
    FacetCol_clean = factor(FacetCol_clean, levels = assign_check$Facet)
  )

fill_list <- as.list(colors_depth_site2[levels(df_all_clean_meta2$FacetCol_clean)])

stopifnot(length(fill_list) == nlevels(df_all_clean_meta2$FacetCol_clean))

plot_site <- function(site_id, add_x_label = FALSE, add_y_label = FALSE) {
  
  df_site <- df_all_clean_meta2 %>%
    filter(Site == site_id) %>%
    mutate(FacetCol_clean = droplevels(FacetCol_clean))  
  
  levs_site <- levels(df_site$FacetCol_clean)
  
  fill_list_site <- as.list(colors_depth_site2[levs_site])
  
  if (any(is.na(unlist(fill_list_site)))) {
    print(levs_site)
    print(colors_depth_site2[levs_site])
    stop("Hay algún FacetCol_clean sin color asignado (NA).")
  }
  
  ggplot(df_site, aes(x = Core_Rep, y = RelAbund, fill = Phylum)) +
    geom_col(width = 0.85) +
    facet_grid2(
      rows   = vars(Site),
      cols   = vars(FacetCol_clean),
      switch = "y",
      scales = "free_x",
      space  = "free_x",
      strip  = strip_themed(
        background_x = elem_list_rect(
          fill   = fill_list_site,  
          colour = "black"
        ),
        text_x = elem_list_text(
          colour = "black",
          face   = "bold",
          size   = 15
        )
      )
    ) +
    scale_fill_manual(values = pal_phylum, na.value = "grey95", name = "Phylum") +
    labs(
      x = if (add_x_label) "Biological Replicate" else NULL,
      y = if (add_y_label) "Relative read abundance" else NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      panel.grid.major.x = element_blank(),
      legend.position = "right",
      axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 15),
      strip.text.y = element_text(size = 15, face = "bold"),
      panel.spacing.x = unit(1.2, "lines"),
      panel.spacing.y = unit(0.8, "lines"),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20)
    )
}

p_sto_meta <- plot_site("STO")
p_cco_meta <- plot_site("CCO", add_y_label = TRUE)
p_clr_meta <- plot_site("CLR", add_x_label = TRUE)

p_all_meta <- p_sto_meta / p_cco_meta / p_clr_meta +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18),   
    legend.text  = element_text(size = 16),   
    legend.key.size = unit(1.2, "lines")      
  )
p_all_meta



df_all_clean_meta2 <- df_all_clean_meta2 %>%
  mutate(
    Gene = factor(Marker, levels = c("COI","18S")),
    FacetID   = FacetCol_clean,                               
    YearLabel = str_extract(as.character(FacetID), "\\d{4}")   
  )
df_all_clean_meta2 <- df_all_clean_meta2 %>%
  mutate(
    Marker = factor(Marker, levels = c("18S","COI")),
    YearLabel = factor(as.character(Year)),
    AgeGroup  = factor(ifelse(Depth %in% c("4","5"), "Recent", "Old"),
                       levels = c("Old","Recent"))
  )

colors_year <- c(
  "1939" = "#8B3A3A",
  "2014" = "#DCA0A0",
  "1952" = "cyan4",
  "2017" = "#66C2C2",
  "1524" = "#CDAD00",
  "1974" = "#E6D46A"
)

df_all_clean_meta2 <- df_all_clean_meta2 %>%
  mutate(
    Marker    = factor(Marker, levels = c("18S","COI")),
    YearLabel = factor(as.character(Year), levels = names(colors_year)),
    Site      = factor(Site, levels = c("STO","CCO","CLR"))
  )

fill_year_list <- as.list(colors_year[levels(df_all_clean_meta2$YearLabel)])

df_all_clean_meta2 <- df_all_clean_meta2 %>%
  group_by(Marker, Site, YearLabel) %>%
  filter(sum(RelAbund) > 0) %>%   
  ungroup()

plot_marker <- function(marker_id, show_y = TRUE) {
  
  df_all_clean_meta2 %>%
    filter(Marker == marker_id) %>%
    ggplot(aes(x = Core_Rep, y = RelAbund, fill = Phylum)) +
    
    geom_col(width = 0.85) +
    
    facet_grid2(
      rows   = vars(Site),
      cols   = vars(YearLabel),
      scales = "free_x",
      space  = "free_x",
      switch = "y",
      drop   = TRUE,
      strip  = strip_themed(
        background_x = elem_list_rect(
          fill   = fill_year_list,
          colour = "black"
        ),
        text_x = elem_list_text(
          colour = "black",
          face   = "bold",
          size   = 13
        )
      )
    ) +
    
    scale_fill_manual(values = pal_phylum, na.value = "grey95", name = "Phylum") +
    
    labs(
      x = "Core",
      y = if (show_y) "Relative read abundance" else NULL,
      title = marker_id
    ) +
    
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      strip.text.y = element_text(size = 13, face = "bold"),
      strip.text.x = element_text(size = 13, face = "bold"),
      legend.position = "right"
    )
}


p_18S <- plot_marker("18S", show_y = TRUE)
p_COI <- plot_marker("COI", show_y = FALSE)

p_final <- p_18S | p_COI +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14),
    legend.key.size = unit(1.1, "lines")
  )

p_final


##### Separated genes

plot_site_marker <- function(site_id, marker_id,
                             add_x_label = FALSE, add_y_label = FALSE) {
  
  df_site <- df_all_clean_meta2 %>%
    filter(Site == site_id, Marker == marker_id) %>%
    mutate(
      Year = factor(Year, levels = names(year_cols)),
      Year = droplevels(Year)
    )
  
  levs_year <- levels(df_site$Year)
  fill_list_year <- as.list(year_cols[levs_year])
  
  if (any(is.na(unlist(fill_list_year)))) {
    print(levs_year)
    print(year_cols[levs_year])
    stop("Hay algún Year sin color asignado (NA).")
  }
  
  ggplot(df_site, aes(x = Core_Rep, y = RelAbund, fill = Phylum)) +
    geom_col(width = 0.85) +
    facet_grid2(
      rows   = vars(Site),
      cols   = vars(Year),
      switch = "y",
      scales = "free_x",
      space  = "free_x",
      strip  = strip_themed(
        background_x = elem_list_rect(
          fill   = fill_list_year,   
          colour = "black"
        ),
        text_x = elem_list_text(
          colour = "black",
          face   = "bold",
          size   = 15
        )
      )
    ) +
    scale_fill_manual(values = pal_phylum, na.value = "grey95", name = "Phylum") +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 15),
      strip.text.y = element_text(size = 15, face = "bold"),
      panel.spacing.x = unit(1.2, "lines"),
      panel.spacing.y = unit(0.8, "lines"),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20)
    )
}


col_coi_body <- (
  plot_site_marker("STO", "COI") /
    plot_site_marker("CCO", "COI") /
    plot_site_marker("CLR", "COI")
)

col_18s_body <- (
  plot_site_marker("STO", "18S") /
    plot_site_marker("CCO", "18S") /
    plot_site_marker("CLR", "18S")
)

col_coi <- wrap_elements(
  full = col_coi_body +
    plot_annotation(title = "COI") &
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.58))
)

col_18s <- wrap_elements(
  full = col_18s_body +
    plot_annotation(title = "18S") &
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.58))
)

p_leg <- plot_site_marker("STO", "COI") +
  theme(legend.position = "right")  # aquí SÍ la queremos

leg <- cowplot::get_legend(p_leg)

panel <- (col_coi | col_18s) +
  plot_annotation(caption = "Biological Replicates") &
  theme(
    plot.caption = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.caption.position = "plot"
  )

p_all <- cowplot::ggdraw() +
  cowplot::draw_label(
    "Relative read abundance",
    x = 0.02, y = 0.5, angle = 90,
    size = 20, fontface = "bold"
  ) +
  cowplot::draw_plot(panel, x = 0.06, y = 0.02, width = 0.80, height = 0.96) +
  cowplot::draw_plot(leg,   x = 0.87, y = 0.20, width = 0.12, height = 0.60)

p_all


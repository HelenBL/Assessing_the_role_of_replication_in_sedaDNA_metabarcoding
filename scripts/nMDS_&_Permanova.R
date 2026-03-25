#NMDS & PERMANOVAS

library(tidyverse)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(permute)
library(vegan)
library(plotly)
library(geometry)   
library(dplyr)
library(patchwork)

### COI
metadata<- read.xlsx("metadata_f.xlsx")
arel_mitjCOI<- read.xlsx("Clean_Data_COI_AbRel.xlsx")
rownames(arel_mitjCOI) <- arel_mitjCOI$X1

arel_mitjCOI1 <- arel_mitjCOI[,c(4:21)]

arel_mitjCOI1_t<- as.data.frame(t(arel_mitjCOI1))

arel_mitjCOI1_t$Sample<- rownames(arel_mitjCOI1_t)

arel_mitjCOI1_t[, 1:(ncol(arel_mitjCOI1_t)-1)] <- lapply(arel_mitjCOI1_t[, 1:(ncol(arel_mitjCOI1_t)-1)], as.numeric)


any(is.na(arel_mitjCOI1_t))     
which(is.na(arel_mitjCOI1_t))

dist_COI <- vegdist(arel_mitjCOI1_t[, -ncol(arel_mitjCOI1_t)], method = "bray")

any(is.na(dist_COI))     
which(is.na(dist_COI))   


#k=2 dimensions
nmds2_COI_result<- metaMDS(dist_COI, k = 2, trymax = 100)

nmds2_COI_result$stress

nmds2_COI_points <- as.data.frame(nmds2_COI_result$points)
nmds2_COI_points$Sample<- rownames(nmds2_COI_points)
nmds2_COI_points<- merge(nmds2_COI_points, metadata, by = "Sample")
nmds2_COI_points$Depth<- as.factor(nmds2_COI_points$Depth) #paso findaries a categories

colors <- c("CCO" = "#8B3A3A", "CLR" = "cyan4", "STO" = "#CDAD00")

stress2_label_COI <- bquote(bold("Stress:") ~ .(round(nmds2_COI_result$stress, 4)))
stress2_label_COI <- paste0("Stress = ", round(nmds2_COI_result$stress, 3))

ggplot(nmds2_COI_points, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = Site, shape = Depth), size = 4, alpha = 0.9) +
  stat_ellipse(
    aes(fill = Site, group = Site),
    geom = "polygon",
    alpha = 0.2,
    color = NA,
    show.legend = FALSE
  ) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = c(15, 16, 17, 18, 3, 4)) +
  theme_classic(base_size = 14) +
  labs(
    title = "NMDS by Sample (COI)",
    color = "Location",
    shape = "Depth",
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 10)
  )+
  annotate("text",
           x = -2, y = 4,
           label = stress2_label_COI,
           hjust = 1.1, vjust = -0.5,
           size = 4.5)


#k=3 dimensions
nmds3_COI_result<- metaMDS(dist_COI, k = 3, trymax = 100)


nmds3_COI_result$stress

nmds3_COI_points <- as.data.frame(nmds3_COI_result$points)
nmds3_COI_points$Sample<- rownames(nmds3_COI_points)
nmds3_COI_points<- merge(nmds3_COI_points, metadata, by = "Sample")
nmds3_COI_points$Depth<- as.factor(nmds3_COI_points$Depth) #paso findaries a categories

colors <- c("CCO" = "#8B3A3A", "CLR" = "cyan4", "STO" = "#CDAD00")


stress3_label_COI <- bquote(bold("Stress:") ~ .(round(nmds3_COI_result$stress, 4)))
stress3_label_COI <- paste0("Stress = ", round(nmds3_COI_result$stress, ))

ggplot(nmds3_COI_points, aes(x = MDS1, y = MDS3)) +
  geom_point(aes(color = Site, shape = Depth), size = 4, alpha = 0.9) +
  stat_ellipse(
    aes(fill = Site, group = Site),
    geom = "polygon",
    alpha = 0.2,
    color = NA,
    show.legend = FALSE
  ) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = c(15, 16, 17, 18, 3, 4)) +
  theme_classic(base_size = 14) +
  labs(
    title = "NMDS by Sample (COI)",
    color = "Location",
    shape = "Depth",
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 10)
  )+
  annotate("text",
           x = -2, y = 4,
           label = stress3_label_COI,
           hjust = 1.1, vjust = -0.5,
           size = 4.5)

base_nmds3 <- function(x, y, title_suffix) {
  ggplot(nmds3_COI_points, aes_string(x = x, y = y)) +
    geom_point(aes(color = Site, shape = Depth), size = 4, alpha = 0.9) +
    stat_ellipse(
      aes(fill = Site, group = Site),
      geom = "polygon",
      alpha = 0.18,
      color = NA,
      show.legend = FALSE
    ) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_shape_manual(values = c(15, 16, 17, 18, 3, 4)) +
    theme_classic(base_size = 20) +
    labs(
      title = paste("NMDS (COI)", title_suffix),
      color = "Site",
      shape = "Depth"
    ) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text  = element_text(size = 9)
    )
}

p12_coi <- base_nmds3("MDS1", "MDS2", "axes 1–2")
p13_coi <- base_nmds3("MDS1", "MDS3", "axes 1–3")
p23 <- base_nmds3("MDS2", "MDS3", "axes 2–3")

p_nmds3_coi <- (p12_coi | p13_coi | p23) +
  plot_annotation(
    title = paste0("3D NMDS of COI communities (", stress3_label_18S, ")"),
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  )

p_nmds3_coi

#### 3D Plot with Plotly
# k = 3 dimensiones
nmds3_COI_result <- metaMDS(dist_COI, k = 3, trymax = 100)

nmds3_COI_result$stress

nmds3_COI_points <- as.data.frame(nmds3_COI_result$points)
nmds3_COI_points$Sample <- rownames(nmds3_COI_points)
nmds3_COI_points <- merge(nmds3_COI_points, metadata, by = "Sample")
nmds3_COI_points$Depth <- as.factor(nmds3_COI_points$Depth)

colors <- c("CCO" = "#8B3A3A", "CLR" = "cyan4", "STO" = "#CDAD00")

symbol_map_depth <- c(
  "4"  = "square",
  "5"  = "circle",
  "45" = "triangle-up",
  "50" = "diamond"
)
nmds3_COI_points$symbol <- symbol_map_depth[nmds3_COI_points$Depth]


colors<- c("CCO" = "darkolivegreen3", "CLR" = "cyan4", "STO" = "darkgoldenrod1")


fill_colors <- c("CCO" = "darkolivegreen3", "CLR" = "cyan4", "STO" = "darkgoldenrod1")

fill_colors_alpha <- paste0(fill_colors, "55")  


make_hull <- function(df, site){
  df_site <- df %>% filter(Site == site)
  pts <- as.matrix(df_site[, c("MDS1", "MDS2", "MDS3")])
  if (nrow(pts) < 4) {
    warning("Site ", site, ": not enough points to build a 3D hull")
    return(NULL)
  }
  ch <- convhulln(pts)        
  if (is.null(dim(ch))) {
    if (length(ch) %% 3 != 0) {
      stop("Convex hull for site ", site, " has unexpected length.")
    }
    ch <- matrix(ch, ncol = 3, byrow = TRUE)
  }
  list(
    x = pts[, 1],
    y = pts[, 2],
    z = pts[, 3],
    i = ch[, 1] - 1,
    j = ch[, 2] - 1,
    k = ch[, 3] - 1,
    type     = "mesh3d",
    color    = fill_colors_alpha[site],
    name     = paste(site, "hull"),
    opacity  = 0.25,
    showscale = FALSE
  )
}


hull_CCO <- make_hull(nmds3_COI_points, "CCO")
hull_CLR <- make_hull(nmds3_COI_points, "CLR")
hull_STO <- make_hull(nmds3_COI_points, "STO")


nmds3_COI_points$MarkerSymbol <- as.character(symbol_map_depth[nmds3_COI_points$Depth])
unique(nmds3_COI_points$MarkerSymbol)


p <- plot_ly()

p <- p %>% add_mesh(
  x = hull_CCO$x, y = hull_CCO$y, z = hull_CCO$z,
  i = hull_CCO$i, j = hull_CCO$j, k = hull_CCO$k,
  color = colors, opacity = 0.25,
  name = "CCO hull", showscale = FALSE
)

p <- p %>% add_mesh(
  x = hull_CLR$x, y = hull_CLR$y, z = hull_CLR$z,
  i = hull_CLR$i, j = hull_CLR$j, k = hull_CLR$k,
  color = colors, opacity = 0.25,
  name = "CLR hull", showscale = FALSE
)

p <- p %>% add_mesh(
  x = hull_STO$x, y = hull_STO$y, z = hull_STO$z,
  i = hull_STO$i, j = hull_STO$j, k = hull_STO$k,
  color = colors, opacity = 0.25,
  name = "STO hull", showscale = FALSE
)

for (d in unique(nmds3_COI_points$Depth)) {
  
  df_sub <- nmds3_COI_points[nmds3_COI_points$Depth == d, ]
  
  p <- p %>% add_markers(
    data = df_sub,
    x = ~MDS1, y = ~MDS2, z = ~MDS3,
    color = ~Site,
    colors = colors,
    symbol = I(symbol_map_depth[as.character(d)]),
    marker = list(size = 8, opacity = 0.9),
    name = paste("Depth", d)     
  )
}

p <- p %>% layout(
  title = paste0("3D NMDS (COI, ", stress3_label_COI, ")"),
  scene = list(
    xaxis = list(title = "MDS1"),
    yaxis = list(title = "MDS2"),
    zaxis = list(title = "MDS3")
  ),
  legend = list(
    title = list(text = "Groups"),
    font = list(family = "Arial")
  )
)

p




########## 18S

###AMB 18SSSS
arel_mitj18S<- read.xlsx("Clean_Data_18S_AbRel.xlsx")

rownames(arel_mitj18S) <- arel_mitj18S$X1
arel_mitj18S1 <- arel_mitj18S[,c(4:21)]

arel_mitj18S1_t<- as.data.frame(t(arel_mitj18S1))

arel_mitj18S1_t$Sample<- rownames(arel_mitj18S1_t)

arel_mitj18S1_t[, 1:(ncol(arel_mitj18S1_t)-1)] <- lapply(arel_mitj18S1_t[, 1:(ncol(arel_mitj18S1_t)-1)], as.numeric)

dist_18S <- vegdist(arel_mitj18S1_t[, -ncol(arel_mitj18S1_t)], method = "bray")


#k=2 dimensions
nmds2_18S_result<- metaMDS(dist_18S, k = 2, trymax = 100)


nmds2_18S_result$stress

nmds2_18S_points <- as.data.frame(nmds2_18S_result$points)
nmds2_18S_points$Sample<- rownames(nmds2_18S_points)
nmds2_18S_points<- merge(nmds2_18S_points, metadata, by = "Sample")
nmds2_18S_points$Depth<- as.factor(nmds2_18S_points$Depth) 

#colors<- c("CCO" = "darkolivegreen3", "CLR" = "cyan4", "STO" = "darkgoldenrod1")
colors <- c("CCO" = "#8B3A3A", "CLR" = "cyan4", "STO" = "#CDAD00")

stress2_label_18S <- bquote(bold("Stress:") ~ .(round(nmds2_18S_result$stress, 4)))
stress2_label_18S <- paste0("Stress = ", round(nmds2_18S_result$stress, ))

ggplot(nmds2_18S_points, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = Site, shape = Depth), size = 4, alpha = 0.9) +
  stat_ellipse(
    aes(fill = Site, group = Site),
    geom = "polygon",
    alpha = 0.2,
    color = NA,
    show.legend = FALSE
  ) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = c(15, 16, 17, 18, 3, 4)) +
  theme_classic(base_size = 14) +
  labs(
    title = "NMDS by Sample (18S)",
    color = "Location",
    shape = "Depth",
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 10)
  )+
  annotate("text",
           x = -2, y = 4,
           label = stress2_label_18S,
           hjust = 1.1, vjust = -0.5,
           size = 4.5)


#k=3 dimensions
nmds3_18S_result<- metaMDS(dist_18S, k = 3, trymax = 100)

#stress nmds 3 dim
nmds3_18S_result$stress

nmds3_18S_points <- as.data.frame(nmds3_18S_result$points)
nmds3_18S_points$Sample<- rownames(nmds3_18S_points)
nmds3_18S_points<- merge(nmds3_18S_points, metadata, by = "Sample")
nmds3_18S_points$Depth<- as.factor(nmds3_18S_points$Depth) 

colors <- c("CCO" = "#8B3A3A", "CLR" = "cyan4", "STO" = "#CDAD00")

stress3_label_18S <- bquote(bold("Stress:") ~ .(round(nmds3_18S_result$stress, 4)))
stress3_label_18S <- paste0("Stress = ", round(nmds3_18S_result$stress, ))

ggplot(nmds3_18S_points, aes(x = MDS2, y = MDS3)) +
  geom_point(aes(color = Site, shape = Depth), size = 4, alpha = 0.9) +
  stat_ellipse(
    aes(fill = Site, group = Site),
    geom = "polygon",
    alpha = 0.2,
    color = NA,
    show.legend = FALSE
  ) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = c(15, 16, 17, 18, 3, 4)) +
  theme_classic(base_size = 14) +
  labs(
    title = "NMDS by Sample (18S)",
    color = "Location",
    shape = "Depth",
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 10)
  )+
  annotate("text",
           x = -1, y = 4,
           label = stress3_label_18S,
           hjust = 1.1, vjust = -0.5,
           size = 4.5)

base_nmds3 <- function(x, y, title_suffix) {
  ggplot(nmds3_18S_points, aes_string(x = x, y = y)) +
    geom_point(aes(color = Site, shape = Depth), size = 4, alpha = 0.9) +
    stat_ellipse(
      aes(fill = Site, group = Site),
      geom = "polygon",
      alpha = 0.18,
      color = NA,
      show.legend = FALSE
    ) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_shape_manual(values = c(15, 16, 17, 18, 3, 4)) +
    theme_classic(base_size = 20) +
    labs(
      title = paste("NMDS (18S)", title_suffix),
      color = "Site",
      shape = "Depth"
    ) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text  = element_text(size = 9)
    )
}
p12_coi + p12_18
p12_18 <- base_nmds3("MDS1", "MDS2", "axes 1–2")
p13_18 <- base_nmds3("MDS1", "MDS3", "axes 1–3")
p23 <- base_nmds3("MDS2", "MDS3", "axes 2–3")

p_nmds3_18S <- (p12_18 | p13_18 | p23) +
  plot_annotation(
    title = paste0("3D NMDS of 18S communities (", stress3_label_18S, ")"),
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  )

p_nmds3_18S



#### 3 D with plotly

colors<- c("CCO" = "darkolivegreen3", "CLR" = "cyan4", "STO" = "darkgoldenrod1")

hull_CCO_18 <- make_hull(nmds3_18S_points, "CCO")
hull_CLR_18 <- make_hull(nmds3_18S_points, "CLR")
hull_STO_18 <- make_hull(nmds3_18S_points, "STO")


nmds3_18S_points$MarkerSymbol <- as.character(symbol_map_depth[nmds3_18S_points$Depth])
unique(nmds3_18S_points$MarkerSymbol)


p <- plot_ly()

p <- p %>% add_mesh(
  x = hull_CCO_18$x, y = hull_CCO_18$y, z = hull_CCO_18$z,
  i = hull_CCO_18$i, j = hull_CCO_18$j, k = hull_CCO_18$k,
  color = colors, opacity = 0.25,
  name = "CCO hull", showscale = FALSE
)

p <- p %>% add_mesh(
  x = hull_CLR_18$x, y = hull_CLR_18$y, z = hull_CLR_18$z,
  i = hull_CLR_18$i, j = hull_CLR_18$j, k = hull_CLR_18$k,
  color = colors, opacity = 0.25,
  name = "CLR hull", showscale = FALSE
)

p <- p %>% add_mesh(
  x = hull_STO_18$x, y = hull_STO_18$y, z = hull_STO_18$z,
  i = hull_STO_18$i, j = hull_STO_18$j, k = hull_STO_18$k,
  color = colors, opacity = 0.25,
  name = "STO hull", showscale = FALSE
)

for (d in unique(nmds3_18S_points$Depth)) {
  
  df_sub <- nmds3_18S_points[nmds3_18S_points$Depth == d, ]
  
  p <- p %>% add_markers(
    data = df_sub,
    x = ~MDS1, y = ~MDS2, z = ~MDS3,
    color = ~Site,
    colors = colors,
    symbol = I(symbol_map_depth[as.character(d)]),
    marker = list(size = 8, opacity = 0.9),
    name = paste("Depth", d)     
  )
}

p <- p %>% layout(
  title = paste0("3D NMDS (18S, ", stress3_label_18S, ")"),
  scene = list(
    xaxis = list(title = "MDS1"),
    yaxis = list(title = "MDS2"),
    zaxis = list(title = "MDS3")
  ),
  legend = list(
    title = list(text = "Groups"),
    font = list(family = "Arial")
  )
)

p



#### PERMANOVA (location)

perm_COI_terms <- adonis2(
  dist_COI ~ Site + Depth + Core_Rep,
  data = metadata,
  permutations = 999,
  by = "margin"
)

perm_18S_terms <- adonis2(
  dist_18S ~ Site + Depth + Core_Rep,
  data = metadata,
  permutations = 999,
  by = "margin"
)

perm_COI_terms
perm_18S_terms

### DISPERSION (betadisper + permutest)

factors_to_check <- c("Site", "Depth", "Core_Rep")

disp_COI <- run_dispersion_tests(dist_COI, metadata, factors_to_check, nperm = 999)
disp_18S <- run_dispersion_tests(dist_18S, metadata, factors_to_check, nperm = 999)

disp_COI$Site$permutest
disp_COI$Depth$permutest
disp_COI$Core_Rep$permutest

disp_18S$Site$permutest
disp_18S$Depth$permutest
disp_18S$Core_Rep$permutest


pairwise_dispersion <- function(betadisp_obj, nperm = 999) {
  permutest(betadisp_obj, permutations = nperm, pairwise = TRUE)
}

pairwise_dispersion(disp_COI$Site$betadisper, nperm = 999)

plot_dispersion <- function(betadisp_obj, group_name = "Group") {
  df <- data.frame(
    group = betadisp_obj$group,
    dist_to_centroid = betadisp_obj$distances
  )
  ggplot(df, aes(x = group, y = dist_to_centroid)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.6) +
    theme_classic(base_size = 14) +
    labs(x = group_name, y = "Distance to centroid (multivariate dispersion)")
}

plot_dispersion(disp_COI$Site$betadisper, "Site (Site)")


# COI
r2_COI <- data.frame(
  Gene   = "COI",
  Factor = rownames(perm_COI_terms),
  R2     = perm_COI_terms$R2 * 100
)

# 18S
r2_18S <- data.frame(
  Gene   = "18S",
  Factor = rownames(perm_18S_terms),
  R2     = perm_18S_terms$R2 * 100
)

r2_all <- rbind(r2_COI, r2_18S) |>
  dplyr::filter(Factor %in% c("Site", "Depth", "Core_Rep", "Residual"))

ggplot(r2_all, aes(x = Factor, y = R2, fill = Gene)) +
  geom_col(position = "dodge", width = 0.7) +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("COI" = "#996fe1ff", "18S" = "#5eb273ff")) +
  scale_x_discrete(labels = c(
    Site   = "Site",
    Depth    = "Depth horizon",
    Core_Rep = "Biological replicate",
    Residual = "Unexplained variation"
  )) +
  labs(
    y = "Variance explained (%)",
    x = "",
    fill = "Marker",
  ) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_text(size = 14),
    
    axis.text.x  = element_text(size = 15),
    axis.text.y  = element_text(size = 15),
    
    axis.title.x = element_text(size = 17, face = "bold"),
    axis.title.y = element_text(size = 17, face = "bold")
  )






####### METAZOA

##### COI
arel_mitjCOI <- arel_mitjCOI %>%
  filter(arel_mitjCOI$kingdom == "Metazoa")
arel_mitjCOI1 <- arel_mitjCOI[,c(4:21)]

arel_mitjCOI1_t<- as.data.frame(t(arel_mitjCOI1))

arel_mitjCOI1_t$Sample<- rownames(arel_mitjCOI1_t)

arel_mitjCOI1_t[, 1:(ncol(arel_mitjCOI1_t)-1)] <- lapply(arel_mitjCOI1_t[, 1:(ncol(arel_mitjCOI1_t)-1)], as.numeric)


any(is.na(arel_mitjCOI1_t))     
which(is.na(arel_mitjCOI1_t))

dist_COI <- vegdist(arel_mitjCOI1_t[, -ncol(arel_mitjCOI1_t)], method = "bray")

any(is.na(dist_COI))     
which(is.na(dist_COI))   


#k=2 dimensions
nmds2_COI_result<- metaMDS(dist_COI, k = 2, trymax = 100)

nmds2_COI_result$stress

nmds2_COI_points <- as.data.frame(nmds2_COI_result$points)
nmds2_COI_points$Sample<- rownames(nmds2_COI_points)
nmds2_COI_points<- merge(nmds2_COI_points, metadata, by = "Sample")
nmds2_COI_points$Depth<- as.factor(nmds2_COI_points$Depth) 

colors <- c("CCO" = "#8B3A3A", "CLR" = "cyan4", "STO" = "#CDAD00")

stress2_label_COI <- bquote(bold("Stress:") ~ .(round(nmds2_COI_result$stress, 4)))
stress2_label_COI <- paste0("Stress = ", round(nmds2_COI_result$stress, 3))

p12_coi<- ggplot(nmds2_COI_points, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = Site, shape = Depth), size = 4, alpha = 0.9) +
  stat_ellipse(
    aes(fill = Site, group = Site),
    geom = "polygon",
    alpha = 0.2,
    color = NA,
    show.legend = FALSE
  ) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = c(15, 16, 17, 18, 3, 4)) +
  theme_classic(base_size = 14) +
  labs(
    title = "NMDS by Sample (COI)",
    color = "Location",
    shape = "Depth",
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 10)
  )+
  annotate("text",
           x = -2, y = 4,
           label = stress2_label_COI,
           hjust = 1.1, vjust = -0.5,
           size = 4.5)
p12_coi

#k=3 dimensions
nmds3_COI_result<- metaMDS(dist_COI, k = 3, trymax = 100)


nmds3_COI_result$stress

nmds3_COI_points <- as.data.frame(nmds3_COI_result$points)
nmds3_COI_points$Sample<- rownames(nmds3_COI_points)
nmds3_COI_points<- merge(nmds3_COI_points, metadata, by = "Sample")
nmds3_COI_points$Depth<- as.factor(nmds3_COI_points$Depth) 

colors <- c("CCO" = "#8B3A3A", "CLR" = "cyan4", "STO" = "#CDAD00")


stress3_label_COI <- bquote(bold("Stress:") ~ .(round(nmds3_COI_result$stress, 4)))
stress3_label_COI <- paste0("Stress = ", round(nmds3_COI_result$stress, ))

ggplot(nmds3_COI_points, aes(x = MDS1, y = MDS3)) +
  geom_point(aes(color = Site, shape = Depth), size = 4, alpha = 0.9) +
  stat_ellipse(
    aes(fill = Site, group = Site),
    geom = "polygon",
    alpha = 0.2,
    color = NA,
    show.legend = FALSE
  ) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = c(15, 16, 17, 18, 3, 4)) +
  theme_classic(base_size = 14) +
  labs(
    title = "NMDS by Sample (COI)",
    color = "Location",
    shape = "Depth",
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 10)
  )+
  annotate("text",
           x = -2, y = 4,
           label = stress3_label_COI,
           hjust = 1.1, vjust = -0.5,
           size = 4.5)

base_nmds3 <- function(x, y, title_suffix) {
  ggplot(nmds3_COI_points, aes_string(x = x, y = y)) +
    geom_point(aes(color = Site, shape = Depth), size = 4, alpha = 0.9) +
    stat_ellipse(
      aes(fill = Site, group = Site),
      geom = "polygon",
      alpha = 0.18,
      color = NA,
      show.legend = FALSE
    ) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_shape_manual(values = c(15, 16, 17, 18, 3, 4)) +
    theme_classic(base_size = 20) +
    labs(
      title = paste("NMDS (COI)", title_suffix),
      color = "Site",
      shape = "Depth"
    ) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text  = element_text(size = 9)
    )
}

p12_coi <- base_nmds3("MDS1", "MDS2", "axes 1–2")
p13_coi <- base_nmds3("MDS1", "MDS3", "axes 1–3")
p23 <- base_nmds3("MDS2", "MDS3", "axes 2–3")

p_nmds3_coi <- (p12_coi | p13_coi | p23) +
  plot_annotation(
    title = paste0("3D NMDS of COI communities (", stress3_label_18S, ")"),
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  )

p_nmds3_coi


#### 3D Plot with Plotly
# k = 3 dimensiones
nmds3_COI_result <- metaMDS(dist_COI, k = 3, trymax = 100)

nmds3_COI_result$stress

nmds3_COI_points <- as.data.frame(nmds3_COI_result$points)
nmds3_COI_points$Sample <- rownames(nmds3_COI_points)
nmds3_COI_points <- merge(nmds3_COI_points, metadata, by = "Sample")
nmds3_COI_points$Depth <- as.factor(nmds3_COI_points$Depth)

symbol_map_depth <- c(
  "4"  = "square",
  "5"  = "circle",
  "45" = "triangle-up",
  "50" = "diamond"
)
nmds3_COI_points$symbol <- symbol_map_depth[nmds3_COI_points$Depth]


colors<- c("CCO" = "darkolivegreen3", "CLR" = "cyan4", "STO" = "darkgoldenrod1")


fill_colors <- c("CCO" = "darkolivegreen3", "CLR" = "cyan4", "STO" = "darkgoldenrod1")

fill_colors_alpha <- paste0(fill_colors, "55")  


make_hull <- function(df, site){
  df_site <- df %>% filter(Site == site)
  pts <- as.matrix(df_site[, c("MDS1", "MDS2", "MDS3")])
  if (nrow(pts) < 4) {
    warning("Site ", site, ": not enough points to build a 3D hull")
    return(NULL)
  }
  ch <- convhulln(pts)        
  if (is.null(dim(ch))) {
    if (length(ch) %% 3 != 0) {
      stop("Convex hull for site ", site, " has unexpected length.")
    }
    ch <- matrix(ch, ncol = 3, byrow = TRUE)
  }
  list(
    x = pts[, 1],
    y = pts[, 2],
    z = pts[, 3],
    i = ch[, 1] - 1,
    j = ch[, 2] - 1,
    k = ch[, 3] - 1,
    type     = "mesh3d",
    color    = fill_colors_alpha[site],
    name     = paste(site, "hull"),
    opacity  = 0.25,
    showscale = FALSE
  )
}


hull_CCO <- make_hull(nmds3_COI_points, "CCO")
hull_CLR <- make_hull(nmds3_COI_points, "CLR")
hull_STO <- make_hull(nmds3_COI_points, "STO")


nmds3_COI_points$MarkerSymbol <- as.character(symbol_map_depth[nmds3_COI_points$Depth])
unique(nmds3_COI_points$MarkerSymbol)


p <- plot_ly()

p <- p %>% add_mesh(
  x = hull_CCO$x, y = hull_CCO$y, z = hull_CCO$z,
  i = hull_CCO$i, j = hull_CCO$j, k = hull_CCO$k,
  color = colors, opacity = 0.25,
  name = "CCO hull", showscale = FALSE
)

p <- p %>% add_mesh(
  x = hull_CLR$x, y = hull_CLR$y, z = hull_CLR$z,
  i = hull_CLR$i, j = hull_CLR$j, k = hull_CLR$k,
  color = colors, opacity = 0.25,
  name = "CLR hull", showscale = FALSE
)

p <- p %>% add_mesh(
  x = hull_STO$x, y = hull_STO$y, z = hull_STO$z,
  i = hull_STO$i, j = hull_STO$j, k = hull_STO$k,
  color = colors, opacity = 0.25,
  name = "STO hull", showscale = FALSE
)

for (d in unique(nmds3_COI_points$Depth)) {
  
  df_sub <- nmds3_COI_points[nmds3_COI_points$Depth == d, ]
  
  p <- p %>% add_markers(
    data = df_sub,
    x = ~MDS1, y = ~MDS2, z = ~MDS3,
    color = ~Site,
    colors = colors,
    symbol = I(symbol_map_depth[as.character(d)]),
    marker = list(size = 8, opacity = 0.9),
    name = paste("Depth", d)     
  )
}

p <- p %>% layout(
  title = paste0("3D NMDS (COI, ", stress3_label_COI, ")"),
  scene = list(
    xaxis = list(title = "MDS1"),
    yaxis = list(title = "MDS2"),
    zaxis = list(title = "MDS3")
  ),
  legend = list(
    title = list(text = "Groups"),
    font = list(family = "Arial")
  )
)

p





###### 18S
arel_mitj18S <- arel_mitj18S %>%
  filter(arel_mitj18S$kingdom == "Metazoa")
arel_mitj18S1 <- arel_mitj18S[,c(4:21)]

arel_mitj18S1_t<- as.data.frame(t(arel_mitj18S1))

arel_mitj18S1_t$Sample<- rownames(arel_mitj18S1_t)

arel_mitj18S1_t[, 1:(ncol(arel_mitj18S1_t)-1)] <- lapply(arel_mitj18S1_t[, 1:(ncol(arel_mitj18S1_t)-1)], as.numeric)

dist_18S <- vegdist(arel_mitj18S1_t[, -ncol(arel_mitj18S1_t)], method = "bray")


#k=2 dimensions
nmds2_18S_result<- metaMDS(dist_18S, k = 2, trymax = 100)


nmds2_18S_result$stress

nmds2_18S_points <- as.data.frame(nmds2_18S_result$points)
nmds2_18S_points$Sample<- rownames(nmds2_18S_points)
nmds2_18S_points<- merge(nmds2_18S_points, metadata, by = "Sample")
nmds2_18S_points$Depth<- as.factor(nmds2_18S_points$Depth) 

#colors<- c("CCO" = "darkolivegreen3", "CLR" = "cyan4", "STO" = "darkgoldenrod1")
colors <- c("CCO" = "#8B3A3A", "CLR" = "cyan4", "STO" = "#CDAD00")

stress2_label_18S <- bquote(bold("Stress:") ~ .(round(nmds2_18S_result$stress, 4)))
stress2_label_18S <- paste0("Stress = ", round(nmds2_18S_result$stress, ))

p12_18S<- ggplot(nmds2_18S_points, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = Site, shape = Depth), size = 4, alpha = 0.9) +
  stat_ellipse(
    aes(fill = Site, group = Site),
    geom = "polygon",
    alpha = 0.2,
    color = NA,
    show.legend = FALSE
  ) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = c(15, 16, 17, 18, 3, 4)) +
  theme_classic(base_size = 14) +
  labs(
    title = "NMDS by Sample (18S)",
    color = "Location",
    shape = "Depth",
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 10)
  )+
  annotate("text",
           x = -2, y = 4,
           label = stress2_label_18S,
           hjust = 1.1, vjust = -0.5,
           size = 4.5)

p12_18S

p12_coi + p12_18S




#k=3 dimensions
nmds3_18S_result<- metaMDS(dist_18S, k = 3, trymax = 100)

#stress nmds 3 dim
nmds3_18S_result$stress

nmds3_18S_points <- as.data.frame(nmds3_18S_result$points)
nmds3_18S_points$Sample<- rownames(nmds3_18S_points)
nmds3_18S_points<- merge(nmds3_18S_points, metadata, by = "Sample")
nmds3_18S_points$Depth<- as.factor(nmds3_18S_points$Depth) 

colors <- c("CCO" = "#8B3A3A", "CLR" = "cyan4", "STO" = "#CDAD00")

stress3_label_18S <- bquote(bold("Stress:") ~ .(round(nmds3_18S_result$stress, 4)))
stress3_label_18S <- paste0("Stress = ", round(nmds3_18S_result$stress, ))

ggplot(nmds3_18S_points, aes(x = MDS2, y = MDS3)) +
  geom_point(aes(color = Site, shape = Depth), size = 4, alpha = 0.9) +
  stat_ellipse(
    aes(fill = Site, group = Site),
    geom = "polygon",
    alpha = 0.2,
    color = NA,
    show.legend = FALSE
  ) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = c(15, 16, 17, 18, 3, 4)) +
  theme_classic(base_size = 14) +
  labs(
    title = "NMDS by Sample (18S)",
    color = "Location",
    shape = "Depth",
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 10)
  )+
  annotate("text",
           x = -1, y = 4,
           label = stress3_label_18S,
           hjust = 1.1, vjust = -0.5,
           size = 4.5)

base_nmds3 <- function(x, y, title_suffix) {
  ggplot(nmds3_18S_points, aes_string(x = x, y = y)) +
    geom_point(aes(color = Site, shape = Depth), size = 4, alpha = 0.9) +
    stat_ellipse(
      aes(fill = Site, group = Site),
      geom = "polygon",
      alpha = 0.18,
      color = NA,
      show.legend = FALSE
    ) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_shape_manual(values = c(15, 16, 17, 18, 3, 4)) +
    theme_classic(base_size = 20) +
    labs(
      title = paste("NMDS (18S)", title_suffix),
      color = "Site",
      shape = "Depth"
    ) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text  = element_text(size = 9)
    )
}
p12_coi + p12_18
p12_18 <- base_nmds3("MDS1", "MDS2", "axes 1–2")
p13_18 <- base_nmds3("MDS1", "MDS3", "axes 1–3")
p23 <- base_nmds3("MDS2", "MDS3", "axes 2–3")

p_nmds3_18S <- (p12_18 | p13_18 | p23) +
  plot_annotation(
    title = paste0("3D NMDS of 18S communities (", stress3_label_18S, ")"),
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  )

p_nmds3_18S



#### 3 D with plotly

colors<- c("CCO" = "darkolivegreen3", "CLR" = "cyan4", "STO" = "darkgoldenrod1")

hull_CCO_18 <- make_hull(nmds3_18S_points, "CCO")
hull_CLR_18 <- make_hull(nmds3_18S_points, "CLR")
hull_STO_18 <- make_hull(nmds3_18S_points, "STO")


nmds3_18S_points$MarkerSymbol <- as.character(symbol_map_depth[nmds3_18S_points$Depth])
unique(nmds3_18S_points$MarkerSymbol)


p <- plot_ly()

p <- p %>% add_mesh(
  x = hull_CCO_18$x, y = hull_CCO_18$y, z = hull_CCO_18$z,
  i = hull_CCO_18$i, j = hull_CCO_18$j, k = hull_CCO_18$k,
  color = colors, opacity = 0.25,
  name = "CCO hull", showscale = FALSE
)

p <- p %>% add_mesh(
  x = hull_CLR_18$x, y = hull_CLR_18$y, z = hull_CLR_18$z,
  i = hull_CLR_18$i, j = hull_CLR_18$j, k = hull_CLR_18$k,
  color = colors, opacity = 0.25,
  name = "CLR hull", showscale = FALSE
)

p <- p %>% add_mesh(
  x = hull_STO_18$x, y = hull_STO_18$y, z = hull_STO_18$z,
  i = hull_STO_18$i, j = hull_STO_18$j, k = hull_STO_18$k,
  color = colors, opacity = 0.25,
  name = "STO hull", showscale = FALSE
)

for (d in unique(nmds3_18S_points$Depth)) {
  
  df_sub <- nmds3_18S_points[nmds3_18S_points$Depth == d, ]
  
  p <- p %>% add_markers(
    data = df_sub,
    x = ~MDS1, y = ~MDS2, z = ~MDS3,
    color = ~Site,
    colors = colors,
    symbol = I(symbol_map_depth[as.character(d)]),
    marker = list(size = 8, opacity = 0.9),
    name = paste("Depth", d)     
  )
}

p <- p %>% layout(
  title = paste0("3D NMDS (18S, ", stress3_label_18S, ")"),
  scene = list(
    xaxis = list(title = "MDS1"),
    yaxis = list(title = "MDS2"),
    zaxis = list(title = "MDS3")
  ),
  legend = list(
    title = list(text = "Groups"),
    font = list(family = "Arial")
  )
)

p





################### PERMAQNOVA -- Metazoa

perm_COI_terms <- adonis2(
  dist_COI ~ Site + Depth + Core_Rep,
  data = metadata,
  permutations = 999,
  by = "margin"
)

perm_18S_terms <- adonis2(
  dist_18S ~ Site + Depth + Core_Rep,
  data = metadata,
  permutations = 999,
  by = "margin"
)

perm_COI_terms
perm_18S_terms

### DISPERSION (betadisper + permutest)

factors_to_check <- c("Site", "Depth", "Core_Rep")

disp_COI <- run_dispersion_tests(dist_COI, metadata, factors_to_check, nperm = 999)
disp_18S <- run_dispersion_tests(dist_18S, metadata, factors_to_check, nperm = 999)

disp_COI$Site$permutest
disp_COI$Depth$permutest
disp_COI$Core_Rep$permutest

disp_18S$Site$permutest
disp_18S$Depth$permutest
disp_18S$Core_Rep$permutest


pairwise_dispersion <- function(betadisp_obj, nperm = 999) {
  permutest(betadisp_obj, permutations = nperm, pairwise = TRUE)
}

pairwise_dispersion(disp_COI$Site$betadisper, nperm = 999)

plot_dispersion <- function(betadisp_obj, group_name = "Group") {
  df <- data.frame(
    group = betadisp_obj$group,
    dist_to_centroid = betadisp_obj$distances
  )
  ggplot(df, aes(x = group, y = dist_to_centroid)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.6) +
    theme_classic(base_size = 14) +
    labs(x = group_name, y = "Distance to centroid (multivariate dispersion)")
}

plot_dispersion(disp_COI$Site$betadisper, "Site (Site)")


# COI
r2_COI <- data.frame(
  Gene   = "COI",
  Factor = rownames(perm_COI_terms),
  R2     = perm_COI_terms$R2 * 100
)

# 18S
r2_18S <- data.frame(
  Gene   = "18S",
  Factor = rownames(perm_18S_terms),
  R2     = perm_18S_terms$R2 * 100
)

r2_all <- rbind(r2_COI, r2_18S) |>
  dplyr::filter(Factor %in% c("Site", "Depth", "Core_Rep", "Residual"))

ggplot(r2_all, aes(x = Factor, y = R2, fill = Gene)) +
  geom_col(position = "dodge", width = 0.7) +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("COI" = "#996fe1ff", "18S" = "#5eb273ff")) +
  scale_x_discrete(labels = c(
    Site   = "Site",
    Depth    = "Depth horizon",
    Core_Rep = "Biological replicate",
    Residual = "Unexplained variation"
  )) +
  labs(
    y = "Variance explained (%)",
    x = "",
    fill = "Marker",
  ) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_text(size = 14),
    
    axis.text.x  = element_text(size = 15),
    axis.text.y  = element_text(size = 15),
    
    axis.title.x = element_text(size = 17, face = "bold"),
    axis.title.y = element_text(size = 17, face = "bold")
  )


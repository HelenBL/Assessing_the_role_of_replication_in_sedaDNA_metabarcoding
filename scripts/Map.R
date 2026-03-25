# ### Mapping using ggOceanmaps
# ##https://mikkovihtakari.github.io/ggOceanMaps/articles/ggOceanMaps.html
# 
library(stars)
library(terra)
library(ggOceanMaps)
library(sf)
library(ggspatial)
library(raster)
library(dplyr)
library(sp)
library(raster)
library(ncdf4)
library(RColorBrewer)
library(sf)
library(tmap)
library(geodata)
library(raster)
library(sf)
library(ggplot2)
library(ggspatial)
library(lattice)
library(marmap)
library(leaflet)
library(leaflet.extras)
library(maptiles)
library(ggspatial)


metadata_loc <- readxl::read_excel(file.choose())

metadata_loc$Lon <- gsub("[–—−]", "-", metadata_loc$Lon)
metadata_loc$Lon <- as.numeric(metadata_loc$Lon)
metadata_sf <- sf::st_as_sf(metadata_loc, coords = c("Lon", "Lat"), crs = 4326)


######## Include Topography 


alt <- elevation_30s(country="Spain", path=tempdir())
alt_r <- raster(alt)



land_raster <- alt_r > 0  


land_poly <- rasterToPolygons(land_raster, dissolve = TRUE)
land_sf <- st_as_sf(land_poly)

e <- extent(-10, 5, 35, 45)
e_sp <- as(e, "SpatialPolygons")
crs(e_sp) <- crs(alt_r)
alt_r <- crop(alt_r, e_sp)



slope <- terrain(alt_r, "slope")
plot(slope)
aspect <- terrain(alt_r, "aspect")
plot(aspect)


hill2 <- hillShade(slope, aspect, angle = 40, direction = 270)
plot(hill2)




bathy <- getNOAA.bathy(lon1 = -10, lon2 = 5, lat1 = 35, lat2 = 45, resolution = 0.01)

bathy_df <- fortify.bathy(bathy)


bathy_r <- bathy_df
bathy_r$z[bathy_r$z >= 0] <- NA  

relief_r <- bathy_df
relief_r$z[relief_r$z < 0] <- NA  

bathy_df <- as.data.frame(bathy_r, xy = TRUE)
colnames(bathy_df) <- c("x", "y", "depth")
bathy_df <- na.omit(bathy_df)

relief_df <- as.data.frame(relief_r, xy = TRUE)
colnames(relief_df) <- c("x", "y", "elevation")
relief_df <- na.omit(relief_df)


bathy_palette <- colorRampPalette(c("royalblue4", "royalblue3", "steelblue2", "lightskyblue1"))
relief_palette <- colorRampPalette(c("grey20","grey40", "grey70", "grey90"))



bathy_df$color <- bathy_palette(100)[
  cut(bathy_df$depth, breaks = 100, labels = FALSE)
]

relief_df$color <- relief_palette(100)[
  cut(relief_df$elevation, breaks = 100, labels = FALSE)
]

colors<- c("CCO" = "darkolivegreen3", "CLR" = "cyan4", "STO" = "darkgoldenrod1")


ggplot() +
  geom_tile(data = bathy_df, aes(x = x, y = y, fill = I(color))) +
  geom_tile(data = relief_df, aes(x = x, y = y, fill = I(color)), alpha = 0.8) +
  geom_sf(data = land_sf, fill = NA, color = "black", size = 1) +
  coord_sf(xlim = c(-10, 5), ylim = c(35, 45), expand = FALSE)+
  labs(
    x = "",
    y = "")+
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14, face="bold"))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  annotation_scale(
    location = "bl",
    style = "bar",     
    pad_x = unit(0.2, "cm"),
    pad_y = unit(0.2, "cm")
  )+
  annotation_north_arrow(
    location = "tl", which_north = "true",
    style = north_arrow_fancy_orienteering,
    height = unit(1, "cm"), width = unit(1, "cm"),
    pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm")
  )
  




bbox_soton <- st_as_sfc(
  st_bbox(c(xmin = -3.55, xmax = -3.35, ymin = 43.36, ymax = 43.48), crs = 4326)
)

tiles <- get_tiles(bbox_soton, provider = "CartoDB.Positron", crop = TRUE, zoom = 13)
tiles_df <- as.data.frame(tiles, xy = TRUE)
names(tiles_df)[3:5] <- c("R", "G", "B")


muestreo_sf <- st_as_sf(data.frame(
  site = c("CCO", "CLR", "STO"),
  lon = c(-3.47851, -3.45429, -3.46607),
  lat = c(43.4170167, 43.4060055, 43.446357),
  tipo = c("CCO", "CLR", "STO")
), coords = c("lon", "lat"), crs = 4326)

rectangulo <- st_as_sfc(
  st_bbox(c(xmin = -3.55, xmax = -3.35, ymin = 43.36, ymax = 43.48), crs = 4326)
)

ggplot() +
  geom_raster(data = tiles_df, aes(x = x, y = y,
                                   fill = rgb(R, G, B, maxColorValue = 255))) +
  scale_fill_identity() +
  geom_sf(data = muestreo_sf, aes(color = tipo), size = 8) +
  scale_color_manual(values = c("CCO" = "darkolivegreen3", "CLR" = "cyan4", "STO"= "goldenrod1")) +
  geom_sf(data = rectangulo, fill = NA, color = "black", size = 2) +
  coord_sf(xlim = c(-3.55, -3.35), ylim = c(43.36, 43.48), expand = FALSE) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title = element_blank()) +
  labs(color = "Sampling Site", face= "bold") +
  annotation_north_arrow(
    location = "tl", which_north = "true",
    style = ggspatial::north_arrow_fancy_orienteering,
    height = unit(1.2, "cm"), width = unit(1.2, "cm"),
    pad_x = unit(0.5, "cm"), pad_y = unit(0.5, "cm")
  ) +
  annotation_scale(
    location = "bl", width_hint = 0.2,
    bar_cols = c("black", "white"), text_cex = 0.8, line_width = 1
  )


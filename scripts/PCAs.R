### PCAs

library(ape)       
library(ggplot2)
library(dplyr)

###COI
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


pcoa_COI <- pcoa(dist_COI)

pcoa_points <- as.data.frame(pcoa_COI$vectors[, 1:3])
colnames(pcoa_points) <- c("PCoA1", "PCoA2", "PCoA3")
pcoa_points$Sample <- rownames(pcoa_points)

var_expl <- 100 * pcoa_COI$values$Relative_eig[1:3]

pcoa_points <- merge(pcoa_points, metadata, by = "Sample")
pcoa_points$Depth <- as.factor(pcoa_points$Depth)

colors <- c("CCO" = "#8B3A3A", "CLR" = "cyan4", "STO" = "#CDAD00")

ggplot(pcoa_points, aes(x = PCoA1, y = PCoA3)) +
  geom_point(aes(color = Site, shape = Depth), size = 5, alpha = 0.9) +
  scale_color_manual(values = colors, name = "Site") +
  theme_classic(base_size = 14) +
  labs(
    title = "PCoA COI (Bray–Curtis)",
    x = paste0("PCoA1 (", round(var_expl[1], 1), "%)"),
    y = paste0("PCoA3 (", round(var_expl[3], 1), "%)"),
    shape = "Depth (cm)"
  ) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold")
  )

pc_coi<-ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = Site, shape = Depth), size = 5, alpha = 0.9) +
  scale_color_manual(values = colors, name = "Site") +
  theme_classic(base_size = 14) +
  labs(
    title = "PCoA COI (Bray–Curtis)",
    x = paste0("PCoA1 (", round(var_expl[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_expl[2], 1), "%)"),
    shape = "Depth (cm)"
  ) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold")
  )



###### 18S

arel_mitj18S<- read.xlsx("Clean_Data_18S_AbRel.xlsx")

rownames(arel_mitj18S) <- arel_mitj18S$X1
arel_mitj18S1 <- arel_mitj18S[,c(4:21)]

arel_mitj18S1_t<- as.data.frame(t(arel_mitj18S1))

arel_mitj18S1_t$Sample<- rownames(arel_mitj18S1_t)

arel_mitj18S1_t[, 1:(ncol(arel_mitj18S1_t)-1)] <- lapply(arel_mitj18S1_t[, 1:(ncol(arel_mitj18S1_t)-1)], as.numeric)

dist_18S <- vegdist(arel_mitj18S1_t[, -ncol(arel_mitj18S1_t)], method = "bray")

pcoa_18 <- pcoa(dist_18S)

pcoa_points_18 <- as.data.frame(pcoa_18$vectors[, 1:3])
colnames(pcoa_points_18) <- c("PCoA1","PCoA2", "PCoA3")
pcoa_points_18$Sample <- rownames(pcoa_points_18)

var_expl <- 100 * pcoa_18$values$Relative_eig[1:2]

pcoa_points_18 <- merge(pcoa_points_18, metadata, by = "Sample")
pcoa_points_18$Depth <- as.factor(pcoa_points_18$Depth)


pc_18<-ggplot(pcoa_points_18, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = Site, shape = Depth), size = 5, alpha = 0.9) +
  scale_color_manual(values = colors, name = "Site") +
  theme_classic(base_size = 14) +
  labs(
    title = "PCoA 18S (Bray–Curtis)",
    x = paste0("PCoA1 (", round(var_expl[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_expl[2], 1), "%)"),
    shape = "Depth (cm)"
  ) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold")
  )
pc_coi + pc_18





############METAZOA

#### COI

arel_mitjCOI <- arel_mitjCOI %>%
  filter(arel_mitjCOI$kingdom == "Metazoa")
arel_mitjCOI1 <- arel_mitjCOI[,c(4:21)]

arel_mitjCOI1_t<- as.data.frame(t(arel_mitjCOI1))

arel_mitjCOI1_t$Sample<- rownames(arel_mitjCOI1_t)

arel_mitjCOI1_t[, 1:(ncol(arel_mitjCOI1_t)-1)] <- lapply(arel_mitjCOI1_t[, 1:(ncol(arel_mitjCOI1_t)-1)], as.numeric)


any(is.na(arel_mitjCOI1_t))     
which(is.na(arel_mitjCOI1_t))

dist_COI <- vegdist(arel_mitjCOI1_t[, -ncol(arel_mitjCOI1_t)], method = "bray")

pcoa_COI <- pcoa(dist_COI)

pcoa_points <- as.data.frame(pcoa_COI$vectors[, 1:3])
colnames(pcoa_points) <- c("PCoA1", "PCoA2", "PCoA3")
pcoa_points$Sample <- rownames(pcoa_points)

var_expl <- 100 * pcoa_COI$values$Relative_eig[1:3]

pcoa_points <- merge(pcoa_points, metadata, by = "Sample")
pcoa_points$Depth <- as.factor(pcoa_points$Depth)


ggplot(pcoa_points, aes(x = PCoA1, y = PCoA3)) +
  geom_point(aes(color = Site, shape = Depth), size = 5, alpha = 0.9) +
  scale_color_manual(values = colors, name = "Site") +
  theme_classic(base_size = 14) +
  labs(
    title = "PCoA COI (Bray–Curtis)",
    x = paste0("PCoA1 (", round(var_expl[1], 1), "%)"),
    y = paste0("PCoA3 (", round(var_expl[3], 1), "%)"),
    shape = "Depth (cm)"
  ) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold")
  )

pc_coi<-ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = Site, shape = Depth), size = 5, alpha = 0.9) +
  scale_color_manual(values = colors, name = "Site") +
  theme_classic(base_size = 14) +
  labs(
    title = "PCoA COI (Bray–Curtis)",
    x = paste0("PCoA1 (", round(var_expl[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_expl[2], 1), "%)"),
    shape = "Depth (cm)"
  ) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold")
  )


### 18S

arel_mitj18S <- arel_mitj18S %>%
  filter(arel_mitj18S$kingdom == "Metazoa")
arel_mitj18S1 <- arel_mitj18S[,c(4:21)]

arel_mitj18S1_t<- as.data.frame(t(arel_mitj18S1))

arel_mitj18S1_t$Sample<- rownames(arel_mitj18S1_t)

arel_mitj18S1_t[, 1:(ncol(arel_mitj18S1_t)-1)] <- lapply(arel_mitj18S1_t[, 1:(ncol(arel_mitj18S1_t)-1)], as.numeric)

dist_18S <- vegdist(arel_mitj18S1_t[, -ncol(arel_mitj18S1_t)], method = "bray")

pcoa_18 <- pcoa(dist_18S)

pcoa_points_18 <- as.data.frame(pcoa_18$vectors[, 1:3])
colnames(pcoa_points_18) <- c("PCoA1","PCoA2", "PCoA3")
pcoa_points_18$Sample <- rownames(pcoa_points_18)

var_expl <- 100 * pcoa_18$values$Relative_eig[1:2]

pcoa_points_18 <- merge(pcoa_points_18, metadata, by = "Sample")
pcoa_points_18$Depth <- as.factor(pcoa_points_18$Depth)


pc_18<-ggplot(pcoa_points_18, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = Site, shape = Depth), size = 5, alpha = 0.9) +
  scale_color_manual(values = colors, name = "Site") +
  theme_classic(base_size = 14) +
  labs(
    title = "PCoA 18S (Bray–Curtis)",
    x = paste0("PCoA1 (", round(var_expl[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_expl[2], 1), "%)"),
    shape = "Depth (cm)"
  ) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold")
  )
pc_coi + pc_18

###dj-RDA

library(vegan)
library(ggplot2)
library(dplyr)

metadata<- read.xlsx("metadata_f.xlsx")
arel_mitjCOI<- read.xlsx("Clean_Data_COI_AbRel.xlsx")
arel_mitjCOI1 <- column_to_rownames(arel_mitjCOI, var = colnames(arel_mitjCOI)[1])

arel_mitjCOI1 <- arel_mitjCOI1[,c(3:20)]

arel_mitjCOI1_t<- as.data.frame(t(arel_mitjCOI1))

arel_mitjCOI1_t$Sample<- rownames(arel_mitjCOI1_t)

arel_mitjCOI1_t[, 1:(ncol(arel_mitjCOI1_t)-1)] <- lapply(arel_mitjCOI1_t[, 1:(ncol(arel_mitjCOI1_t)-1)], as.numeric)

any(is.na(arel_mitjCOI1_t))     
dist_COI <- vegdist(arel_mitjCOI1_t[, -ncol(arel_mitjCOI1_t)], method = "bray")
      

dbRDA_COI <- capscale(dist_COI ~ CoreID + Depth, data = metadata)

anova_global   <- anova(dbRDA_COI, permutations = 999)
anova_terms    <- anova(dbRDA_COI, permutations = 999, by = "terms")
anova_axes     <- anova(dbRDA_COI, permutations = 999, by = "axis")

anova_global
anova_terms
anova_axes

scores_sites <- as.data.frame(scores(dbRDA_COI, display = "sites"))
scores_sites$Sample <- rownames(scores_sites)

scores_sites <- merge(scores_sites, metadata, by = "Sample")


var_exp <- round(100 * dbRDA_COI$CCA$eig / sum(dbRDA_COI$CCA$eig), 1)
lab_x <- paste0("dbRDA1 (", var_exp[1], "%)")
lab_y <- paste0("dbRDA2 (", var_exp[2], "%)")


metadata$AgeGroup <- ifelse(metadata$Depth %in% c(4, 5), "recent", "old")
metadata$AgeNum   <- ifelse(metadata$AgeGroup == "recent", 1, -1)

metadata$RegionGroup <- ifelse(metadata$CoreID == "STO", "west", "east")
metadata$RegionNum   <- ifelse(metadata$RegionGroup == "west", 1, -1)

ef <- envfit(dbRDA_COI ~ AgeNum + RegionNum, data = metadata, permutations = 999)

vecs <- as.data.frame(scores(ef, display = "vectors"))
vecs   

arrow_df <- rbind(
  data.frame(
    label = "recent",
    CAP1  = vecs["AgeNum", "CAP1"],
    CAP2  = vecs["AgeNum", "CAP2"]
  ),
  data.frame(
    label = "old",
    CAP1  = -vecs["AgeNum", "CAP1"],
    CAP2  = -vecs["AgeNum", "CAP2"]
  ),
  data.frame(
    label = "west",
    CAP1  = vecs["RegionNum", "CAP1"],
    CAP2  = vecs["RegionNum", "CAP2"]
  ),
  data.frame(
    label = "east",
    CAP1  = -vecs["RegionNum", "CAP1"],
    CAP2  = -vecs["RegionNum", "CAP2"]
  )
)


colors <- c("CCO"="#8B3A3A", "CLR"="cyan4", "STO"="#CDAD00")

p_dbRDA <- ggplot(scores_sites, aes(x = CAP1, y = CAP2)) +
  geom_point(
    aes(color = CoreID, shape = as.factor(Depth)),
    size = 6, alpha = 0.9
  ) +
  stat_ellipse(
    aes(color = CoreID),
    type = "t", level = 0.68, linewidth = 1
  ) +
  scale_color_manual(values = colors, name = "Site") +
  scale_shape_manual(name = "Depth (cm)", values = c(16, 17, 15, 18)) +
  
  geom_segment(
    data = arrow_df,
    aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
    arrow = arrow(length = unit(0.25, "cm")),
    linewidth = 0.8,
    color = "black"
  ) +
  geom_text(
    data = arrow_df,
    aes(x = CAP1 * 1.1, y = CAP2 * 1.1, label = label),
    fontface = "bold",
    size = 6
  ) +
  
  labs(
    title = "COI",
    x = lab_x,
    y = lab_y
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold")
  )

p_dbRDA






################### 18S


arel_mitj18S<- read.xlsx("Clean_Data_18S_AbRel.xlsx")

arel_mitj18S1 <- column_to_rownames(arel_mitj18S, var = colnames(arel_mitj18S)[1])
arel_mitj18S1 <- arel_mitj18S1[,c(3:20)]

arel_mitj18S1_t<- as.data.frame(t(arel_mitj18S1))

arel_mitj18S1_t$Sample<- rownames(arel_mitj18S1_t)

arel_mitj18S1_t[, 1:(ncol(arel_mitj18S1_t)-1)] <- lapply(arel_mitj18S1_t[, 1:(ncol(arel_mitj18S1_t)-1)], as.numeric)

dist_18S <- vegdist(arel_mitj18S1_t[, -ncol(arel_mitj18S1_t)], method = "bray")



dbRDA_18s <- capscale(dist_18S ~ CoreID + Depth, data = metadata)

anova_global   <- anova(dbRDA_18s, permutations = 999)
anova_terms    <- anova(dbRDA_18s, permutations = 999, by = "terms")
anova_axes     <- anova(dbRDA_18s, permutations = 999, by = "axis")

anova_global
anova_terms
anova_axes

scores_sites_18S <- as.data.frame(scores(dbRDA_18s, display = "sites"))
scores_sites_18S$Sample <- rownames(scores_sites_18S)

scores_sites_18S <- merge(scores_sites_18S, metadata, by = "Sample")

var_exp <- round(100 * dbRDA_18s$CCA$eig / sum(dbRDA_18s$CCA$eig), 1)
lab_x <- paste0("dbRDA1 (", var_exp[1], "%)")
lab_y <- paste0("dbRDA2 (", var_exp[2], "%)")


metadata$AgeGroup <- ifelse(metadata$Depth %in% c(4, 5), "recent", "old")
metadata$AgeNum   <- ifelse(metadata$AgeGroup == "recent", 1, -1)

metadata$RegionGroup <- ifelse(metadata$CoreID == "STO", "west", "east")
metadata$RegionNum   <- ifelse(metadata$RegionGroup == "west", 1, -1)

ef <- envfit(dbRDA_18s ~ AgeNum + RegionNum, data = metadata, permutations = 999)

vecs <- as.data.frame(scores(ef, display = "vectors"))
vecs   

arrow_df <- rbind(
  data.frame(
    label = "recent",
    CAP1  = vecs["AgeNum", "CAP1"],
    CAP2  = vecs["AgeNum", "CAP2"]
  ),
  data.frame(
    label = "old",
    CAP1  = -vecs["AgeNum", "CAP1"],
    CAP2  = -vecs["AgeNum", "CAP2"]
  ),
  data.frame(
    label = "west",
    CAP1  = vecs["RegionNum", "CAP1"],
    CAP2  = vecs["RegionNum", "CAP2"]
  ),
  data.frame(
    label = "east",
    CAP1  = -vecs["RegionNum", "CAP1"],
    CAP2  = -vecs["RegionNum", "CAP2"]
  )
)



p_dbRDA_18S <- ggplot(scores_sites_18S, aes(x = CAP1, y = CAP2)) +
  geom_point(
    aes(color = CoreID, shape = as.factor(Depth)),
    size = 6, alpha = 0.9
  ) +
  stat_ellipse(
    aes(color = CoreID),
    type = "t", level = 0.68, linewidth = 1
  ) +
  scale_color_manual(values = colors, name = "Site") +
  scale_shape_manual(name = "Depth (cm)", values = c(16, 17, 15, 18)) +
  
  geom_segment(
    data = arrow_df,
    aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
    arrow = arrow(length = unit(0.25, "cm")),
    linewidth = 0.8,
    color = "black"
  ) +
  geom_text(
    data = arrow_df,
    aes(x = CAP1 * 1.1, y = CAP2 * 1.1, label = label),
    fontface = "bold",
    size = 6
  ) +
  
  labs(
    title = "18S",
    x = lab_x,
    y = lab_y
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold")
  )

p_dbRDA_18S

p_dbRDA + p_dbRDA_18S

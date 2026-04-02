library(readxl)
library(dplyr)
library(ggplot2)
library(vegan)
library(scales)
library(colorspace)
library(openxlsx) 
library (tibble)


#### COI 
reps_out <- read.xlsx("Clean_Data_COI_Reps_AbRel.xlsx")
metadata_reps <- read.xlsx ("metadata_repliques_f.xlsx")

rownames(reps_out) <- reps_out$X1
reps_out <- reps_out[,c(2:127)]


reps_t <- as.data.frame(t(reps_out))
reps_t$Sample <- rownames(reps_t)

reps_t[, 1:(ncol(reps_t)-1)] <-
  lapply(reps_t[, 1:(ncol(reps_t)-1)], as.numeric)

asv_mat_rep <- reps_t[, !(names(reps_t) %in% "Sample")]

alpha_rep <- data.frame(
  Sample   = reps_t$Sample,
  Richness = specnumber(asv_mat_rep),
  Shannon  = diversity(asv_mat_rep, index = "shannon")
)

alpha_rep_meta_clean <- merge(metadata_reps, alpha_rep, by = "Sample")



alpha_rep_meta_clean$X <- paste(alpha_rep_meta_clean$Site, alpha_rep_meta_clean$Core_Rep, sep = "_")
alpha_rep_meta_clean$X <- factor(alpha_rep_meta_clean$X, levels = unique(alpha_rep_meta_clean$X))


colors_site <- c("CCO" = "#8B3A3A", "CLR" = "cyan4", "STO" = "#CDAD00")

colors_depth_site <- c(
  "CCO_5"  = "#DCA0A0",  
  "CCO_45" = "#8B3A3A",  
  
  "CLR_4"  = "#66C2C2",  
  "CLR_45" = "cyan4",    
  
  "STO_5"  = "#E6D46A",  
  "STO_50" = "#CDAD00"  
)



ggplot(alpha_rep_meta_clean, aes(x = X, y = Richness)) +
  
  geom_boxplot(
    aes(fill = Site),
    alpha = 0.7,
    width = 0.6,
    outlier.shape = NA
  ) +
  
  geom_point(
    color = "black",     
    size = 2.3,
    position = position_jitter(width = 0.07), 
    alpha = 0.9
  ) +
  
  scale_fill_manual(values = colors_site, name = "Site") +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = "Site and biological replicate (core)",
    y = "ASV richness"
  ) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(fill = NA)
  )


alpha_rep_meta_clean$Depth  <- as.character(alpha_rep_meta_clean$Depth)
alpha_rep_meta_clean$Site <- as.character(alpha_rep_meta_clean$Site)

alpha_rep_meta_clean$SiteDepth <- paste0(alpha_rep_meta_clean$Site, "_", alpha_rep_meta_clean$Depth)
unique(alpha_rep_meta_clean$SiteDepth)

alpha_rep_meta_clean$X <- paste(alpha_rep_meta_clean$Site, alpha_rep_meta_clean$Core_Rep, sep = "_")
alpha_rep_meta_clean$X <- factor(alpha_rep_meta_clean$X, levels = unique(alpha_rep_meta_clean$X))

ggplot(alpha_rep_meta_clean, aes(x = X, y = Richness)) +
  
  geom_boxplot(
    aes(fill = SiteDepth),
    alpha = 0.8,
    width = 0.6,
    outlier.shape = NA
  ) +
  
  geom_point(aes(fill = SiteDepth),
             color = "black",
             size = 2.3,
             alpha = 0.9,
             position = position_dodge(width = 0.6)   
  ) +
  
  scale_fill_manual(values = colors_depth_site, name = "Site & depth") +
  
  theme_classic(base_size = 14) +
  labs(
    x = "Site and biological replicate (core)",
    y = "ASV richness"
  ) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.title  = element_text(face = "bold"),
    panel.border  = element_rect(fill = NA)
  )


ggplot(alpha_rep_meta_clean, aes(x = X, y = Shannon)) +
  
  geom_boxplot(
    aes(fill = SiteDepth),
    alpha = 0.8,
    width = 0.6,
    outlier.shape = NA
  ) +
  
  geom_point(aes(fill = SiteDepth),
             color = "black",
             size = 2.3,
             alpha = 0.9,
             position = position_dodge(width = 0.6)   
  ) +
  
  scale_fill_manual(values = colors_depth_site, name = "Site & depth") +
  
  theme_classic(base_size = 14) +
  labs(
    x = "Site and biological replicate (core)",
    y = "ASV Shannon"
  ) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.title  = element_text(face = "bold"),
    panel.border  = element_rect(fill = NA)
  )


#write.xlsx(alpha_rep_meta_clean, "alpha_diver_clean_coi.xlsx")





### 18S

reps_out <- read.xlsx("Clean_Data_18S_Reps_AbRel.xlsx")

rownames(reps_out) <- reps_out$X1
reps_out <- reps_out[,c(2:141)]

reps_t <- as.data.frame(t(reps_out))
reps_t$Sample <- rownames(reps_t)

reps_t[, 1:(ncol(reps_t)-1)] <-
  lapply(reps_t[, 1:(ncol(reps_t)-1)], as.numeric)

asv_mat_rep <- reps_t[, !(names(reps_t) %in% "Sample")]

alpha_rep <- data.frame(
  Sample   = reps_t$Sample,
  Richness = specnumber(asv_mat_rep),
  Shannon  = diversity(asv_mat_rep, index = "shannon")
)

alpha_rep_meta_clean <- merge(metadata_reps, alpha_rep, by = "Sample")


alpha_rep_meta_clean$X <- paste(alpha_rep_meta_clean$Site, alpha_rep_meta_clean$Core_Rep, sep = "_")
alpha_rep_meta_clean$X <- factor(alpha_rep_meta_clean$X, levels = unique(alpha_rep_meta_clean$X))



ggplot(alpha_rep_meta_clean, aes(x = X, y = Richness)) +
  
  geom_boxplot(
    aes(fill = Site),
    alpha = 0.7,
    width = 0.6,
    outlier.shape = NA
  ) +
  
  geom_point(
    color = "black",     
    size = 2.3,
    position = position_jitter(width = 0.07),  
    alpha = 0.9
  ) +
  
  scale_fill_manual(values = colors_site, name = "Site") +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = "Site and biological replicate (core)",
    y = "ASV richness"
  ) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(fill = NA)
  )


alpha_rep_meta_clean$Depth  <- as.character(alpha_rep_meta_clean$Depth)
alpha_rep_meta_clean$Site <- as.character(alpha_rep_meta_clean$Site)

alpha_rep_meta_clean$SiteDepth <- paste0(alpha_rep_meta_clean$Site, "_", alpha_rep_meta_clean$Depth)
unique(alpha_rep_meta_clean$SiteDepth)

alpha_rep_meta_clean$X <- paste(alpha_rep_meta_clean$Site, alpha_rep_meta_clean$Core_Rep, sep = "_")
alpha_rep_meta_clean$X <- factor(alpha_rep_meta_clean$X, levels = unique(alpha_rep_meta_clean$X))

ggplot(alpha_rep_meta_clean, aes(x = X, y = Richness)) +

    geom_boxplot(
    aes(fill = SiteDepth),
    alpha = 0.8,
    width = 0.6,
    outlier.shape = NA
  ) +
  
  geom_point(aes(fill = SiteDepth),
             color = "black",
             size = 2.3,
             alpha = 0.9,
             position = position_dodge(width = 0.6)   
  ) +
  
  scale_fill_manual(values = colors_depth_site, name = "Site & depth") +
  
  theme_classic(base_size = 14) +
  labs(
    x = "Site and biological replicate (core)",
    y = "ASV richness"
  ) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.title  = element_text(face = "bold"),
    panel.border  = element_rect(fill = NA)
  )


ggplot(alpha_rep_meta_clean, aes(x = X, y = Shannon)) +
  
  geom_boxplot(
    aes(fill = SiteDepth),
    alpha = 0.8,
    width = 0.6,
    outlier.shape = NA
  ) +
  
  geom_point(aes(fill = SiteDepth),
             color = "black",
             size = 2.3,
             alpha = 0.9,
             position = position_dodge(width = 0.6)   
  ) +
  
  scale_fill_manual(values = colors_depth_site, name = "Site & depth") +
  
  theme_classic(base_size = 14) +
  labs(
    x = "Site and biological replicate (core)",
    y = "ASV Shannon"
  ) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.title  = element_text(face = "bold"),
    panel.border  = element_rect(fill = NA)
  )


#write.xlsx(alpha_rep_meta_clean, "alpha_diver_clean_18s.xlsx")


################# METAZOA
reps_out <- read.xlsx("Clean_Data_COI_Reps_AbRel.xlsx")
tax_COI <- read.xlsx("Clean_Data_COI_AVS_taxonomy.xlsx")
colnames(reps_out)[colnames(reps_out) == "X1"] <- "ASV"

metazoa_reps_coi <- reps_out %>%
  left_join(tax_COI, by = "ASV") %>%
  filter(kingdom == "Metazoa")

rownames(metazoa_reps_coi) <- metazoa_reps_coi$ASV
metazoa_reps_coi <- metazoa_reps_coi[,c(2:127)]


reps_t <- as.data.frame(t(metazoa_reps_coi))
reps_t$Sample <- rownames(reps_t)

reps_t[, 1:(ncol(reps_t)-1)] <-
  lapply(reps_t[, 1:(ncol(reps_t)-1)], as.numeric)

asv_mat_rep <- reps_t[, !(names(reps_t) %in% "Sample")]

alpha_rep <- data.frame(
  Sample   = reps_t$Sample,
  Richness = specnumber(asv_mat_rep),
  Shannon  = diversity(asv_mat_rep, index = "shannon")
)

alpha_rep_meta_clean <- merge(metadata_reps, alpha_rep, by = "Sample")



alpha_rep_meta_clean$X <- paste(alpha_rep_meta_clean$Site, alpha_rep_meta_clean$Core_Rep, sep = "_")
alpha_rep_meta_clean$X <- factor(alpha_rep_meta_clean$X, levels = unique(alpha_rep_meta_clean$X))


colors_site <- c("CCO" = "#8B3A3A", "CLR" = "cyan4", "STO" = "#CDAD00")

colors_depth_site <- c(
  "CCO_5"  = "#DCA0A0",  
  "CCO_45" = "#8B3A3A",  
  
  "CLR_4"  = "#66C2C2",  
  "CLR_45" = "cyan4",    
  
  "STO_5"  = "#E6D46A",  
  "STO_50" = "#CDAD00"  
)



ggplot(alpha_rep_meta_clean, aes(x = X, y = Richness)) +
  
  geom_boxplot(
    aes(fill = Site),
    alpha = 0.7,
    width = 0.6,
    outlier.shape = NA
  ) +
  
  geom_point(
    color = "black",     
    size = 2.3,
    position = position_jitter(width = 0.07), 
    alpha = 0.9
  ) +
  
  scale_fill_manual(values = colors_site, name = "Site") +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = "Site and biological replicate (core)",
    y = "ASV richness"
  ) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(fill = NA)
  )


alpha_rep_meta_clean$Depth  <- as.character(alpha_rep_meta_clean$Depth)
alpha_rep_meta_clean$Site <- as.character(alpha_rep_meta_clean$Site)

alpha_rep_meta_clean$SiteDepth <- paste0(alpha_rep_meta_clean$Site, "_", alpha_rep_meta_clean$Depth)
unique(alpha_rep_meta_clean$SiteDepth)

alpha_rep_meta_clean$X <- paste(alpha_rep_meta_clean$Site, alpha_rep_meta_clean$Core_Rep, sep = "_")
alpha_rep_meta_clean$X <- factor(alpha_rep_meta_clean$X, levels = unique(alpha_rep_meta_clean$X))

ggplot(alpha_rep_meta_clean, aes(x = X, y = Richness)) +
  
  geom_boxplot(
    aes(fill = SiteDepth),
    alpha = 0.8,
    width = 0.6,
    outlier.shape = NA
  ) +
  
  geom_point(aes(fill = SiteDepth),
             color = "black",
             size = 2.3,
             alpha = 0.9,
             position = position_dodge(width = 0.6)   
  ) +
  
  scale_fill_manual(values = colors_depth_site, name = "Site & depth") +
  
  theme_classic(base_size = 14) +
  labs(
    x = "Site and biological replicate (core)",
    y = "ASV richness"
  ) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.title  = element_text(face = "bold"),
    panel.border  = element_rect(fill = NA)
  )


ggplot(alpha_rep_meta_clean, aes(x = X, y = Shannon)) +
  
  geom_boxplot(
    aes(fill = SiteDepth),
    alpha = 0.8,
    width = 0.6,
    outlier.shape = NA
  ) +
  
  geom_point(aes(fill = SiteDepth),
             color = "black",
             size = 2.3,
             alpha = 0.9,
             position = position_dodge(width = 0.6)   
  ) +
  
  scale_fill_manual(values = colors_depth_site, name = "Site & depth") +
  
  theme_classic(base_size = 14) +
  labs(
    x = "Site and biological replicate (core)",
    y = "ASV Shannon"
  ) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.title  = element_text(face = "bold"),
    panel.border  = element_rect(fill = NA)
  )



###18S
reps_18 <- read.xlsx("Clean_Data_18S_Reps_AbRel.xlsx")
colnames(reps_18)[colnames(reps_18) == "X1"] <- "ASV"
tax_18 <- read.xlsx("Clean_Data_18S_AVS_taxonomy.xlsx")

metazoa_reps_18 <- reps_18 %>%
  left_join(tax_18, by = "ASV") %>%
  filter(kingdom == "Metazoa")


rownames(metazoa_reps_18) <- metazoa_reps_18$ASV
metazoa_reps_18 <- metazoa_reps_18[,c(2:141)]

reps_t <- as.data.frame(t(metazoa_reps_18))
reps_t$Sample <- rownames(reps_t)

reps_t[, 1:(ncol(reps_t)-1)] <-
  lapply(reps_t[, 1:(ncol(reps_t)-1)], as.numeric)

asv_mat_rep <- reps_t[, !(names(reps_t) %in% "Sample")]

alpha_rep <- data.frame(
  Sample   = reps_t$Sample,
  Richness = specnumber(asv_mat_rep),
  Shannon  = diversity(asv_mat_rep, index = "shannon")
)

alpha_rep_meta_clean <- merge(metadata_reps, alpha_rep, by = "Sample")


alpha_rep_meta_clean$X <- paste(alpha_rep_meta_clean$Site, alpha_rep_meta_clean$Core_Rep, sep = "_")
alpha_rep_meta_clean$X <- factor(alpha_rep_meta_clean$X, levels = unique(alpha_rep_meta_clean$X))



ggplot(alpha_rep_meta_clean, aes(x = X, y = Richness)) +
  
  geom_boxplot(
    aes(fill = Site),
    alpha = 0.7,
    width = 0.6,
    outlier.shape = NA
  ) +
  
  geom_point(
    color = "black",     
    size = 2.3,
    position = position_jitter(width = 0.07),  
    alpha = 0.9
  ) +
  
  scale_fill_manual(values = colors_site, name = "Site") +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = "Site and biological replicate (core)",
    y = "ASV richness"
  ) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(fill = NA)
  )


alpha_rep_meta_clean$Depth  <- as.character(alpha_rep_meta_clean$Depth)
alpha_rep_meta_clean$Site <- as.character(alpha_rep_meta_clean$Site)

alpha_rep_meta_clean$SiteDepth <- paste0(alpha_rep_meta_clean$Site, "_", alpha_rep_meta_clean$Depth)
unique(alpha_rep_meta_clean$SiteDepth)

alpha_rep_meta_clean$X <- paste(alpha_rep_meta_clean$Site, alpha_rep_meta_clean$Core_Rep, sep = "_")
alpha_rep_meta_clean$X <- factor(alpha_rep_meta_clean$X, levels = unique(alpha_rep_meta_clean$X))

ggplot(alpha_rep_meta_clean, aes(x = X, y = Richness)) +
  
  geom_boxplot(
    aes(fill = SiteDepth),
    alpha = 0.8,
    width = 0.6,
    outlier.shape = NA
  ) +
  
  geom_point(aes(fill = SiteDepth),
             color = "black",
             size = 2.3,
             alpha = 0.9,
             position = position_dodge(width = 0.6)   
  ) +
  
  scale_fill_manual(values = colors_depth_site, name = "Site & depth") +
  
  theme_classic(base_size = 14) +
  labs(
    x = "Site and biological replicate (core)",
    y = "ASV richness"
  ) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.title  = element_text(face = "bold"),
    panel.border  = element_rect(fill = NA)
  )


ggplot(alpha_rep_meta_clean, aes(x = X, y = Shannon)) +
  
  geom_boxplot(
    aes(fill = SiteDepth),
    alpha = 0.8,
    width = 0.6,
    outlier.shape = NA
  ) +
  
  geom_point(aes(fill = SiteDepth),
             color = "black",
             size = 2.3,
             alpha = 0.9,
             position = position_dodge(width = 0.6)   
  ) +
  
  scale_fill_manual(values = colors_depth_site, name = "Site & depth") +
  
  theme_classic(base_size = 14) +
  labs(
    x = "Site and biological replicate (core)",
    y = "ASV Shannon"
  ) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.title  = element_text(face = "bold"),
    panel.border  = element_rect(fill = NA)
  )

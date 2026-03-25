### Mixed Models

library(lme4)
library(reshape)
library(dplyr)
library(readxl)
library(MASS)
library(tidyr)
library(data.table)

reps <- read.xlsx("Clean_Data_COI_Reps_AbRel.xlsx")
colnames(reps)[colnames(reps) == "X1"] <- "ASV"
metadata_reps <- read.xlsx ("metadata_repliques_f.xlsx")
taxo <- read.xlsx("Clean_Data_COI_AVS_taxonomy.xlsx")
colnames(taxo)[colnames(taxo) == "X1"] <- "ASV"

rownames(reps) <- reps$ASV
reps_all <- reps[,c(2:127)]


reps_all_t <- as.data.frame(t(reps_all))
reps_all_t$Sample <- rownames(reps_all_t)

reps_all_t[, 1:(ncol(reps_all_t)-1)] <-
  lapply(reps_all_t[, 1:(ncol(reps_all_t)-1)], as.numeric)

asv_all <- reps_all_t[, !(names(reps_all_t) %in% "Sample")]


alpha_rep_coi <- data.frame(
  Sample   = reps_all_t$Sample,
  Richness = specnumber(asv_all),
  Shannon  = diversity(asv_all, index = "shannon")
)

alpha_rep_coi$Marker <- "COI"

###18S

reps_out_18 <- read.xlsx("Clean_Data_18S_Reps_AbRel.xlsx")
colnames(reps_out_18)[colnames(reps_out_18) == "X1"] <- "ASV"
metadata_reps <- read.xlsx ("metadata_repliques_f.xlsx")
taxo <- read.xlsx("Clean_Data_18S_AVS_taxonomy.xlsx")
colnames(taxo)[colnames(taxo) == "X1"] <- "ASV"

rownames(reps_out_18) <- reps_out_18$ASV
reps_out_18 <- reps_out_18[,c(2:146)]


reps_out_18t <- as.data.frame(t(reps_out_18))
reps_out_18t$Sample <- rownames(reps_out_18t)

reps_out_18t[, 1:(ncol(reps_out_18t)-1)] <-
  lapply(reps_out_18t[, 1:(ncol(reps_out_18t)-1)], as.numeric)

asv_all_18 <- reps_out_18t[, !(names(reps_out_18t) %in% "Sample")]


alpha_rep_18 <- data.frame(
  Sample   = reps_out_18t$Sample,
  Richness = specnumber(asv_all_18),
  Shannon  = diversity(asv_all_18, index = "shannon")
)

alpha_rep_18$Marker <- "18S"


alpha_rep_all<- rbind(alpha_rep_coi, alpha_rep_18)

alpha_rep_all <- merge(alpha_rep_all, metadata_reps, by="Sample")


alpha_rep_all$Core_Rep <- factor(alpha_rep_all$Core_Rep)
alpha_rep_all$Site   <- factor(alpha_rep_all$Site)
alpha_rep_all$Marker   <- factor(alpha_rep_all$Marker)
alpha_rep_all$AgeGroup <- ifelse(alpha_rep_all$Depth %in% c(4, 5), "recent", "old")

alpha_rep_all$AgeGroup <- factor(alpha_rep_all$AgeGroup, levels=c("old","recent"))







### Reps tech

alpha_rep_all_18s <- alpha_rep_all%>%
  filter(Marker=="18S")

alpha_rep_all_coi <- alpha_rep_all%>%
  filter(Marker=="COI")

## 18S
m_full <- glmer.nb(Richness ~ Site + AgeGroup +
                     (1 | Site:Core_Rep) +
                     (1 | Site:Core_Rep:TechRep),
                   data = alpha_rep_all_18s)


VarCorr(m_full)
vc <- as.data.frame(VarCorr(m_full))
vc[, c("grp","vcov","sdcor")]
vc <- as.data.frame(VarCorr(m_full))

total_var <- sum(vc$vcov)

vc$prop <- vc$vcov / total_var

vc

m_full_s <- lmer(Shannon ~ Site + AgeGroup +
                   (1 | Site:Core_Rep) +
                   (1 | Site:Core_Rep:TechRep),
                 data = alpha_rep_all_18s,
                 REML = TRUE)


VarCorr(m_full_s)
vc <- as.data.frame(VarCorr(m_full_s))
vc[, c("grp","vcov","sdcor")]

vc <- as.data.frame(VarCorr(m_full_s))

total_var <- sum(vc$vcov)

vc$prop <- vc$vcov / total_var

vc


##COI

m_full <- glmer.nb(Richness ~ Site + AgeGroup +
                     (1 | Site:Core_Rep) +
                     (1 | Site:Core_Rep:TechRep),
                   data = alpha_rep_all_coi)


vc <- as.data.frame(VarCorr(m_full))
vc[, c("grp","vcov","sdcor")]

vc <- as.data.frame(VarCorr(m_full))

total_var <- sum(vc$vcov)

vc$prop <- vc$vcov / total_var

vc



m_full_s <- lmer(Shannon ~ Site + AgeGroup +
                   (1 | Site:Core_Rep) +
                   (1 | Site:Core_Rep:TechRep),
                 data = alpha_rep_all_coi,
                 REML = TRUE)


vc <- as.data.frame(VarCorr(m_full_s))
vc[, c("grp","vcov","sdcor")]

vc <- as.data.frame(VarCorr(m_full_s))

total_var <- sum(vc$vcov)

vc$prop <- vc$vcov / total_var

vc




######### METAZOA
################## METAZOA
##### COI

reps_out_metazoa <- reps_out %>%
  filter(ASV %in% taxo$ASV[taxo$kingdom == "Metazoa"])

rownames(reps_out_metazoa) <- reps_out_metazoa$ASV
reps_out_metazoa <- reps_out_metazoa[,c(2:127)]


reps_out_metazoa_t <- as.data.frame(t(reps_out_metazoa))
reps_out_metazoa_t$Sample <- rownames(reps_out_metazoa_t)

reps_out_metazoa_t[, 1:(ncol(reps_out_metazoa_t)-1)] <-
  lapply(reps_out_metazoa_t[, 1:(ncol(reps_out_metazoa_t)-1)], as.numeric)

asv_mat_rep_metazoa <- reps_out_metazoa_t[, !(names(reps_out_metazoa_t) %in% "Sample")]


alpha_rep_metazoa <- data.frame(
  Sample   = reps_out_metazoa_t$Sample,
  Richness = specnumber(asv_mat_rep_metazoa),
  Shannon  = diversity(asv_mat_rep_metazoa, index = "shannon")
)

alpha_rep_metazoa$Marker <- "COI"


###18S

reps_out_metazoa_18 <- reps_out_18 %>%
  filter(ASV %in% taxo$ASV[taxo$kingdom == "Metazoa"])

rownames(reps_out_metazoa_18) <- reps_out_metazoa_18$ASV
reps_out_metazoa_18 <- reps_out_metazoa_18[,c(2:146)]

reps_out_metazoa_t_18 <- as.data.frame(t(reps_out_metazoa_18))
reps_out_metazoa_t_18$Sample <- rownames(reps_out_metazoa_t_18)

reps_out_metazoa_t_18[, 1:(ncol(reps_out_metazoa_t_18)-1)] <-
  lapply(reps_out_metazoa_t_18[, 1:(ncol(reps_out_metazoa_t_18)-1)], as.numeric)

asv_mat_rep_metazoa_18 <- reps_out_metazoa_t_18[, !(names(reps_out_metazoa_t_18) %in% "Sample")]


alpha_rep_metazoa_18 <- data.frame(
  Sample   = reps_out_metazoa_t_18$Sample,
  Richness = specnumber(asv_mat_rep_metazoa_18),
  Shannon  = diversity(asv_mat_rep_metazoa_18, index = "shannon")
)

alpha_rep_metazoa_18$Marker <- "18S"


alpha_rep_meta_all_metazoa<- rbind(alpha_rep_metazoa, alpha_rep_metazoa_18)

alpha_rep_all_metazoa <- merge(alpha_rep_meta_all_metazoa, metadata_reps, by="Sample")


alpha_rep_all_metazoa$Core_Rep <- factor(alpha_rep_all_metazoa$Core_Rep)
alpha_rep_all_metazoa$Site   <- factor(alpha_rep_all_metazoa$Site)
alpha_rep_all_metazoa$Marker   <- factor(alpha_rep_all_metazoa$Marker)
alpha_rep_all_metazoa$AgeGroup <- ifelse(alpha_rep_all_metazoa$Depth %in% c(4, 5), "recent", "old")

alpha_rep_all_metazoa$AgeGroup <- factor(alpha_rep_all_metazoa$AgeGroup, levels=c("old","recent"))



#write.xlsx(alpha_rep_all_metazoa, "Diversities_Metazoa_BothGENES.XLSX")


### Rep Techs Metazoa

alpha_rep_all_metazoa_18s <- alpha_rep_all_metazoa%>%
  filter(Marker=="18S")

alpha_rep_all_metazoa_coi <- alpha_rep_all_metazoa%>%
  filter(Marker=="COI")

## 18S
m_full_meta <- glmer.nb(Richness ~ Site + AgeGroup +
                          (1 | Site:AgeGroup:Core_Rep) +
                          (1 | Site:Core_Rep:TechRep),
                        data = alpha_rep_all_metazoa_18s)


VarCorr(m_full_meta)
vc <- as.data.frame(VarCorr(m_full_meta))
vc[, c("grp","vcov","sdcor")]

vc <- as.data.frame(VarCorr(m_full_meta))

total_var <- sum(vc$vcov)

vc$prop <- vc$vcov / total_var

vc

m_full_meta_s <- lmer(Shannon ~ Site + AgeGroup +
                        (1 | Site:AgeGroup:Core_Rep) +
                        (1 | Site:Core_Rep:TechRep),
                      data = alpha_rep_all_metazoa_18s,
                      REML = TRUE)


VarCorr(m_full_meta_s)
vc <- as.data.frame(VarCorr(m_full_meta_s))
vc[, c("grp","vcov","sdcor")]
vc <- as.data.frame(VarCorr(m_full_meta_s))

total_var <- sum(vc$vcov)

vc$prop <- vc$vcov / total_var

vc


##COI

m_full_meta <- glmer.nb(Richness ~ Site + AgeGroup +
                          (1 | Site:Core_Rep) +
                          (1 | Site:Core_Rep:TechRep),
                        data = alpha_rep_all_metazoa_coi)


VarCorr(m_full_meta)
vc <- as.data.frame(VarCorr(m_full_meta))
vc[, c("grp","vcov","sdcor")]
vc <- as.data.frame(VarCorr(m_full_meta))

total_var <- sum(vc$vcov)

vc$prop <- vc$vcov / total_var

vc

m_full_meta_s <- lmer(Shannon ~ Site + AgeGroup +
                        (1 | Site:Core_Rep) +
                        (1 | Site:Core_Rep:TechRep),
                      data = alpha_rep_all_metazoa_coi,
                      REML = TRUE)


VarCorr(m_full_meta_s)
vc <- as.data.frame(VarCorr(m_full_meta_s))
vc[, c("grp","vcov","sdcor")]
vc <- as.data.frame(VarCorr(m_full_meta_s))

total_var <- sum(vc$vcov)

vc$prop <- vc$vcov / total_var

vc









###################### Rep Bio

### Differences between Biological Replicates
### Eukariota
m_core1 <- lm(
  Shannon ~ Marker * CoreID * AgeGroup * Core_Rep,
  data = alpha_rep_all
)

summary(m_core1)
emmeans(m_core1, pairwise ~ Core_Rep | Marker + CoreID + AgeGroup)


#### Richness
m_core2 <- glm.nb(
  Richness ~ Marker * CoreID * AgeGroup * Core_Rep,
  data = alpha_rep_all
)

emmeans(m_core2, pairwise ~ Core_Rep | Marker + CoreID + AgeGroup)



### Metazoa

m_core1 <- lm(
  Shannon ~ Marker * CoreID * AgeGroup * Core_Rep,
  data = alpha_rep_meta_all
)

summary(m_core1)
emmeans(m_core1, pairwise ~ Core_Rep | Marker + CoreID + AgeGroup)


#### Richness
m_core2 <- glm.nb(
  Richness ~ Marker * CoreID * AgeGroup * Core_Rep,
  data = alpha_rep_meta_all
)
library(emmeans)

emmeans(m_core2, pairwise ~ Core_Rep | Marker + CoreID + AgeGroup)


#### Differences in Sites, Age Group and Marker within Biological Replicates
### Shannon
m_core4 <- lmer(
  Shannon ~ Marker * CoreID * AgeGroup +
    (1 | Marker:CoreID:Core_Rep),
  data = alpha_rep_all
)

summary(m_core4)
anova(m_core4)
VarCorr(m_core4)

emmeans(m_core4, pairwise ~ CoreID | AgeGroup + Marker)


## Richness

m_core3 <- glmer.nb(
  Richness ~ Marker * CoreID * AgeGroup +
    (1 | Marker:CoreID:Core_Rep),
  data = alpha_rep_all
)

summary(m_core3)
anova(m_core3)
VarCorr(m_core3)
emmeans(m_core3, pairwise ~ Site | AgeGroup + Marker)



#####METAZOA

### Shannon
m_core4 <- lmer(
  Shannon ~ Marker * CoreID * AgeGroup +
    (1 | Marker:CoreID:Core_Rep),
  data = alpha_rep_meta_all
)

summary(m_core4)
anova(m_core4)
VarCorr(m_core4)

emmeans(m_core4, pairwise ~ CoreID | AgeGroup + Marker)


## Richness

m_core3 <- glmer.nb(
  Richness ~ Marker * CoreID * AgeGroup +
    (1 | Marker:CoreID:Core_Rep),
  data = alpha_rep_meta_all
)

summary(m_core3)
anova(m_core3)
VarCorr(m_core3)
emmeans(m_core3, pairwise ~ Site | AgeGroup + Marker)




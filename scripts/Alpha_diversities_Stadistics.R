############# Alpha diversities Stadistics

library(lme4)
library(reshape)
library(dplyr)
library(readxl)
library(MASS)
library(tidyr)
library(data.table)
library(performance)
library(emmeans)


### COI

reps <- read.xlsx("Clean_Data_COI_Reps_AbRel.xlsx")
colnames(reps)[colnames(reps) == "X1"] <- "ASV"
metadata_reps <- read.xlsx ("metadata_repliques_f.xlsx")

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

reps_18 <- read.xlsx("Clean_Data_18S_Reps_AbRel.xlsx")
colnames(reps_18)[colnames(reps_18) == "X1"] <- "ASV"

rownames(reps_18) <- reps_18$ASV
reps_18_all <- reps_18[,c(2:142)]

reps_18_all_t <- as.data.frame(t(reps_18_all))
reps_18_all_t$Sample <- rownames(reps_18_all_t)

reps_18_all_t[, 1:(ncol(reps_18_all_t)-1)] <-
  lapply(reps_18_all_t[, 1:(ncol(reps_18_all_t)-1)], as.numeric)

asv_all_18 <- reps_18_all_t[, !(names(reps_18_all_t) %in% "Sample")]

alpha_rep_18 <- data.frame(
  Sample   = reps_18_all_t$Sample,
  Richness = specnumber(asv_all_18),
  Shannon  = diversity(asv_all_18, index = "shannon")
)

alpha_rep_18$Marker <- "18S"


##### Pool Markers and merge with metadata

alpha_all<- rbind(alpha_rep_coi, alpha_rep_18)

alpha_all <- merge(alpha_all, metadata_reps, by="Sample")


alpha_all$Core_Rep <- factor(alpha_all$Core_Rep)
alpha_all$TechRep <- factor(alpha_all$TechRep)
alpha_all$Site   <- factor(alpha_all$Site)
alpha_all$Marker   <- factor(alpha_all$Marker)
alpha_all$AgeGroup <- ifelse(alpha_all$Depth %in% c(4, 5), "recent", "old")
alpha_all$AgeGroup <- factor(alpha_all$AgeGroup, levels=c("old","recent"))
alpha_all$Depth   <- factor(alpha_all$Depth)






### Stadistics Technical Replicates - Mixed model
## Response variable - Diversity
## Fixed Effects - Site and AgeGroup and TechRep
## Random Effects - (1 | Site:Core_Rep:AgeGroup)


alpha_all_18s <- alpha_all%>%
  filter(Marker=="18S")

alpha_all_coi <- alpha_all%>%
  filter(Marker=="COI")

## 18S


m_rich <- glmer.nb(
  Richness ~ Site + AgeGroup + TechRep +
    (1 | Site:Core_Rep:AgeGroup),
  data = alpha_all_18s
)

m_null <- glmer.nb(
  Richness ~ Site + AgeGroup +
    (1 | Site:Core_Rep:AgeGroup),
  data = alpha_all_18s
)

anova(m_null, m_rich)


summary(m_rich)
VarCorr(m_rich)
vc <- as.data.frame(VarCorr(m_rich))
vc[, c("grp","vcov","sdcor")]

vc <- as.data.frame(VarCorr(m_rich))

total_var <- sum(vc$vcov)

vc$prop <- vc$vcov / total_var

vc
icc(m_rich)


m_shannon <- lmer(
  Shannon ~ Site + AgeGroup + TechRep +
    (1 | Site:Core_Rep:AgeGroup),
  data = alpha_all_18s,
  REML = FALSE
)

m_null_s <- lmer(
  Shannon ~ Site + AgeGroup +
    (1 | Site:Core_Rep:AgeGroup),
  data = alpha_all_18s,
  REML = FALSE
)

anova(m_null_s, m_shannon)

summary(m_shannon)
VarCorr(m_shannon)
vc <- as.data.frame(VarCorr(m_shannon))
vc[, c("grp","vcov","sdcor")]

vc <- as.data.frame(VarCorr(m_shannon))

total_var <- sum(vc$vcov)

vc$prop <- vc$vcov / total_var

vc
icc(m_shannon)


##COI


m_rich <- glmer.nb(
  Richness ~ Site + AgeGroup + TechRep +
    (1 | Site:Core_Rep:AgeGroup),
  data = alpha_all_coi
)

m_null <- glmer.nb(
  Richness ~ Site + AgeGroup +
    (1 | Site:Core_Rep:AgeGroup),
  data = alpha_all_coi
)

anova(m_null, m_rich)


summary(m_rich)
VarCorr(m_rich)
vc <- as.data.frame(VarCorr(m_rich))
vc[, c("grp","vcov","sdcor")]

vc <- as.data.frame(VarCorr(m_rich))

total_var <- sum(vc$vcov)

vc$prop <- vc$vcov / total_var

vc
icc(m_rich)


m_shannon <- lmer(
  Shannon ~ Site + AgeGroup + TechRep +
    (1 | Site:Core_Rep:AgeGroup),
  data = alpha_all_coi,
  REML = FALSE
)

m_null_s <- lmer(
  Shannon ~ Site + AgeGroup +
    (1 | Site:Core_Rep:AgeGroup),
  data = alpha_all_coi,
  REML = FALSE
)

anova(m_null_s, m_shannon)

summary(m_shannon)
VarCorr(m_shannon)
vc <- as.data.frame(VarCorr(m_shannon))
vc[, c("grp","vcov","sdcor")]

vc <- as.data.frame(VarCorr(m_shannon))

total_var <- sum(vc$vcov)

vc$prop <- vc$vcov / total_var

vc
icc(m_shannon)



######### METAZOA
################## METAZOA
##### COI

tax_COI <- read.xlsx("Clean_Data_COI_AVS_taxonomy.xlsx")

metazoa_reps_coi <- reps %>%
  left_join(tax_COI, by = "ASV") %>%
  filter(kingdom == "Metazoa")


rownames(metazoa_reps_coi) <- metazoa_reps_coi$ASV
metazoa_reps_coi <- metazoa_reps_coi[,c(2:127)]


metazoa_reps_coi_t <- as.data.frame(t(metazoa_reps_coi))
metazoa_reps_coi_t$Sample <- rownames(metazoa_reps_coi_t)

metazoa_reps_coi_t[, 1:(ncol(metazoa_reps_coi_t)-1)] <-
  lapply(metazoa_reps_coi_t[, 1:(ncol(metazoa_reps_coi_t)-1)], as.numeric)

asv_mat_rep_metazoa <- metazoa_reps_coi_t[, !(names(metazoa_reps_coi_t) %in% "Sample")]



alpha_metazoa <- data.frame(
  Sample   = metazoa_reps_coi_t$Sample,
  Richness = specnumber(asv_mat_rep_metazoa),
  Shannon  = diversity(asv_mat_rep_metazoa, index = "shannon")
)

alpha_metazoa$Marker <- "COI"

###18S
tax_18 <- read.xlsx("Clean_Data_18S_AVS_taxonomy.xlsx")

metazoa_reps_18 <- reps_18 %>%
  left_join(tax_18, by = "ASV") %>%
  filter(kingdom == "Metazoa")


rownames(metazoa_reps_18) <- metazoa_reps_18$ASV
metazoa_reps_18 <- metazoa_reps_18[,c(2:142)]


metazoa_reps_18_t <- as.data.frame(t(metazoa_reps_18))
metazoa_reps_18_t$Sample <- rownames(metazoa_reps_18_t)

metazoa_reps_18_t[, 1:(ncol(metazoa_reps_18_t)-1)] <-
  lapply(metazoa_reps_18_t[, 1:(ncol(metazoa_reps_18_t)-1)], as.numeric)

asv_mat_rep_metazoa_18 <- metazoa_reps_18_t[, !(names(metazoa_reps_18_t) %in% "Sample")]


alpha_metazoa_18 <- data.frame(
  Sample   = metazoa_reps_18_t$Sample,
  Richness = specnumber(asv_mat_rep_metazoa_18),
  Shannon  = diversity(asv_mat_rep_metazoa_18, index = "shannon")
)

alpha_metazoa_18$Marker <- "18S"

alpha_meta_all_metazoa<- rbind(alpha_metazoa, alpha_metazoa_18)

alpha_all_metazoa <- merge(alpha_meta_all_metazoa, metadata_reps, by="Sample")


alpha_all_metazoa$Core_Rep <- factor(alpha_all_metazoa$Core_Rep)
alpha_all_metazoa$TechRep <- factor(alpha_all_metazoa$TechRep)
alpha_all_metazoa$Site   <- factor(alpha_all_metazoa$Site)
alpha_all_metazoa$Marker   <- factor(alpha_all_metazoa$Marker)
alpha_all_metazoa$AgeGroup <- ifelse(alpha_all_metazoa$Depth %in% c(4, 5), "recent", "old")
alpha_all_metazoa$AgeGroup <- factor(alpha_all_metazoa$AgeGroup, levels=c("old","recent"))
alpha_all_metazoa$Depth   <- factor(alpha_all_metazoa$Depth)


#write.xlsx(alpha_rep_all_metazoa, "Diversities_Metazoa_BothGENES.XLSX")


### Rep Techs

alpha_rep_all_metazoa_18s <- alpha_rep_all_metazoa%>%
  filter(Marker=="18S")

alpha_rep_all_metazoa_coi <- alpha_rep_all_metazoa%>%
  filter(Marker=="COI")

## 18S

m_rich <- glmer.nb(
  Richness ~ Site + AgeGroup + TechRep +
    (1 | Site:Core_Rep:AgeGroup),
  data = alpha_rep_all_metazoa_18s
)

m_null <- glmer.nb(
  Richness ~ Site + AgeGroup +
    (1 | Site:Core_Rep:AgeGroup),
  data = alpha_rep_all_metazoa_18s
)

anova(m_null, m_rich)


summary(m_rich)
VarCorr(m_rich)
vc <- as.data.frame(VarCorr(m_rich))
vc[, c("grp","vcov","sdcor")]

vc <- as.data.frame(VarCorr(m_rich))

total_var <- sum(vc$vcov)

vc$prop <- vc$vcov / total_var

vc
icc(m_rich)


m_shannon <- lmer(
  Shannon ~ Site + AgeGroup + TechRep +
    (1 | Site:Core_Rep:AgeGroup),
  data = alpha_rep_all_metazoa_18s,
  REML = FALSE
)

m_null_s <- lmer(
  Shannon ~ Site + AgeGroup +
    (1 | Site:Core_Rep:AgeGroup),
  data = alpha_rep_all_metazoa_18s,
  REML = FALSE
)

anova(m_null_s, m_shannon)

summary(m_shannon)
VarCorr(m_shannon)
vc <- as.data.frame(VarCorr(m_shannon))
vc[, c("grp","vcov","sdcor")]

vc <- as.data.frame(VarCorr(m_shannon))

total_var <- sum(vc$vcov)

vc$prop <- vc$vcov / total_var

vc
icc(m_shannon)


##COI


m_rich <- glmer.nb(
  Richness ~ Site + AgeGroup + TechRep +
    (1 | Site:Core_Rep:AgeGroup),
  data = alpha_rep_all_metazoa_coi
)

m_null <- glmer.nb(
  Richness ~ Site + AgeGroup +
    (1 | Site:Core_Rep:AgeGroup),
  data = alpha_rep_all_metazoa_coi
)

anova(m_null, m_rich)


summary(m_rich)
VarCorr(m_rich)
vc <- as.data.frame(VarCorr(m_rich))
vc[, c("grp","vcov","sdcor")]

vc <- as.data.frame(VarCorr(m_rich))

total_var <- sum(vc$vcov)

vc$prop <- vc$vcov / total_var

vc
icc(m_rich)


m_shannon <- lmer(
  Shannon ~ Site + AgeGroup + TechRep +
    (1 | Site:Core_Rep:AgeGroup),
  data = alpha_rep_all_metazoa_coi,
  REML = FALSE
)

m_null_s <- lmer(
  Shannon ~ Site + AgeGroup +
    (1 | Site:Core_Rep:AgeGroup),
  data = alpha_rep_all_metazoa_coi,
  REML = FALSE
)

anova(m_null_s, m_shannon)

summary(m_shannon)
VarCorr(m_shannon)
vc <- as.data.frame(VarCorr(m_shannon))
vc[, c("grp","vcov","sdcor")]

vc <- as.data.frame(VarCorr(m_shannon))

total_var <- sum(vc$vcov)

vc$prop <- vc$vcov / total_var

vc
icc(m_shannon)







###################### Rep Bio
### Differences across Bio Replicates??

### Eukaryota Community
##Shannon


m_core1 <- lm(
  Shannon ~ Marker * Site * AgeGroup + 
    Marker:Site:Core_Rep:AgeGroup,
  data = alpha_all
)

summary(m_core1)
emmeans(m_core1, pairwise ~ Core_Rep | Marker + Site + AgeGroup)


##Richness

m_core2 <- glm.nb(
  Richness ~ Marker * Site * AgeGroup + 
    Marker:Site:Core_Rep:AgeGroup,
  data = alpha_all
)

emmeans(m_core2, pairwise ~ Core_Rep | Marker + Site + AgeGroup)




### Metazoan Community
##Shannon

m_core3 <- lm(
  Shannon ~ Marker * Site * AgeGroup + 
    Marker:Site:Core_Rep:AgeGroup,
  data = alpha_all_metazoa
)

summary(m_core3)
emmeans(m_core3, pairwise ~ Core_Rep | Marker + Site + AgeGroup)


##Richness

m_core4 <- glm.nb(
  Richness ~ Marker * Site * AgeGroup + 
    Marker:Site:Core_Rep:AgeGroup,
  data = alpha_all_metazoa
)

emmeans(m_core4, pairwise ~ Core_Rep | Marker + Site + AgeGroup)



### Differences Across Sites? 

## Eukaryota Community
##Shannon

m_site1 <- lmer(
  Shannon ~ Marker * Site * AgeGroup +
    (1 | Site:Core_Rep),
  data = alpha_all
)

emmeans(m_site1, pairwise ~ Site | Marker + AgeGroup, adjust = "tukey")

## Richness
m_site2 <- glmer.nb(
  Richness ~ Marker * Site * AgeGroup +
    (1 | Site:Core_Rep),
  data = alpha_all
)
emmeans(m_site2, pairwise ~ Site | Marker + AgeGroup, adjust = "tukey")


## Metazoan Community
##Shannon

m_site3 <- lmer(
  Shannon ~ Marker * Site * AgeGroup +
    (1 | Site:Core_Rep),
  data = alpha_all_metazoa
)

emmeans(m_site3, pairwise ~ Site | Marker + AgeGroup, adjust = "tukey")

## Richness
m_site2 <- glmer.nb(
  Richness ~ Marker * Site * AgeGroup +
    (1 | Site:Core_Rep),
  data = alpha_all_metazoa
)
emmeans(m_site2, pairwise ~ Site | Marker + AgeGroup, adjust = "tukey")




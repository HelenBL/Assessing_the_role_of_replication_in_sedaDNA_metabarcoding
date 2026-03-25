## HMSC 

library(openxlsx)
library(dplyr)
library(tidyr)
library(tibble)
library(Hmsc)
library(coda)



metadata <- read.xlsx("metadata_f.xlsx")

## COI
arel_mitjCOI_kp <- read.xlsx("Clean_Data_COI_AbRel.xlsx")
rownames(arel_mitjCOI_kp) <- arel_mitjCOI_kp$X1
arel_mitjCOI_kp <- arel_mitjCOI_kp[, 2:21]

arel_mitjCOI_kp_long <- pivot_longer(
  arel_mitjCOI_kp,
  cols = matches("^CCO|^CLR|^STO"),
  names_to = "Sample",
  values_to = "ReadsCOI"
)

abdCOI_meta <- left_join(arel_mitjCOI_kp_long, metadata, by = "Sample")
abdCOI_meta$Core_Rep <- as.factor(abdCOI_meta$Core_Rep)
abdCOI_meta$Depth <- as.factor(abdCOI_meta$Depth)
abdCOI_meta$Site <- as.factor(abdCOI_meta$Site)

## 18S
arel_mitj18S_kp <- read.xlsx("Clean_Data_18S_AbRel.xlsx")
rownames(arel_mitj18S_kp) <- arel_mitj18S_kp$X1
arel_mitj18S_kp <- arel_mitj18S_kp[, 2:21]

arel_mitj18S_kp_long <- pivot_longer(
  arel_mitj18S_kp,
  cols = matches("^CCO|^CLR|^STO"),
  names_to = "Sample",
  values_to = "Reads18S"
)

abd18S_meta <- left_join(arel_mitj18S_kp_long, metadata, by = "Sample")
abd18S_meta$Core_Rep <- as.factor(abd18S_meta$Core_Rep)
abd18S_meta$Depth <- as.factor(abd18S_meta$Depth)
abd18S_meta$Site <- as.factor(abd18S_meta$Site)

## Extra Functions
add_agegroup <- function(df) {
  if ("AgeGroup" %in% names(df)) {
    df$AgeGroup <- factor(df$AgeGroup, levels = c("old", "recent"))
    return(df)
  }
  
  if ("Depth" %in% names(df)) {
    df <- df %>%
      mutate(
        Depth_chr = as.character(Depth),
        Depth_num = suppressWarnings(as.numeric(Depth_chr)),
        AgeGroup = ifelse(Depth_num <= 5, "recent", "old"),
        AgeGroup = factor(AgeGroup, levels = c("old", "recent"))
      ) %>%
      select(-Depth_chr, -Depth_num)
    return(df)
  }
  
  if ("Year" %in% names(df)) {
    recent_years <- c("2014", "2017", "1974")
    df <- df %>%
      mutate(
        AgeGroup = ifelse(as.character(Year) %in% recent_years, "recent", "old"),
        AgeGroup = factor(AgeGroup, levels = c("old", "recent"))
      )
    return(df)
  }
  
  stop("No AgeGroup, Depth o Year to define age groups.")
}


prepare_hmsc_data <- function(df,
                              value_col,
                              tax_col = "phylum",
                              sample_col = "Sample",
                              site_col = "Site",
                              core_col = "Core_Rep",
                              keep_only_metazoa = FALSE,
                              kingdom_col = "kingdom") {
  
  df <- add_agegroup(df)
  
  if (keep_only_metazoa) {
    df <- df %>% filter(.data[[kingdom_col]] == "Metazoa")
  }
  
  df2 <- df %>%
    select(
      sample   = all_of(sample_col),
      site     = all_of(site_col),
      core_rep = all_of(core_col),
      AgeGroup,
      taxon    = all_of(tax_col),
      value    = all_of(value_col)
    ) %>%
    mutate(
      sample   = as.character(sample),
      site     = factor(site),
      core_rep = factor(core_rep),
      core_id  = factor(interaction(site, core_rep, drop = TRUE)),
      taxon    = as.character(taxon)
    ) %>%
    filter(!is.na(taxon), taxon != "")
  
  comm <- df2 %>%
    group_by(sample, taxon) %>%
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = taxon, values_from = value, values_fill = 0)
  
  meta <- df2 %>%
    distinct(sample, site, core_rep, core_id, AgeGroup) %>%
    arrange(sample)
  
  comm <- comm %>%
    slice(match(meta$sample, sample))
  
  stopifnot(all(comm$sample == meta$sample))
  
  Y <- comm %>%
    select(-sample) %>%
    as.matrix()
  
  Y <- ifelse(Y > 0, 1, 0)
  storage.mode(Y) <- "numeric"
  
  keep_taxa <- apply(Y, 2, function(x) length(unique(x)) > 1)
  Y <- Y[, keep_taxa, drop = FALSE]
  
  XData <- meta %>%
    select(site, AgeGroup) %>%
    mutate(
      site = factor(site),
      AgeGroup = factor(AgeGroup, levels = c("old", "recent"))
    )
  
  studyDesign <- data.frame(
    sample = factor(meta$sample),
    site   = factor(meta$site),
    core   = factor(meta$core_id)
  )
  
  return(list(
    Y = Y,
    XData = XData,
    studyDesign = studyDesign,
    meta = meta
  ))
}


fit_hmsc_pa <- function(prepped,
                        nChains = 2,
                        samples = 1000,
                        transient = 500,
                        thin = 10,
                        nParallel = 2,
                        verbose = 100) {
  
  rL.site <- HmscRandomLevel(units = prepped$studyDesign$site)
  rL.core <- HmscRandomLevel(units = prepped$studyDesign$core)
  
  m <- Hmsc(
    Y = prepped$Y,
    XData = prepped$XData,
    XFormula = ~ site * AgeGroup,
    studyDesign = prepped$studyDesign,
    ranLevels = list(site = rL.site, core = rL.core),
    distr = "probit"
  )
  
  m <- sampleMcmc(
    m,
    thin = thin,
    samples = samples,
    transient = transient,
    nChains = nChains,
    nParallel = nParallel,
    verbose = verbose
  )
  
  return(m)
}


summarise_beta <- function(model) {
  getPostEstimate(model, parName = "Beta")
}

variance_partition <- function(model) {
  cn <- colnames(model$X)
  print(cn)  
  
  group <- rep(NA_integer_, length(cn))
  
  group[cn == "(Intercept)"] <- 1
  
  group[grepl("^site", cn)] <- 2
  group[grepl("^AgeGroup", cn)] <- 3
  
  group[grepl("site:AgeGroup|AgeGroup:site", cn)] <- 4
  
  if (any(is.na(group))) {
    cat("Columns without assigned group:\n")
    print(cn[is.na(group)])
  }
  
  present_groups <- sort(unique(group[!is.na(group)]))
  
  group_reindexed <- match(group, present_groups)
  
  all_names <- c("Intercept", "Site", "AgeGroup", "Site:AgeGroup")
  groupnames_reindexed <- all_names[present_groups]
  
  keep <- !is.na(group_reindexed)
  
  VP <- computeVariancePartitioning(
    model,
    group = group_reindexed[keep],
    groupnames = groupnames_reindexed
  )
  
  return(VP)
}

check_hmsc <- function(model) {
  mpost <- convertToCodaObject(model)
  print(summary(mpost$Beta))
  invisible(mpost)
}


plotVariancePartitioning_custom <- function(
    hM, VP,
    cols = NULL,
    main = "Variance Partitioning",
    cex.names = 0.7,
    cex.legend = 0.9,
    ...) {
  
  ng <- dim(VP$vals)[1]
  
  if (is.null(cols)) {
    cols <- c(
      "#1b9e77",  
      "#d95f02",  
      "#7570b3",  
      "#e7298a",  
      "#66a61e",  
      "#e6ab02",  
      "#a6761d",  
      "#666666"   
    )
    
    if (length(cols) < ng) {
      cols <- grDevices::rainbow(ng)
    } else {
      cols <- cols[1:ng]
    }
  }
  
  leg <- VP$groupnames
  
  for (r in seq_len(hM$nr)) {
    leg <- c(leg, paste("Random: ", hM$rLNames[r], sep = ""))
  }
  
  means <- round(100 * rowMeans(VP$vals), 1)
  for (i in seq_len(min(length(means), length(leg)))) {
    leg[i] <- paste(leg[i], " (mean = ", toString(means[i]), ")", sep = "")
  }
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mar = c(10, 4, 4, 2) + 0.1)
  
  barplot(
    VP$vals,
    main = main,
    xlab = "Phylum",
    ylab = "Variance proportion",
    names.arg = colnames(VP$vals),  
    las = 2,                        
    cex.names = cex.names,          
    legend.text = leg,
    args.legend = list(
      x = "topright",
      cex = cex.legend,
      bty = "n"
    ),
    col = cols,
    ...
  )
}

######## DATASETS

# Eukaryota COI
prep_euk_coi <- prepare_hmsc_data(
  df = abdCOI_meta,
  value_col = "ReadsCOI",
  tax_col = "phylum",
  sample_col = "Sample",
  site_col = "Site",
  core_col = "Core_Rep",
  keep_only_metazoa = FALSE
)

# Eukaryota 18S
prep_euk_18s <- prepare_hmsc_data(
  df = abd18S_meta,
  value_col = "Reads18S",
  tax_col = "phylum",
  sample_col = "Sample",
  site_col = "Site",
  core_col = "Core_Rep",
  keep_only_metazoa = FALSE
)

# Metazoa COI
prep_met_coi <- prepare_hmsc_data(
  df = abdCOI_meta,
  value_col = "ReadsCOI",
  tax_col = "phylum",
  sample_col = "Sample",
  site_col = "Site",
  core_col = "Core_Rep",
  keep_only_metazoa = TRUE
)

# Metazoa 18S
prep_met_18s <- prepare_hmsc_data(
  df = abd18S_meta,
  value_col = "Reads18S",
  tax_col = "phylum",
  sample_col = "Sample",
  site_col = "Site",
  core_col = "Core_Rep",
  keep_only_metazoa = TRUE
)

######### MODELS

set.seed(123)

m_euk_coi <- fit_hmsc_pa(prep_euk_coi)
m_euk_18s <- fit_hmsc_pa(prep_euk_18s)
m_met_coi <- fit_hmsc_pa(prep_met_coi)
m_met_18s <- fit_hmsc_pa(prep_met_18s)

###########OUTPUTS

# Fixed Effects
beta_euk_coi <- summarise_beta(m_euk_coi)
beta_euk_18s <- summarise_beta(m_euk_18s)
beta_met_coi <- summarise_beta(m_met_coi)
beta_met_18s <- summarise_beta(m_met_18s)

beta_euk_coi$mean
beta_euk_18s$mean
beta_met_coi$mean
beta_met_18s$mean
beta_euk_coi
beta_euk_18s
beta_met_coi
beta_met_18s

# Basic Diagnostic
check_hmsc(m_euk_coi)
check_hmsc(m_euk_18s)
check_hmsc(m_met_coi)
check_hmsc(m_met_18s)

# Variance partition
VP_euk_coi <- variance_partition(m_euk_coi)
VP_euk_18s <- variance_partition(m_euk_18s)
VP_met_coi <- variance_partition(m_met_coi)
VP_met_18s <- variance_partition(m_met_18s)
VP_euk_coi$vals
VP_euk_18s$vals
VP_met_coi$vals
VP_met_18s$vals

plotVariancePartitioning_custom(m_euk_coi, VP = VP_euk_coi)
plotVariancePartitioning_custom(m_euk_18s, VP = VP_euk_18s)
plotVariancePartitioning_custom(m_met_coi, VP = VP_met_coi)
plotVariancePartitioning_custom(m_met_18s, VP = VP_met_18s)

### Model adjustment
pred_euk_coi <- computePredictedValues(m_euk_coi)
fit_euk_coi <- evaluateModelFit(hM = m_euk_coi, predY = pred_euk_coi)

pred_euk_18s <- computePredictedValues(m_euk_18s)
fit_euk_18s <- evaluateModelFit(hM = m_euk_18s, predY = pred_euk_18s)

pred_met_coi <- computePredictedValues(m_met_coi)
fit_met_coi <- evaluateModelFit(hM = m_met_coi, predY = pred_met_coi)

pred_met_18s <- computePredictedValues(m_met_18s)
fit_met_18s <- evaluateModelFit(hM = m_met_18s, predY = pred_met_18s)

fit_euk_coi
fit_euk_18s
fit_met_coi
fit_met_18s



#### Final Table
vp_to_long <- function(VP, dataset_name) {
  
  df <- as.data.frame(t(VP$vals))
  df$Phylum <- rownames(df)
  
  df_long <- df %>%
    pivot_longer(
      cols = -Phylum,
      names_to = "Effect",
      values_to = "Value"
    ) %>%
    mutate(
      Dataset = dataset_name,
      Metric = "VariancePartitioning"
    )
  
  return(df_long)
}

beta_to_long <- function(beta_obj, dataset_name) {
  
  df <- as.data.frame(beta_obj$mean)
  df$Coefficient = rownames(df)
  
  df_long <- df %>%
    pivot_longer(
      cols = -Coefficient,
      names_to = "Phylum",
      values_to = "Value"
    ) %>%
    mutate(
      Dataset = dataset_name,
      Metric = "Beta",
      Effect = Coefficient
    ) %>%
    select(Dataset, Phylum, Metric, Effect, Value)
  
  return(df_long)
}

fit_to_long <- function(fit_obj, phyla_names, dataset_name) {
  
  df <- data.frame(
    Phylum = phyla_names,
    AUC = fit_obj$AUC,
    TjurR2 = fit_obj$TjurR2,
    RMSE = fit_obj$RMSE
  )
  
  df_long <- df %>%
    pivot_longer(
      cols = -Phylum,
      names_to = "Effect",
      values_to = "Value"
    ) %>%
    mutate(
      Dataset = dataset_name,
      Metric = "ModelFit"
    )
  
  return(df_long)
}

table_all <- bind_rows(
  vp_to_long(VP_euk_coi, "Euk_COI"),
  vp_to_long(VP_euk_18s, "Euk_18S"),
  vp_to_long(VP_met_coi, "Met_COI"),
  vp_to_long(VP_met_18s, "Met_18S"),
  
  beta_to_long(beta_euk_coi, "Euk_COI"),
  beta_to_long(beta_euk_18s, "Euk_18S"),
  beta_to_long(beta_met_coi, "Met_COI"),
  beta_to_long(beta_met_18s, "Met_18S"),
  
  fit_to_long(fit_euk_coi, colnames(VP_euk_coi$vals), "Euk_COI"),
  fit_to_long(fit_euk_18s, colnames(VP_euk_18s$vals), "Euk_18S"),
  fit_to_long(fit_met_coi, colnames(VP_met_coi$vals), "Met_COI"),
  fit_to_long(fit_met_18s, colnames(VP_met_18s$vals), "Met_18S")
)

write.xlsx(table_all, "HMSC_full_results_long.xlsx")


vp_summary <- function(VP, dataset_name) {
  means <- rowMeans(VP$vals)
  
  data.frame(
    Dataset = dataset_name,
    Site = means["Site"],
    AgeGroup = means["AgeGroup"],
    Random_site = means["Random: site"],
    Random_core = means["Random: core"]
  )
}

table_main <- bind_rows(
  vp_summary(VP_euk_coi, "Euk COI"),
  vp_summary(VP_euk_18s, "Euk 18S"),
  vp_summary(VP_met_coi, "Met COI"),
  vp_summary(VP_met_18s, "Met 18S")
)

table_main

table_long <- table_main %>%
  pivot_longer(
    cols = -Dataset,
    names_to = "Effect",
    values_to = "Value"
  )
table_long

ggplot(table_long, aes(x = Dataset, y = Value, fill = Effect)) +
  geom_col() +
  scale_fill_manual(values = c(
    "Site" = "#1b9e77",
    "AgeGroup" = "#d95f02",
    "Random_site" = "#7570b3",
    "Random_core" = "#e7298a"
  )) +
  labs(
    y = "Explained variance",
    x = "",
    fill = "Effect"
  ) +
  theme_classic(base_size = 14)



model_fit_summary_short <- function(fit_obj, dataset_name) {
  data.frame(
    Dataset = dataset_name,
    Mean_AUC = mean(fit_obj$AUC, na.rm = TRUE),
    Mean_TjurR2 = mean(fit_obj$TjurR2, na.rm = TRUE),
    Mean_RMSE = mean(fit_obj$RMSE, na.rm = TRUE)
  )
}

table_fit_short <- bind_rows(
  model_fit_summary_short(fit_euk_coi, "Euk COI"),
  model_fit_summary_short(fit_euk_18s, "Euk 18S"),
  model_fit_summary_short(fit_met_coi, "Met COI"),
  model_fit_summary_short(fit_met_18s, "Met 18S")
)

table_fit_short





model_fit_long <- function(fit_obj, phyla_names, dataset_name) {
  data.frame(
    Dataset = dataset_name,
    Phylum = phyla_names,
    AUC = fit_obj$AUC,
    TjurR2 = fit_obj$TjurR2,
    RMSE = fit_obj$RMSE
  )
}

table_fit_full <- bind_rows(
  model_fit_long(fit_euk_coi, colnames(prep_euk_coi$Y), "Euk COI"),
  model_fit_long(fit_euk_18s, colnames(prep_euk_18s$Y), "Euk 18S"),
  model_fit_long(fit_met_coi, colnames(prep_met_coi$Y), "Met COI"),
  model_fit_long(fit_met_18s, colnames(prep_met_18s$Y), "Met 18S")
)

table_fit_full

write.xlsx(
  list(
    "HMSC_fit_summary" = table_main,
    "HMSC_fit_full" = table_fit_full
  ),
  file = "HMSC_model_fit_tables.xlsx"
)

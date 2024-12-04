library(tidyverse)
library(data.table)

library(meta)
library(parallel)

# Rscript mets_run_metas_w_correction.R ../metadata_20241017/metadata_CRC_staging_20241017.tsv inp_dats.txt out_metas/ comparisons_metas.txt

# This script assumes that the minimum number of samples in at least one of the two categories to compare is 10
# and that the features is present in both cases and controls

# INP_METAD = "../data_March23/metadata_formatted_Mar23.tsv"
# INP_TABLE = "out_metas_divs/divs_formatted.tsv"
# OUT_DIR = "out_metas_divs/"
# TOT_COMPARISONS = "comparisons_metas.txt"

args = commandArgs(trailingOnly=TRUE)
# 
INP_METAD = args[1]
INP_DATS = args[2]
OR_OUT_DIR = args[3]
TOT_COMPARISONS = args[4]
MAX_POOL = 20

# 
# INP_METAD = "../metadata_20241017/metadata_CRC_staging_20241017.tsv"
# INP_DATS = "../Figure1/inp_dats.txt"
# OR_OUT_DIR = "out_metas/"
# TOT_COMPARISONS = "../Figure1/comparisons_metas.txt"
# MAX_POOL = 20

print(INP_METAD)
print(INP_DATS)
print(OR_OUT_DIR)
print(TOT_COMPARISONS)
print(MAX_POOL)


# stop("Here")

########################
# Need to implement the filters already present in the original code!!!!
########################


tot_dats <- fread(INP_DATS, data.table = FALSE, header = FALSE, sep = ":")

if (! dir.exists(OR_OUT_DIR) ) {
  dir.create(OR_OUT_DIR, showWarnings = FALSE)
}

inp_metadata <- fread(INP_METAD, data.table = FALSE)
dim(inp_metadata)
inp_metadata <- inp_metadata[
  which(
    (inp_metadata$sex != "nd") &
      (!is.na(inp_metadata$age)) &
      (!is.na(inp_metadata$BMI))
  ), ]
dim(inp_metadata)

compute_met <- function(mean_diff, n.e, n.c, sd.e, sd.c) {
  
  npn <- function(x) {
    ##
    ## Check for non-positive values in vector
    ##
    selNA <- is.na(x)
    res <- selNA
    if (sum(!selNA) > 0)
      x[!selNA] <- x[!selNA] <= 0
    ##
    res
  }
  
  npn.n <- npn(n.e) | npn(n.c)
  
  N <- n.e + n.c
  
  S.within <- sqrt(((n.e - 1) * sd.e^2 + (n.c - 1) * sd.c^2) / (N - 2))
  smd <- ifelse(npn.n, NA, (mean_diff) / S.within)
  
  J <- function(x)
    exp(lgamma(x / 2) - log(sqrt(x / 2)) - lgamma((x - 1) / 2))
  ##
  K <- function(x) 1 - (x - 2) / (x * J(x)^2)
  
  TE <- J(N - 2) * smd
  seTE <- ifelse(npn.n, NA, sqrt(1 / n.e + 1 / n.c + TE^2 * K(N - 2)))
  
  return(c(TE, seTE))
  
  
}

estimate_effect_sizes <- function(
    tmp_df, # this must contain the feature to test, the target class, and the covariates to correct for + the study name column
    covariates = c("age", "sex", "BMI"),
    comparison = "",
    feature = ""
    ) {
  
  inp_dat <- data.frame()
  j <- 1
  
  for (study in sort(unique(tmp_df$study_name))) {
    
    # print(study)
    
    ttmp <- tmp_df[which(tmp_df$study_name == study), ]
    
    covariates_to_consider <- c()
    
    for (cov in covariates) {
      
      if (length(unique(ttmp[, cov])) >= 2) {
        covariates_to_consider <- c(covariates_to_consider, cov)
      }

    }
    
    if (length( unique(ttmp$tar_col) ) == 2) {
      
      
      # print(paste0(feature, "~", paste(c("tar_col", covariates_to_consider), collapse = "+")))
      # print(paste0(feature, "~", paste(covariates_to_consider, collapse = "+")))
      
      formula_all_covs <- paste0(feature, "~", paste(c("tar_col", covariates_to_consider), collapse = "+"))
      formula_only_covs <- paste0(feature, "~", paste(covariates_to_consider, collapse = "+"))
      
      fit <- glm(formula(formula_all_covs), data = ttmp, family = "gaussian")
      fit_no_tar_class <- glm(formula(formula_only_covs), data = ttmp, family = "gaussian")
      
      ttmp$residuals_no_tar_class <- residuals(fit_no_tar_class)
      
      inp_dat[j, "study_name"] <- study
      inp_dat[j, "feat_name"] <- feature
      
      inp_dat[j, "formula_all_covs"] <- as.character(formula_all_covs)
      inp_dat[j, "formula_only_covs"] <- formula_only_covs
      
      inp_dat[j, "mean.e"] <- mean(ttmp[which(ttmp$tar_col  == "class2"), feature])
      inp_dat[j, "sd.e"] <- sd(ttmp[which(ttmp$tar_col  == "class2"), feature])
      inp_dat[j, "n.e"] <- length(ttmp[which(ttmp$tar_col  == "class2"), feature])
      
      inp_dat[j, "mean.c"] <- mean(ttmp[which(ttmp$tar_col  == "class1"), feature])
      inp_dat[j, "sd.c"] <- sd(ttmp[which(ttmp$tar_col  == "class1"), feature])
      inp_dat[j, "n.c"] <- length(ttmp[which(ttmp$tar_col  == "class1"), feature])
      
      
      inp_dat[j, "estimate"] <- coef(summary(fit))[2, "Estimate"]
      inp_dat[j, "sd_error"] <- coef(summary(fit))[2, "Std. Error"]
      
      inp_dat[j, "sd.e.glm"] <- sd(ttmp[which(ttmp$tar_col  == "class2"), "residuals_no_tar_class"])
      inp_dat[j, "sd.c.glm"] <- sd(ttmp[which(ttmp$tar_col  == "class1"), "residuals_no_tar_class"])
      
      inp_dat[j, "pval"] <- coef(summary(fit))[2, "Pr(>|t|)"]
      
      inp_dat[j, "n"] <- dim(ttmp)[1]
      
      j <- j + 1
      
    }
    
  }
  
  inp_dat[, c("Hedges_g", "seTE")] <- t(apply(inp_dat[, c("estimate", "n.e", "n.c", "sd.e.glm", "sd.c.glm")], 1, function(x) {
    mean_diff <- x[1]
    n.e <- x[2]
    n.c <- x[3]
    sd.e <- x[4]
    sd.c <- x[5]
    compute_met(mean_diff, n.e, n.c, sd.e, sd.c)
  }))
  
  for (tmp_row in rownames(inp_dat)) {
    
    tmp_smd <- inp_dat[tmp_row, "Hedges_g", drop = TRUE]
    tmp_se <- inp_dat[tmp_row, "seTE", drop = TRUE]
    
    confidence_interval <- ci(tmp_smd, tmp_se, level = 0.95, df = NULL, null.effect = 0)
    
    inp_dat[tmp_row, "lower95ci"] <- confidence_interval$lower
    inp_dat[tmp_row, "upper95ci"] <- confidence_interval$upper
    inp_dat[tmp_row, "p"] <- confidence_interval$p
    
  }
  
  # inp_dat$N <- inp_dat$n.e + inp_dat$n.c
  inp_dat$sd <- inp_dat$seTE * sqrt(inp_dat$n.e + inp_dat$n.c)
  
  # inp_dat[, c("Hedges_g", "seTE")]
  
  # ci_est <- ci(inp_dat_single_study$)
  
  return(inp_dat)
}

perform_usual_metaanalysis <- function(inp_dat) {
  
  # inp_dat <- inp_dat_single_study
  
  m.cont <- metacont(n.e = n.e,
                     mean.e = mean.e,
                     sd.e = sd.e,
                     n.c = n.c,
                     mean.c = mean.c,
                     sd.c = sd.c,
                     studlab = study_name,
                     data = inp_dat,
                     sm = "SMD",
                     method.smd = "Hedges", # Cohen
                     fixed = FALSE,
                     random = TRUE,
                     method.tau = "REML", # DL
                     hakn = TRUE,
                     title = "", control=list(stepadj=0.5, maxiter=10000))
  
  object <-summary(m.cont)
  
  ci.r <- data.frame(summary_measure = object$TE.random,
                     lower95ci = object$lower.random,
                     upper95ci = object$upper.random,
                     p = object$pval.random,
                     heterogeneity = object$I2,
                     pval_heterogeneity = object$pval.Q,
                     w_per_random = NA,
                     study = "model"
  )
  
  ci.study <- data.frame(summary_measure = object$TE,
                         lower95ci = object$lower,
                         upper95ci = object$upper,
                         p = object$pval,
                         heterogeneity = NA,
                         pval_heterogeneity = NA,
                         w_per_random = 100 * object$w.random / sum(object$w.random, na.rm = TRUE),
                         study = object$studlab
  )
  
  out_meta <- rbind(ci.study, ci.r)
  out_meta$sm = object$sm
  
  return(list(object, out_meta))
  
}

for (row in rownames(tot_dats)) {
  print(row)
  
  print(tot_dats[row, ])  
  
  CURR_NAME = tot_dats[row, "V1"]
  INP_TABLE = tot_dats[row, "V2"]
  
  OUT_DIR = file.path(OR_OUT_DIR, CURR_NAME)
  
  inp_profiles = fread(INP_TABLE, data.table = FALSE); inp_profiles[1:3,1:3]
  
  rownames(inp_profiles) <- gsub("\\|", "___", inp_profiles$feat_name); inp_profiles$feat_name <- NULL
  inp_profiles <- as.data.frame(t(inp_profiles))
  inp_profiles$sample_id <- rownames(inp_profiles)
  
  print(OUT_DIR)
  
  if (! dir.exists(OUT_DIR) ) {
    dir.create(OUT_DIR, showWarnings = FALSE)
  }
  
  comparisons <- fread(TOT_COMPARISONS, data.table = FALSE, header = FALSE, sep = ":")
  
  
  
  for (comp_row in rownames(comparisons)) {
    
    study_col <- comparisons[comp_row, "V1"]
    tar_col <- comparisons[comp_row, "V2"]
    comb <- comparisons[comp_row, "V3"]
    MIN_N_STUDES <- comparisons[comp_row, "V3"]
      
      print(comb)
    
    class1 <- unlist(strsplit(unlist(strsplit(comb, split = "_vs_"))[1], split = "__"))
    class2 <- unlist(strsplit(unlist(strsplit(comb, split = "_vs_"))[2], split = "__"))
    
    tmp_metadata <- inp_metadata[
      which(inp_metadata[, tar_col] %in% c(class1, class2)), 
      ]
    
    colnames(tmp_metadata)[which(colnames(tmp_metadata) == tar_col)] <- "tar_col"
    
    tmp_metadata[
      which(tmp_metadata$tar_col %in% c(class1)), "tar_col"
    ] <- "class1"
    
    tmp_metadata[
      which(tmp_metadata$tar_col %in% c(class2)), "tar_col"
    ] <- "class2"
    
    dim(tmp_metadata)
    dim(inp_profiles)
    
    tmp_metadata <- inner_join(tmp_metadata, inp_profiles, by = "sample_id")
    
    features_to_test <- setdiff(colnames(inp_profiles), "sample_id")
    
    tot_out_mets <- list()
    tot_out_mets_w_corrections <- list()
    
    j <- 1
    
    for (feat in features_to_test) {

      # feat <- "oral2gut_score"
      
      # tmp_df <- tmp_metadata
      # covariates <- c("age", "sex", "BMI")
      # comparison <- comb
      # feature <- feat
        
        valid_studies <- c()
        
        for (study in unique(tmp_metadata$study)) {
            
            check1 <- tmp_metadata[
                which(
                    (tmp_metadata$study == study) & 
                    (tmp_metadata$tar_col == "class1")
                ), ]
            pos_hits1 <- length(which(check1[, feat] > 0))
            prev_feat1 <- length(which(check1[, feat] > 0)) / dim(check1)[1]
            
            check2 <- tmp_metadata[
                which(
                    (tmp_metadata$study == study) & 
                    (tmp_metadata$tar_col == "class2")
                ), ]
            pos_hits2 <- length(which(check2[, feat] > 0))
            prev_feat2 <- length(which(check2[, feat] > 0)) / dim(check2)[1]
            
            
            cond1 <- (dim(check1)[1] >= 10) & (dim(check2)[1] >= 10)
            cond2 <- (prev_feat1 >= 0.1) | (prev_feat2 >= 0.1)
            cond3 <- (pos_hits1 >= 5) | (pos_hits2 >= 5)
            
            if (cond1 & cond2 & cond3) {
                valid_studies <- c(valid_studies, study)
            }
            
            
        }
        
    if (length(valid_studies) >= 3) {
        
        tmp_metadata <- tmp_metadata[
            which(tmp_metadata$study %in% valid_studies), 
        ]
        
        inp_dat_single_study <- estimate_effect_sizes(
            tmp_metadata, ### 
            covariates = c("age", "sex", "BMI"),
            comparison = comb,
            feature = feat
            )

          ### Usual meta-analysis

          out_meta <- perform_usual_metaanalysis(inp_dat_single_study)
          out_meta[[2]]$feat_name <- feat
          out_meta[[2]]$comp <- comb

          ###

          out <- metamean(
            n,
            Hedges_g,
            sd, study_name, data = inp_dat_single_study,
            fixed = FALSE,
            random = TRUE,
            method.tau = "REML",
            hakn = TRUE, null.effect = 0,
            sm = "MRAW",
            control=list(stepadj=0.5, maxiter=10000)
            )
          objs <- summary(out)

          to_join <- as.data.frame(objs)[, c("studlab", "w.random")]
          colnames(to_join)[which(colnames(to_join) == "studlab")] <- "study_name"

          out_metas_w_correction <- as.data.frame(inp_dat_single_study) %>% inner_join(to_join, by = "study_name")

          colnames(out_metas_w_correction)[which(colnames(out_metas_w_correction) == "estimate")] <- "estimate.glm"
          colnames(out_metas_w_correction)[which(colnames(out_metas_w_correction) == "sd_error")] <- "sd_error.glm"
          colnames(out_metas_w_correction)[which(colnames(out_metas_w_correction) == "pval")] <- "pval.glm"

          colnames(out_metas_w_correction)[which(colnames(out_metas_w_correction) == "Hedges_g")] <- "summary_measure"
          colnames(out_metas_w_correction)[which(colnames(out_metas_w_correction) == "seTE")] <- "se"


          metas_corrected <- data.frame(
            "study_name" = "model",
            "feat_name" = feat,
            "formula_all_covs" = NA,
            "formula_only_covs" = NA,
            "n.e" = NA,
            "n.c" = NA,
            "n" = NA,
            "summary_measure" = objs$TE.random,
            "seTE" = NA,
            "lower95ci" = objs$lower.random,
            "upper95ci" = objs$upper.random,
            "p" = objs$pval.random,

            "heterogeneity" = objs$I2,
            "pval_heterogeneity" = objs$pval.Q
          )

          out_metas_w_correction <- dplyr::bind_rows(out_metas_w_correction, metas_corrected)

          out_metas_w_correction$sm <- "Hedge's g"
          out_metas_w_correction$comp <- comb

          # Output to save

          tot_out_mets[[j]] <- out_meta[[2]]
          tot_out_mets_w_corrections[[j]] <- out_metas_w_correction

          j <- j + 1
        
        
        
    }
      
      
      
      
      # break
      
    }
    
    tot_out_mets <- dplyr::bind_rows(tot_out_mets)
    tot_out_mets_w_corrections <- dplyr::bind_rows(tot_out_mets_w_corrections)
    
    tot_out_mets %>% write_tsv(file.path(OUT_DIR, paste0("out_metas_nocovariates__", comb, ".tsv")))
    tot_out_mets_w_corrections %>% write_tsv(file.path(OUT_DIR, paste0("out_metas_withcovariates__", comb, ".tsv")))
    
    dim(tmp_metadata)
    
    
    # break
    
  }
  
  
}

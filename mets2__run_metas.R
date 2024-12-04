suppressMessages(library("tidyverse"))
suppressMessages(library("meta"))

suppressMessages(library("data.table"))

library("parallel")

args = commandArgs(trailingOnly=TRUE)

# Rscript run_metas2.R out_metas/dat_met.tsv out_metas/
# Rscript run_metas2.R out_metas_divs/dat_met.tsv out_metas_divs/


inp_dat <- args[1] # "out_metas/dat_met.tsv"
out_dir <- args[2] # "out_metas/"
TOT_COMPARISONS <- args[3] # 3 # the number of minimum number of studies

print(inp_dat)
print(out_dir)

# stop("Here")

# SIGNIF_THR = 0.01
# FDR_THR = 0.2
# TOP_LOW_RANKED = 15

PREV_THR = 0.1
# MIN_N_STUDES=5
MIN_N_CLASS = 10
MIN_N_POSITIVE_HITS = 5

N_CORES <- 10

n_dat_met <- fread(inp_dat, data.table=FALSE)     

print("Here")

TOT_COMPARISONS_df <- fread(TOT_COMPARISONS, data.table = FALSE, sep = ":")
            
print(TOT_COMPARISONS_df)
# stop("Here")
                      
#########################################
# BLOCK 2
#########################################

# Filter based on the sample size
n_dat_met <- n_dat_met[which((n_dat_met$n.e >= MIN_N_CLASS) & (n_dat_met$n.c >= MIN_N_CLASS)),]
# Filter based on the prevalence
n_dat_met <- n_dat_met[which((n_dat_met$prev.c >= PREV_THR) | (n_dat_met$prev.e >= PREV_THR)),]
# Filter based on the minimum number of positive hits in at least one sample
n_dat_met <- n_dat_met[which((n_dat_met$hits.c >= MIN_N_POSITIVE_HITS) | (n_dat_met$hits.e >= MIN_N_POSITIVE_HITS)),]

# Filter based on the prevalence
# n_dat_met <- n_dat_met[which((n_dat_met$sd.e > 0) & (n_dat_met$sd.c > 0)),]

out_metas_metacont <- data.frame()
out_metas_metabin <- data.frame()

warnings_errors_metacont <- data.frame()
warnings_errors_metabin <- data.frame()

return_error_warning_out_metas <- function(
    comp,
    index,
    n_curr_tmp,
    err,
    message
) {
    
    out <- data.frame(
        "summary_measure" = NA,
        "lower95ci" = NA,
        "upper95ci" = NA,
        "p" = NA,
        "heterogeneity" = NA,
        "pval_heterogeneity" = NA,
        "w_per_random" = NA,
        "study" = "model",
        "comp" = comp,
        "feat" = index,
        "message" = message,
        "or_dim" = dim(n_curr_tmp)[1],
        "extended_message" = as.character(err),
        "sm" = NA
    )
    
    return(out)
    
}

run_metacont <- function(list_input) {
    
    # list("index" = "", "curr_tmp" = "", "MIN_N_STUDES" = 3)
    
    index <- list_input[["index"]]
    curr_tmp <- list_input[["curr_tmp"]]
    MIN_N_STUDES <- list_input[["MIN_N_STUDES"]]
    
    n_curr_tmp <- curr_tmp[which(curr_tmp$feat_name == index), ]
    n_curr_tmp_metacont <- n_curr_tmp[which((n_curr_tmp$sd.e > 0) & (n_curr_tmp$sd.c > 0)),]
        
    if (dim(n_curr_tmp_metacont)[1] >= MIN_N_STUDES) {

        # Perform metacont

        tryCatch(
                expr = {
                    m.cont <- metacont(n.e = n.e,
                               mean.e = mean.e,
                               sd.e = sd.e,
                               n.c = n.c,
                               mean.c = mean.c,
                               sd.c = sd.c,
                               studlab = study_name,
                               data = n_curr_tmp_metacont,
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
                    out_meta$comp <- comp
                    out_meta$feat <- index
                    out_meta$sm = object$sm
                    
                    out_meta$message = "correct"
                    out_meta$or_dim = dim(n_curr_tmp_metacont)[1]
                    out_meta$extended_message = as.character("correct")

                    # out_metas_metacont <- rbind(out_metas_metacont, out_meta)
                    return(out_meta)

                },
                error = function(e){ 
                    # cat("\nError [metacont]: ", index)
                    # warnings_errors_metacont[index, "type"] <- "error"
                    # warnings_errors_metacont[index, "message"] <- e
                    
                    out_meta <- return_error_warning_out_metas(
                        comp,
                        index,
                        n_curr_tmp,
                        e,
                        "error"
                    )

                    return(out_meta)
                    
                },
                warning = function(w){
                    # cat("\nWarning [metacont]: ", index)
                    # warnings_errors_metacont[index, "type"] <- "warning"
                    # warnings_errors_metacont[index, "message"] <- w
                    
                    out_meta <- return_error_warning_out_metas(
                        comp,
                        index,
                        n_curr_tmp,
                        w,
                        "warning"
                    )

                    return(out_meta)
                }

            )
        }
    
     else {
        
        out_meta <- return_error_warning_out_metas(
                comp,
                index,
                n_curr_tmp,
                "N. dataset < MIN N dataset",
                "error"
            )
        
        return(out_meta)
        
    }
    
    
}

run_metabin <- function(list_input) {
    
    # list_input <- list("index" = "", "curr_tmp" = "", "MIN_N_STUDES" = 3)
    
    index <- list_input[["index"]]
    curr_tmp <- list_input[["curr_tmp"]]
    MIN_N_STUDES <- list_input[["MIN_N_STUDES"]]
    
    n_curr_tmp <- curr_tmp[which(curr_tmp$feat_name == index), ]
    n_curr_tmp_metacont <- n_curr_tmp[which((n_curr_tmp$sd.e > 0) & (n_curr_tmp$sd.c > 0)),]
    
    if (dim(n_curr_tmp)[1] >= MIN_N_STUDES) {

        # Perform metabin

        tryCatch(
                expr = {
                    m.bin <- metabin(event.e = hits.e,
                                      n.e = n.e,
                                      event.c = hits.c,
                                      n.c = n.c,
                                      studlab = study_name,
                                      data = n_curr_tmp,
                                      sm = "RR",
                                      fixed = FALSE,
                                      random = TRUE,
                                      title = "", control=list(stepadj=0.5, maxiter=10000))
                object <-summary(m.bin)


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
                                       study = object$studlab)

                out_meta <- rbind(ci.study, ci.r)
                out_meta$comp <- comp
                out_meta$feat <- index
                out_meta$sm = object$sm
                
                out_meta$message = "correct"
                out_meta$or_dim = dim(n_curr_tmp)[1]
                out_meta$extended_message = as.character("correct")
                    
                return(out_meta)

                # out_metas_metabin <- rbind(out_metas_metabin, out_meta)

                },
                error = function(e){ 
                    # cat("\nError [metabin]: ", index)
                    # warnings_errors_metacont[index, "type"] <- "error"
                    # warnings_errors_metacont[index, "message"] <- e
                    
                    out_meta <- return_error_warning_out_metas(
                        comp,
                        index,
                        n_curr_tmp,
                        e,
                        "error"
                    )

                    return(out_meta)
                    
                },
                warning = function(w){
                    # cat("\nWarning [metabin]: ", index)
                    # warnings_errors_metacont[index, "type"] <- "warning"
                    # warnings_errors_metacont[index, "message"] <- w
                    
                    out_meta <- return_error_warning_out_metas(
                        comp,
                        index,
                        n_curr_tmp,
                        w,
                        "warning"
                    )
                    
                }

            )
        } else {
        
        out_meta <- return_error_warning_out_metas(
                comp,
                index,
                n_curr_tmp,
                "N. dataset < MIN N dataset",
                "error"
            )
        
        return(out_meta)
        
    }
    
    
}


for (comp in unique(n_dat_met$comparison)) {
    
    cat("\n --- --- --- --- --- \n", comp, "\n --- --- --- --- --- \n")
    
    # This line is important to detect the number of minimum studies to allow in the meta-analysis
    MIN_N_STUDES <- TOT_COMPARISONS_df[which(TOT_COMPARISONS_df$V3 == comp), "V4", drop = TRUE]
    
    cat("\n --> MIN_N_STUDES: ", MIN_N_STUDES, "\n")
    
    curr_tmp <- n_dat_met[which(n_dat_met$comparison == comp), ]
    
    out_metas_metacont <- mclapply(
        unique(curr_tmp$feat_name), function(x) { 
            list("feat_name" = x, "output" = run_metacont(list("index" = x, "curr_tmp" = curr_tmp, "MIN_N_STUDES" = MIN_N_STUDES)))
        }, mc.cores = N_CORES)
    # out_metas_metabin <- mclapply(
    #     unique(curr_tmp$feat_name), function(x) { 
    #         list("feat_name" = x, "output" = run_metabin(list("index" = x, "curr_tmp" = curr_tmp, "MIN_N_STUDES" = MIN_N_STUDES)))
    #     }, mc.cores = N_CORES)
    
    out_metas_metacont <- lapply(out_metas_metacont, function(x) {return(x[["output"]])})
    # out_metas_metabin <- lapply(out_metas_metabin, function(x) {return(x[["output"]])})    
    
    out_metas_metacont <- do.call(rbind, out_metas_metacont)
    # out_metas_metabin <- do.call(rbind, out_metas_metabin)
    
    out_metas_metacont_warnings_errors <- out_metas_metacont[which(out_metas_metacont$message != "correct"), ]
    # out_metas_metabin_warnings_errors <- out_metas_metabin[which(out_metas_metabin$message != "correct"), ]
    
    out_metas_metacont <- out_metas_metacont[which(out_metas_metacont$message == "correct"), ]
    # out_metas_metabin <- out_metas_metabin[which(out_metas_metabin$message == "correct"), ]
    
    # Perform FDR correction for the metacont meta-analysis

    for (study in unique(out_metas_metacont$study)) {

        out_metas_metacont[
            which(
                (out_metas_metacont$study == study)
            ), 
            "FDR"] <- p.adjust(
            out_metas_metacont[
                which(
                (out_metas_metacont$study == study)
                ), "p"], method="BH")

    }

    write_tsv(out_metas_metacont, paste0(out_dir, "out_metas_metacont___", comp, ".tsv"))

    # Perform FDR correction for the metabin meta-analysis

#     for (study in unique(out_metas_metabin$study)) {

#         out_metas_metabin[
#             which(
#                 (out_metas_metabin$study == study)
#             ), 
#             "FDR"] <- p.adjust(
#             out_metas_metabin[
#                 which(
#                 (out_metas_metabin$study == study)
#                 ), "p"], method="BH")
#     }
    
#     write_tsv(out_metas_metabin, paste0(out_dir, "out_metas_metabin___", comp, ".tsv"))

    # Write the tables with warnings and errors

    write_tsv(warnings_errors_metacont, paste0(out_dir, "warnings_errors_metacont___", comp, ".tsv"))
    # write_tsv(warnings_errors_metabin, paste0(out_dir, "warnings_errors_metabin___", comp, ".tsv"))
    

}


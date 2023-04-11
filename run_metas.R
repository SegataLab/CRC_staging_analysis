suppressMessages(library("tidyverse"))
suppressMessages(library("meta"))

library("parallel")

library("optparse")

parser <- OptionParser()
parser <- add_option(parser, c("-i", "--input"), type="character", 
                     default="../data_March23/metaph4_formatted_asinsqrt.tsv")
parser <- add_option(parser, c("-n", "--name"), type="character", 
                     default="mpa4")
parser <- add_option(parser, c("-m", "--metadata"), type="character", 
                     default="../data_March23/metadata_formatted_Mar23.tsv")
# parser <- add_option(parser, c("-c", "--comparisons"), type="character", 
#                      default="")
parser <- add_option(parser, c("-o", "--output_folder"), type="character", 
                     default="out_metas/")

arguments <- parse_args(parser)

# "mpa4__": "../data_Jan23/metaph4_formatted_sgbs.tsv",
# "mpa4_oralsgbs__": "../data_Jan23/metaph4_formatted_sgbs_filt_oral.tsv",
# "hnn36_ec_table__": "../data_Jan23/hnn36_ec_table_formatted_renamed.tsv",
# "hnn36_pathways__": "../data_Jan23/hnn36_pathway_formatted.tsv"

# Metaphlan-4
# out_name <- "mpa4"
# inp_path <- "../data_March23/metaph4_formatted_asinsqrt.tsv"
# inp_met_path <- "../data_March23/metadata_formatted_Mar23.tsv"

# Metaphlan-4 with other levels
# out_name <- "mpa4_wspecies_level"
# inp_path <- "../data_March23/metaph4_formatted_asinsqrt_wspecies_level.tsv"
# inp_met_path <- "../data_March23/metadata_formatted_Mar23.tsv"

# Gene families - EC [UR90]
# out_name <- "hnn36_ec_table"
# inp_path <- "../data_March23/hnn36_ecs_formatted_asinsqrt.tsv"
# inp_met_path <- "../data_March23/metadata_formatted_Mar23_4functional_profiles.tsv"

# Pathway abundances [UR90]
# out_name <- "hnn36_pathways"
# inp_path <- "../data_March23/hnn36_paths_formatted_asinsqrt.tsv"
# inp_met_path <- "../data_March23/metadata_formatted_Mar23_4functional_profiles.tsv"



# out_dir <- "out_metas/"

SIGNIF_THR = 0.01
FDR_THR = 0.2
TOP_LOW_RANKED = 15
PREV_THR = 0.1

MIN_N_STUDES=3

MIN_N_CLASS = 10

inp_path <- arguments$input
out_name <- arguments$name
inp_met_path <- arguments$metadata
# comparisons <- unlist(strsplit(arguments$comparisons, split=";"))
out_dir <- arguments$output_folder

if (!dir.exists(out_dir)) {
    
    dir.create(out_dir)
    
}

inp_profiles <- as.data.frame(read_tsv(inp_path, show_col_types=TRUE))
metadata_table <- as.data.frame(read_tsv(inp_met_path, show_col_types=FALSE))

comparisons <- c(
    "study_name:stage_class:control_vs_0__I__II__III__IV__crc_missing_stage",
    "study_name:stage_class:0__I__II_vs_III:earlyCRC_vs_stageIII",
    "study_name:stage_class:0__I__II_vs_IV:earlyCRC_vs_stageIV",
    "study_name:stage_class:0__I__II_vs_III__IV:earlyCRC_vs_lateCRC",
    "study_name:stage_class:0__I__II__III_vs_IV:nonmetCRC_vs_metCRC",
    "study_name:tumor_location_class:right_sided_vs_left_sided:right_vs_left",
    "study_name:tumor_location_class:right_sided__transverse_vs_left_sided:rightANDtransverse_vs_left",
    "study_name:tumor_location_class:right_sided_vs_transverse__left_sided:right_vs_transverseANDleft"
    )

print("---- Here1")
#########################################
# BLOCK 1
#########################################

dat_met <- data.frame()

for (comparison in comparisons) {
    
    print(comparison)
    
    study_col <- unlist(strsplit(comparison, split=":"))[1]
    tar_col <- unlist(strsplit(comparison, split=":"))[2]
    comb <- unlist(strsplit(comparison, split=":"))[3]
    
    comb1 <- unlist(
        strsplit(
            unlist(
                strsplit(
                    comb, split="_vs_"
                )
            )[1], split="__")
        )
    comb2 <- unlist(
        strsplit(
            unlist(
                strsplit(
                    comb, split="_vs_"
                )
            )[2], split="__")
        )

    tmp_pmeta <- metadata_table[which(metadata_table[, tar_col] %in% c(comb1, comb2) ), ]
    
    if (tar_col == "tumor_location_class") {
        
        tmp_pmeta <- tmp_pmeta[which(tmp_pmeta$stage_class %in% c(
            "0", "I", "II", "III", "IV", "crc_missing_stage"
        )), ]
    }

    for (dat in sort(unique(tmp_pmeta[, study_col, drop=TRUE]))) {
        samples1 <- tmp_pmeta[which(
            (tmp_pmeta[, study_col] == dat) &
            (tmp_pmeta[, tar_col] %in% comb1)), "sample_id"]
        samples2 <- tmp_pmeta[which(
            (tmp_pmeta[, study_col] == dat) &
            (tmp_pmeta[, tar_col] %in% comb2)), "sample_id"]
        
        samples1 <- intersect(samples1, colnames(inp_profiles))
        samples2 <- intersect(samples2, colnames(inp_profiles))
        
        
        if ( (length(samples1) > MIN_N_CLASS) & (length(samples2) > MIN_N_CLASS)) {
            

#             cat("\n\n", dat, length(samples1), length(samples2))

            o_features_to_consider <- inp_profiles[, c(samples1, samples2)]
            
            o_features_to_consider1 <- inp_profiles[, c("feat_name", samples1)]
            features_to_consider1 <- o_features_to_consider1[
                apply(o_features_to_consider1, 1, function(x) sum(x > 0)/length(x) > PREV_THR), "feat_name", drop=TRUE]
            
            o_features_to_consider2 <- inp_profiles[, c("feat_name", samples2)]
            features_to_consider2 <- o_features_to_consider2[
                apply(o_features_to_consider2, 1, function(x) sum(x > 0)/length(x) > PREV_THR), "feat_name", drop=TRUE]
            
            features_to_consider <- unique(c(features_to_consider1, features_to_consider2))
            
            cat("\n\n--------------\n",
                "\nDataset:", dat,
                "\nN. samples1:", length(samples1),
                "\nN. samples2:", length(samples2),
                "\nN. feats to consider 1:", length(features_to_consider1)[1],
                "\nN. feats to consider 2:", length(features_to_consider2)[1],
                "\nN. tot feats to consider:", length(features_to_consider))
                      
                      
            tmp_out_df <- mclapply(features_to_consider, function(x) {
                
                tryCatch(
                    {
                        n.e <- length(samples2)
                        mean.e <- mean( as.matrix(inp_profiles[which(inp_profiles$feat_name == x), samples2]) )
                        sd.e <- sd( as.matrix(inp_profiles[which(inp_profiles$feat_name == x), samples2]) )

                        n.c <- length(samples1)
                        mean.c <- mean( as.matrix(inp_profiles[which(inp_profiles$feat_name == x), samples1]) )
                        sd.c <- sd( as.matrix(inp_profiles[which(inp_profiles$feat_name == x), samples1]) )

                        dat_met <- data.frame(
                            "feat_name"=x,
                            "n.e"=n.e,
                            "mean.e"=mean.e,
                            "sd.e"=sd.e,
                            "n.c"=n.c,
                            "mean.c"=mean.c,
                            "sd.c"=sd.c,
                            "study_name"=dat,
                            "comparison"=comb)

                        return(dat_met)
                    
                    },
                    warning = function(w) {
                        dat_met <- data.frame(
                            "feat_name"=x,
                            "n.e"="warning",
                            "mean.e"="warning",
                            "sd.e"="warning",
                            "n.c"="warning",
                            "mean.c"="warning",
                            "sd.c"="warning",
                            "study_name"=dat,
                            "comparison"=comb)
                        return(dat_met)
                    },
                    error = function(e) {
                        dat_met <- data.frame(
                            "feat_name"=x,
                            "n.e"="error",
                            "mean.e"="error",
                            "sd.e"="error",
                            "n.c"="error",
                            "mean.c"="error",
                            "sd.c"="error",
                            "study_name"=dat,
                            "comparison"=comb)
                        return(dat_met)
                    }
                )

                

            }, mc.cores = 20)
            dat_met <- rbind(dat_met, do.call("rbind", tmp_out_df))
            # break
        }
    }
}
                      
print("---- Here2")
                      
                      
#########################################
# BLOCK 2
#########################################
                      
n_dat_met <- na.omit(dat_met)
n_dat_met <- n_dat_met[which((n_dat_met$sd.e > 0) & (n_dat_met$sd.c > 0)),]

out_metas <- data.frame()

for (comp in unique(n_dat_met$comparison)) {
    
    print(comp)
    
    curr_tmp <- n_dat_met[which(n_dat_met$comparison == comp), ]
    
    # df <- Reduce(rbind, listOfDataFrames)

    for (index in unique(curr_tmp$feat_name)) {

        # print(index)
        
        n_curr_tmp <- curr_tmp[which(curr_tmp$feat_name == index), ]
        
        if (dim(n_curr_tmp)[1] >= MIN_N_STUDES) {
            tryCatch(
                    expr = {
                        m.cont <- metacont(n.e = n.e,
                                   mean.e = mean.e,
                                   sd.e = sd.e,
                                   n.c = n.c,
                                   mean.c = mean.c,
                                   sd.c = sd.c,
                                   studlab = study_name,
                                   data = n_curr_tmp,
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
                                           study = "model")

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

                        out_metas <- rbind(out_metas, out_meta)

                    },
                    error = function(e){ 
                        cat("\nError: ", index)
                    },
                    warning = function(w){
                        cat("\nWarning: ", index)
                    }
                    # ,
                    # finally = {
                    #     cat("\nLast iteration: ", index)
                    # }
                )
            }
        }
}
                      
# out_metas$p[is.na(out_metas$p)]
out_metas[which(out_metas$study == "model"), "FDR"] <- p.adjust(out_metas[which(out_metas$study == "model"), "p"], method="BH")
out_metas[which(out_metas$study != "model"), "FDR"] <- p.adjust(out_metas[which(out_metas$study != "model"), "p"], method="BH")

print("Here")

out_metas$signif <- unlist(lapply(out_metas$FDR, function(x) {
    if (x <= 0.2) {
        return("yes")
    } else {
        return("no")
    }
}))

write_tsv(dat_met, paste0(out_dir, "dat_met_", out_name, ".tsv"))
write_tsv(out_metas, paste0(out_dir, "out_metas_", out_name, ".tsv"))
library("optparse")

library(tidyverse)
library(data.table)
library(parallel)

library(vegan)

# Example usage:
# conda activate gpiccinno_oncobiome
# Rscript ../master_scripts/PERMANOVA_from_distance.R -i bray_curtis_metaph4_only_sgbs__0_1.rds -m ../metadata_20241017/metadata_CRC_staging_20241017.tsv -c config_permanova.tsv -n 20 -o out_PERMANOVA__SGBs

read_params <- function () {

    option_list = list(
        make_option(c("-i", "--distance_object"), type="character", default=NULL, 
                    help="", metavar="character"),
        make_option(c("-m", "--metadata_table"), type="character", default=NULL, 
                    help="", metavar="character"),
        # make_option(c("-s", "--level_prefix"), type="character", default="t__", 
        #             help="", metavar="character"),
        # make_option(c("-d", "--distance"), type="character", default="bray", 
        #             help="", metavar="character"),
        make_option(c("-c", "--config_file"), type="character", default=NULL, 
                    help="", metavar="character"),
        make_option(c("-n", "--ncores"), type="numeric", default=1, 
                    help=""),
        make_option(c("-o", "--output_folder"), type="character", default="curr_tmp_output_folder", 
                    help="", metavar="character")
    ); 

    opt_parser = OptionParser(option_list=option_list);
    opt = parse_args(opt_parser);

    return(opt)
    
}

main <- function() {
    
    params <- read_params()
    
    print(params)
    
    curr_dist_profiles <- as.matrix(readRDS(params$distance_object))
    metadata_table <- fread(params$metadata_table, data.table = FALSE)
    
    ### This is just for testing
    metadata_table <- metadata_table[which(grepl("Wirbel|Zeller", metadata_table$study_name)), ]
    ### To remove after testing

    classes_to_test <- fread(params$config_file, data.table = FALSE)$target_comps
    
    print(classes_to_test)
    
    ## Check that the distance is among our options
    
    # if (! distance %in% c("bray", "jaccard")) {
    #   stop("Provided distance not in allowed options!")
    # }
    
    ####
    
    set.seed(42)    

    tot_permanovas <- list()
    
    cat("\n\n ---- START PERMANOVAs")

    for (class_to_compare in classes_to_test) {


        print(class_to_compare)
        target_col <- gsub(":.*", "", class_to_compare)
        
        print(target_col) #; stop("Here")

        class1 <- unlist(
            strsplit(
                    unlist(
                        strsplit(
                            gsub(".*:", "", class_to_compare), split="_vs_"
                        )
                    )[1], split="__")
                )
        class2 <- unlist(
            strsplit(
                unlist(
                    strsplit(
                        gsub(".*:", "", class_to_compare), split="_vs_"
                    )
                )[2], split="__")
        )

        cat("\n ---------- \n", class_to_compare, "-", target_col, "-", class1, "-", class2)

        ### PERFORM the PERMANOVA without covariates
        
        met_to_test <- as.data.frame(metadata_table)
        rownames(met_to_test) <- met_to_test$sample_id

        met_to_test[which(met_to_test[, target_col] %in% class1), target_col] <- "class1"
        met_to_test[which(met_to_test[, target_col] %in% class2), target_col] <- "class2"

        met_to_test <- met_to_test[which(met_to_test[, target_col] %in% c("class1", "class2")), ]
        print(table(met_to_test[, target_col]))
        
        n_curr_dist_profiles <- as.dist(curr_dist_profiles[rownames(met_to_test), rownames(met_to_test)])

        print(table(met_to_test[, target_col]))

        print(dim(met_to_test))
        # print(as.matrix(n_curr_dist_profiles))
        
        colnames(met_to_test)[which(colnames(met_to_test) == target_col)] <- "target_col"

        perm <- how(nperm = 999)
        setBlocks(perm) <- with(met_to_test, study_name)


        cat("\n\n --- --- --- ", class_to_compare, " --- --- --- \n")
        set.seed(42)
        permanova <- adonis2(n_curr_dist_profiles ~ target_col,
            data = met_to_test,
            permutations=perm,
            method = distance, parallel = params$ncores)
        
        permanova <- as.data.frame(permanova)[1, ]

        permanova$betadisper_P <- as.data.frame(anova(betadisper(n_curr_dist_profiles, met_to_test$target_col)))[1, "Pr(>F)"]
        
        permanova$w_covariates <- "no"
        
        print(permanova)

        
        
        ### PERFORM the PERMANOVA with covariates
        
        met_to_test <- as.data.frame(metadata_table[
          which(
            (!(is.na(metadata_table$age))) &
              (metadata_table$sex != "nd") &
              (!(is.na(metadata_table$BMI)))
            ), 
        ])
        rownames(met_to_test) <- met_to_test$sample_id
        
        met_to_test[which(met_to_test[, target_col] %in% class1), target_col] <- "class1"
        met_to_test[which(met_to_test[, target_col] %in% class2), target_col] <- "class2"
        
        met_to_test <- met_to_test[which(met_to_test[, target_col] %in% c("class1", "class2")), ]
        print(table(met_to_test[, target_col]))
        
        n_curr_dist_profiles <- as.dist(curr_dist_profiles[rownames(met_to_test), rownames(met_to_test)])
        
        print(table(met_to_test[, target_col]))
        
        print(dim(met_to_test))
        # print(as.matrix(n_curr_dist_profiles))
        
        colnames(met_to_test)[which(colnames(met_to_test) == target_col)] <- "target_col"
        
        perm <- how(nperm = 999)
        setBlocks(perm) <- with(met_to_test, study_name)
        
        
        cat("\n\n --- --- --- ", class_to_compare, " --- --- --- \n")
        set.seed(42)
        permanova2 <- adonis2(n_curr_dist_profiles ~ age + sex + BMI + target_col,
                             data = met_to_test,
                             permutations=perm,
                             method = distance, parallel = params$ncores)
        
        permanova2 <- as.data.frame(permanova2)[1:4, ]
        
        permanova2$betadisper_P <- NA
        # permanova2$betadisper_P <- as.data.frame(anova(betadisper(n_curr_dist_profiles, met_to_test$target_col)))["target_col", "Pr(>F)"]
        
        permanova2$w_covariates <- "yes"
        
        print(permanova2)
        
        
        # save the output
        
        tmp_to_save <- rbind(permanova, permanova2)
        tmp_to_save$comp <- class_to_compare
        
        tmp_to_save$target_col <- rownames(tmp_to_save)
        
        
        tot_permanovas[[class_to_compare]] <- tmp_to_save
        
        # break

    }
    
    tot_permanovas <- dplyr::bind_rows(tot_permanovas)
    
    # out_r2_pvals_permanova <- data.frame()
    # 
    # j <- 1
    # for (name in names(tot_permanovas)) {
    # 
    #     out_r2_pvals_permanova[j, "comp"] <- gsub(".*:", "", name)
    #     out_r2_pvals_permanova[j, "r2"] <- tot_permanovas[[name]][1, "R2"]
    #     out_r2_pvals_permanova[j, "F"] <- tot_permanovas[[name]][1, "F"]
    #     out_r2_pvals_permanova[j, "pval"] <- tot_permanovas[[name]][1, "Pr(>F)"]
    # 
    #     j <- j + 1
    # 
    # }

    tot_permanovas %>% write_tsv(paste0(params$output_folder, "out_permanova_analysis.tsv"))
    
}



main()
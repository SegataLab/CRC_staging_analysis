library("optparse")

library(tidyverse)
library(data.table)
library(parallel)

library(vegan)

# Example usage:
# conda activate gpiccinno_oncobiome
# Rscript PERMANOVA.R -i ../formatted_data_20241028/metaph4_only_sgbs__0_1.tsv -m ../metadata_20241017/metadata_CRC_staging_20241017.tsv -c config_permanova.tsv -n 20 -o out_PERMANOVA__SGBs

read_params <- function () {
    
    
    option_list = list(
        make_option(c("-i", "--merged_table_profiles_and_strains"), type="character", default=NULL, 
                    help="", metavar="character"),
        make_option(c("-m", "--metadata_table"), type="character", default=NULL, 
                    help="", metavar="character"),
        make_option(c("-s", "--level_prefix"), type="character", default="t__", 
                    help="", metavar="character"),
        make_option(c("-d", "--distance"), type="character", default="bray", 
                    help="", metavar="character"),
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
    
    metaphlan_table_w_strains <- fread(params$merged_table_profiles_and_strains, data.table = FALSE)
    metadata_table <- fread(params$metadata_table, data.table = FALSE)
    
    ### This is just for testing
    metadata_table <- metadata_table[which(grepl("Wirbel|Zeller", metadata_table$study_name)), ]
    ### To remove after testing
    
    print(dim(metaphlan_table_w_strains))
    
    only_strains <- metaphlan_table_w_strains[which(grepl(params$level_prefix, metaphlan_table_w_strains$feat_name)), ]
    rownames(only_strains) <- only_strains$feat_name; only_strains$feat_name <- NULL
    
    profiles_to_test <- as.data.frame(t(only_strains))
    print(dim(profiles_to_test))
    profiles_to_test <- profiles_to_test[which(rowSums(profiles_to_test) > 0), which(colSums(profiles_to_test) > 0)]
    print(dim(profiles_to_test))
    
    print(profiles_to_test[1:10,1:10])
    
    classes_to_test <- fread(params$config_file, data.table = FALSE)$target_comps
    
    print(classes_to_test)
    
    distance <- params$distance
    
    ## Check that the distance is among our options
    
    if (! distance %in% c("bray", "jaccard")) {
      stop("Provided distance not in allowed options!")
    }
    
    
    curr_dist_profiles <- vegdist(t(profiles_to_test), method = distance)
    
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

        met_to_test <- as.data.frame(metadata_table)
        rownames(met_to_test) <- met_to_test$sample_id
        # head(met_to_test)

        met_to_test[which(met_to_test[, target_col] %in% class1), target_col] <- "class1"
        met_to_test[which(met_to_test[, target_col] %in% class2), target_col] <- "class2"

        met_to_test <- met_to_test[which(met_to_test[, target_col] %in% c("class1", "class2")), ]
        print(table(met_to_test[, target_col]))
        
        # print(head(profiles_to_test))
        
        print(profiles_to_test[1:3,1:3])
        print(rownames(met_to_test)[1:10])
        print(rownames(profiles_to_test)[1:10])

        met_to_test <- met_to_test[which(rownames(met_to_test) %in% rownames(profiles_to_test)), ]
        tmp_profiles_to_test <- as.data.frame(profiles_to_test[rownames(met_to_test), ])

        print(table(met_to_test[, target_col]))

        print(dim(met_to_test))
        print(dim(profiles_to_test))
        
        colnames(met_to_test)[which(colnames(met_to_test) == target_col)] <- "target_col"


        perm <- how(nperm = 99)
        setBlocks(perm) <- with(met_to_test, study_name)


        cat("\n\n --- --- --- ", class_to_compare, " --- --- --- \n")
        # set.seed(42)
        # permanova <- adonis2(profiles_to_test ~ target_col,
        #     data = met_to_test,
        #     permutations=perm,
        #     method = distance, parallel = params$ncores)
        # 
        # print(permanova)
        # 
        # tot_permanovas[[class_to_compare]] <- permanova
        
        ## Compute betadisper
        
        # print(as.data.frame(anova(betadisper(dis, groups))))
        
        break

    }
    
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
    # 
    # out_r2_pvals_permanova %>% write_tsv(paste0(params$output_folder, "out_permanova_analysis.tsv"))
    
}



main()
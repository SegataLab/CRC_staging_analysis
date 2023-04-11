library("data.table")
library("tidyverse")
library("Maaslin2")
library("optparse")

parser <- OptionParser()
parser <- add_option(parser, c("-i", "--input"), type="character", 
                     default="../data_March23/metaph4_formatted_asinsqrt.tsv")
parser <- add_option(parser, c("-n", "--name"), type="character", 
                     default="mpa4")
parser <- add_option(parser, c("-m", "--metadata"), type="character", 
                     default="../data_March23/metadata_formatted_Mar23.tsv")
parser <- add_option(parser, c("-c", "--comparisons"), type="character", 
                     default="study_name:stage_class:control__adenoma__0__I__II__III__IV;study_name:tumor_location_class_binary:right_sided__left_sided")
parser <- add_option(parser, c("-o", "--output_folder"), type="character", 
                     default="out_maaslin2")

arguments <- parse_args(parser)

# To run: OMP_NUM_THREADS=20 Rscript scr03_maaslin2.R

PREV_THR = 0.1
MIN_N_CLASS = 10

inp_path <- arguments$input
out_name <- arguments$name
inp_met_path <- arguments$metadata
comparisons <- unlist(strsplit(arguments$comparisons, split=";"))
out_dir <- arguments$output_folder

cat("\n\n----\n")
print(comparisons)
cat("----\n\n")
# Metaphlan-4
# out_name <- "mpa4"
# inp_path <- "../data_March23/metaph4_formatted_asinsqrt.tsv"
# inp_met_path <- "../data_March23/metadata_formatted_Mar23.tsv"

# Metaphlan-4 with other levels
# out_name <- "mpa4_wspecies_level"r
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

# out_dir <- "out_maaslin2/"

# comparisons <- c(
#     "study_name:stage_class:control__adenoma__0__I__II__III__IV",
#     "study_name:tumor_location_class_binary:right_sided__left_sided"
#     )

if (!dir.exists(out_dir)) {
    
    dir.create(out_dir)
    
}

# Maaslin2

inp_profiles <- fread(inp_path, data.table=FALSE)
metadata_table <- fread(inp_met_path, data.table=FALSE)

#########################################
# BLOCK 1
#########################################

for (comparison in comparisons) {
    
    study_col <- unlist(strsplit(comparison, split=":"))[1]
    tar_col <- unlist(strsplit(comparison, split=":"))[2]
    comb <- unlist(strsplit(comparison, split=":"))[3]
    
    target_classes <- unlist(strsplit(comb, split="__"))
    reference_class <- target_classes[1]
    
    print(reference_class)
#     stop("Here")
    
    n_metadata <- as.data.frame(metadata_table[which(
        metadata_table[, tar_col, drop=TRUE] %in% target_classes
    ), ])
    rownames(n_metadata) <- n_metadata$sample_id; n_metadata$sample_id <- NULL
    
#     if (tar_col == "tumor_location_class") {
        
#         n_metadata <- n_metadata[which(n_metadata$stage_class %in% c(
#             "0", "I", "II", "III", "IV", "crc_missing_stage"
#         )), ]
#     }
    
    n_profiles <- as.data.frame(inp_profiles)
    rownames(n_profiles) <- n_profiles$feat_name; n_profiles$feat_name <- NULL    
    n_profiles <- n_profiles[, intersect(rownames(n_metadata), colnames(n_profiles))]
    n_metadata <- n_metadata[intersect(rownames(n_metadata), colnames(n_profiles)), ]
    
    filt_metadata <- n_metadata[which(
            (!is.na(n_metadata$age)) & 
            (n_metadata$gender != "nd") &
            (!is.na(n_metadata$BMI))
        ), ]
    
    filt_profiles <- as.data.frame(inp_profiles)
    rownames(filt_profiles) <- filt_profiles$feat_name; filt_profiles$feat_name <- NULL
    filt_profiles <- filt_profiles[, intersect(rownames(filt_metadata), colnames(filt_profiles))]
    
    # Run Maaslin2 without correction
    
    n_out_dir <- paste0(out_dir, "out_maaslin2__", out_name, "__", tar_col, "__all/")
    
    n_metadata[, tar_col]  <- relevel(factor(
        n_metadata[, tar_col, drop=TRUE], levels = c(target_classes)), ref = reference_class)
    n_metadata[, study_col] <- as.factor(n_metadata[, study_col, drop=TRUE])
    
    print(unique(n_metadata[, tar_col]))
    print(unique(n_metadata[, study_col]))

    fit_data <- Maaslin2(
        n_profiles, n_metadata,
        n_out_dir, transform = "NONE",
        fixed_effects = c(tar_col),
        random_effects = c(study_col),
        normalization = 'NONE',
        standardize = FALSE, min_prevalence=0.1, cores=20,
        reference=paste0(tar_col, ",", reference_class)
    )
    
    write_tsv(n_metadata, paste0(n_out_dir, "inp_metadata.tsv"))
    write_tsv(n_profiles, paste0(n_out_dir, "inp_profiles.tsv"))
    
    # Run Maaslin2 with correction for age, gender and BMI
    
#     filt_metadata
#     filt_profiles
    
    n_out_dir <- paste0(out_dir, "out_maaslin2__", out_name, "__", tar_col, "__corrected/")
    
    filt_metadata[, tar_col]  <- relevel(factor(
        filt_metadata[, tar_col, drop=TRUE], levels = c(target_classes)), ref = reference_class)
    filt_metadata[, study_col] <- as.factor(filt_metadata[, study_col, drop=TRUE])
    filt_metadata$age <- as.numeric(filt_metadata$age)
    filt_metadata$gender <- as.factor(filt_metadata$gender)
    filt_metadata$BMI <- as.numeric(filt_metadata$BMI)

    fit_data <- suppressMessages(
        Maaslin2(
            filt_profiles, filt_metadata,
            n_out_dir, transform = "NONE",
            fixed_effects = c(tar_col),
            random_effects = c(study_col),
            normalization = 'NONE',
            standardize = FALSE, min_prevalence=0.1, cores=20))
    
    filt_metadata$sample_id <- rownames(filt_metadata)
    filt_profiles$feat_name <- rownames(filt_profiles)
    
    write_tsv(filt_metadata, paste0(n_out_dir, "inp_metadata.tsv"))
    write_tsv(filt_profiles, paste0(n_out_dir, "inp_profiles.tsv"))
    
}
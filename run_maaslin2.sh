mkdir -p out_maaslin2/

Rscript run_maaslin2.R \
    -i ../data_March23/metaph4_formatted_asinsqrt.tsv \
    -n mpa4 \
    -m ../data_March23/metadata_formatted_Mar23.tsv \
    -o out_maaslin2/ 2>&1 | tee out_maaslin2/out_maaslin2__mpa4.log;

Rscript run_maaslin2.R \
    -i ../data_March23/metaph4_formatted_asinsqrt_wspecies_level.tsv \
    -n mpa4_wspecies_level \
    -m ../data_March23/metadata_formatted_Mar23.tsv \
    -o out_maaslin2/ 2>&1 | tee out_maaslin2/out_maaslin2__mpa4_wspecies_level.log;

Rscript run_maaslin2.R \
    -i ../data_March23/hnn36_ecs_formatted_asinsqrt.tsv \
    -n hnn36_ec_table \
    -m ../data_March23/metadata_formatted_Mar23.tsv \
    -o out_maaslin2/ 2>&1 | tee out_maaslin2/out_maaslin2__hnn36_ec_table.log;

Rscript run_maaslin2.R \
    -i ../data_March23/hnn36_paths_formatted_asinsqrt.tsv \
    -n hnn36_pathways \
    -m ../data_March23/metadata_formatted_Mar23.tsv \
    -o out_maaslin2/ 2>&1 | tee out_maaslin2/out_maaslin2__hnn36_pathways.log;
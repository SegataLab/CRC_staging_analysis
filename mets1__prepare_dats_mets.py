import os
from copy import copy
import pandas as pd

import numpy as np

from multiprocessing import Pool

from itertools import product

import sys
import subprocess as sp

# python prepare_dats_mets.py ../data_March23/metadata_formatted_Feb24.tsv inp_dats.txt out_metas/ comparisons_metas.txt

# This script assumes that the minimum number of samples in at least one of the two categories to compare is 10
# and that the features is present in both cases and controls

# INP_METAD = "../data_March23/metadata_formatted_Mar23.tsv"
# INP_TABLE = "out_metas_divs/divs_formatted.tsv"
# OUT_DIR = "out_metas_divs/"
# TOT_COMPARISONS = "comparisons_metas.txt"

CURR_Rscript = "/shares/CIBIO-Storage/CM/scratch/users/gianmarco.piccinno/tools/conda_installation/envs/gpiccinno_oncobiome/bin/Rscript"
PAHT_run_metas = "/shares/CIBIO-Storage/CM/scratch/projects/gpiccinno_onc_crc_staging/20241017_revision_Piccinno_et_al/master_scripts/mets2__run_metas.R"


INP_METAD = sys.argv[1]
INP_DATS = sys.argv[2]
OR_OUT_DIR = sys.argv[3]
TOT_COMPARISONS = sys.argv[4]
MAX_POOL = 20

with open(INP_DATS, "r") as handle:
    tot_dats = [el.strip().split(":") for el in handle.readlines()]

if not os.path.isdir(OR_OUT_DIR):
    os.mkdir(OR_OUT_DIR)
    
def compute_stats(comb):

    comp = comb[0]
    study = comb[1]

    study_col = comp.split(":")[0]
    tar_col = comp.split(":")[1]
    comb = comp.split(":")[2]
    
    MIN_N_STUDES = comp.split(":")[3]

    vals = sum([el.split("__") for el in comb.split("_vs_")], [])

    class1 = comb.split("_vs_")[0]
    class2 = comb.split("_vs_")[1]

    tmp_df2 = copy(inp_metadata.loc[
        (inp_metadata[tar_col].isin(vals)) & (inp_metadata[study_col] == study)
    ])

    samples_class1 = list(tmp_df2.loc[tmp_df2[tar_col].isin(class1.split("__")), "sample_id"].values)
    samples_class2 = list(tmp_df2.loc[tmp_df2[tar_col].isin(class2.split("__")), "sample_id"].values)

    tmp_subdf1 = copy(inp_profiles[["feat_name"] + samples_class1])
    tmp_subdf1.set_index("feat_name", inplace=True)
    tmp_subdf2 = copy(inp_profiles[["feat_name"] + samples_class2])
    tmp_subdf2.set_index("feat_name", inplace=True)
    tot_subdf = copy(inp_profiles[["feat_name"] + samples_class1 + samples_class2])
    tot_subdf.set_index("feat_name", inplace=True)

    ################################################################################
    # Start:
    # Make the summaries
    ################################################################################

    summaries1 = tmp_subdf1.apply(lambda x: [
        np.mean(x),
        np.std(x, ddof = 1),
        sum([1 if el > 0 else 0 for el in x])/len(x),
        sum([1 if el > 0 else 0 for el in x]),
        len(x)
    ],
                                  axis = 1,
                                  result_type="expand").rename(
        columns = {
            0: "mean.c",
            1: "sd.c",
            2: "prev.c",
            3: "hits.c",
            4: "n.c"
        })

    summaries2 = tmp_subdf2.apply(lambda x: [
        np.mean(x),
        np.std(x, ddof = 1),
        sum([1 if el > 0 else 0 for el in x])/len(x),
        sum([1 if el > 0 else 0 for el in x]),
        len(x)
    ],
                                  axis = 1,
                                  result_type="expand").rename(
        columns = {
            0: "mean.e",
            1: "sd.e",
            2: "prev.e",
            3: "hits.e",
            4: "n.e"
        })

    summaries_tot = tot_subdf.apply(lambda x: [
        np.mean(x),
        np.std(x, ddof = 1),
        sum([1 if el > 0 else 0 for el in x])/len(x),
        sum([1 if el > 0 else 0 for el in x]),
        len(x)
    ],
                                  axis = 1,
                                  result_type="expand").rename(
        columns = {
            0: "mean.tot",
            1: "sd.tot",
            2: "prev.tot",
            3: "hits.tot",
            4: "n.tot"
        })

    ################################################################################
    # Start:
    # Make the summaries
    ################################################################################

    n_tmp_df = pd.concat([
        summaries1,
        summaries2,
        summaries_tot
    ], axis = 1)

    n_tmp_df["study_name"] = study
    n_tmp_df["comparison"] = comb

    return n_tmp_df

    
for dat in tot_dats:
    
    CURR_NAME = dat[0]
    INP_TABLE = dat[-1]
    
    OUT_DIR = os.path.join(OR_OUT_DIR, CURR_NAME)

    inp_profiles = pd.read_csv(INP_TABLE, sep="\t")
    inp_metadata = pd.read_csv(INP_METAD, sep="\t")

    if not os.path.isdir(OUT_DIR):
        os.mkdir(OUT_DIR)

    with open(TOT_COMPARISONS, "r") as handle:
        comparisons = [el.strip() for el in handle.readlines()]

    tot_combs = list(product(comparisons, inp_metadata["study_name"].unique()))


    tot_dfs = []

    with Pool(MAX_POOL) as p:
        tot_dfs = list(p.map(compute_stats, tot_combs))

    final_tot_dfs = pd.concat(tot_dfs)

    print(f"{final_tot_dfs.shape[0]:,}")
    final_tot_dfs.dropna(inplace=True)
    final_tot_dfs.reset_index(inplace=True)
    print(f"{final_tot_dfs.shape[0]:,}")

    final_tot_dfs.to_csv(os.path.join(OUT_DIR, "dat_met.tsv"), sep="\t", index=False)

    CMD = [CURR_Rscript, PAHT_run_metas, os.path.join(OUT_DIR, "dat_met.tsv"), OUT_DIR+"/", TOT_COMPARISONS]

    print(CMD)
    sp.run(CMD)
    print("Here")
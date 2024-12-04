from copy import copy
import pandas as pd
import numpy as np

import os
import re

from itertools import combinations, product

import time
from multiprocessing import Pool

import argparse as ap
import json

# Example to run
# python step4_prepare_profiles4classification.py --input_config config_prepare4classification.json --output_folder inp_classification/

def get_time():
    return datetime.datetime.now().strftime("%I:%M:%m%p on %B %d, %Y")

def printt(*args):
    print(get_time(), " --", end = "") # , args
    
    for arg in args:
        sys.stdout.write(" ")
        sys.stdout.write(repr(arg))
        
    sys.stdout.write("\n")
    
    return

def read_params():
    p = ap.ArgumentParser(description=(""),
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)
    
    p.add_argument('-i', '--input_config', type=str, default=None, help="")
    p.add_argument('-o', '--output_folder', type=str, default=None, help=(""))
    
    return p.parse_args()


args = read_params()

with open(args.input_config, "r") as handle:
    inputs = json.load(handle)
    
comparisons = inputs["comparisons"]

inp_tables = {}

for key in inputs["inp_tables"].keys():
    
    if os.path.isfile(os.path.join(args.output_folder.split("/")[0], inputs["inp_tables"][key])):
        inp_tables[key] = os.path.join(args.output_folder.split("/")[0], inputs["inp_tables"][key])
    else:
        inp_tables[key] = inputs["inp_tables"][key]
    

# inp_tables = {key: os.path.join(args.output_folder.split("/")[0], inputs["inp_tables"][key]) for key in inputs["inp_tables"] }

path_metadata = inputs["path_metadata"]
out_path = args.output_folder

if not os.path.isdir(out_path):
    os.mkdir(out_path)
    

first_col_inp_tables = "feat_name"
first_col_met_table = "sample_id"
study_col_met_table = "study_name"
l_thr = 15


metadata = pd.read_csv(path_metadata, sep="\t")
metadata.fillna("", inplace=True)

plain_comparisons = sum([[(el, val) for val in comparisons[el]] for el in comparisons.keys()], [])

crosspred_combs = sum([[
    ("crosspred", tuple([el[0]]), tuple([el[1]])), 
    ("crosspred", tuple([el[1]]), tuple([el[0]]))] 
    for el in combinations(set(metadata[study_col_met_table].values), 2)], [])
print(len(crosspred_combs))

lodo_combs = [("lodo", tuple(set(metadata[study_col_met_table].values) - set([el])), tuple([el])) 
              for el in set(metadata[study_col_met_table].values)]
print(len(lodo_combs))

cvperdat_combs = [("cvperdat", tuple([el])) for el in set(metadata[study_col_met_table].values)]
print(len(cvperdat_combs))

tot_tasks = crosspred_combs + lodo_combs + cvperdat_combs

# Make total tasks for metadata
tot_tasks_combs = list(product(
    inp_tables,
    plain_comparisons,
    tot_tasks))

max_len = len(str(len(tot_tasks_combs)))
print(max_len)

tot_tasks_combs_codes = {j: {
    "taskn": j,
    "inp_table": inp_tables[el[0]],
    "comparison": el[1][1],
    "class1": el[1][1].split("_vs_")[0].split("__"),
    "class2": el[1][1].split("_vs_")[1].split("__"),
    "task_descr": el[2][0],
    "set1": el[2][1],
    "set2": el[2][-1],
    "prefix": el[0],
    "target_col": el[1][0]
    
} for j, el in enumerate(tot_tasks_combs)}
print(len(tot_tasks_combs))
print(len(tot_tasks_combs_codes))

out_summary = pd.DataFrame(list(set([(el[0], el[1][0], el[1][1], el[2][0]) for el in tot_tasks_combs]))).rename(
    columns = {0: "target_profile", 1: "target_col", 2: "target_comparison", 3: "classification_type"}
)

out_summary.to_csv(f"{out_path}/out_summary_of_comparisons.tsv", sep="\t", index=False)

# display(out_summary.head())
# display(out_summary.shape)

def make_classes(taskn):
    
    tmp_met = copy(metadata)
    descr_table = copy(tot_tasks_combs_codes[taskn])
    
    task = tot_tasks_combs_codes[taskn]
    
    inp_table = task["inp_table"]
    comparison = task["comparison"]
    class1 = task["class1"]
    class2 = task["class2"]
    task_descr = task["task_descr"]
    set1 = task["set1"]
    set2 = task["set2"]
    
    prefix = task["prefix"]
    target_col = task["target_col"]
    
    tmp_met[f"task{taskn}_toconsider"] = "nd"
    tmp_met[f"task{taskn}_classes"] = "nd"
    
    samples_set1_class1 = set(tmp_met.loc[((tmp_met[study_col_met_table].isin(set1)) & (tmp_met[target_col].isin(class1))), first_col_met_table].values)
    samples_set1_class2 = set(tmp_met.loc[((tmp_met[study_col_met_table].isin(set1)) & (tmp_met[target_col].isin(class2))), first_col_met_table].values)
    samples_set2_class1 = set(tmp_met.loc[((tmp_met[study_col_met_table].isin(set2)) & (tmp_met[target_col].isin(class1))), first_col_met_table].values)
    samples_set2_class2 = set(tmp_met.loc[((tmp_met[study_col_met_table].isin(set2)) & (tmp_met[target_col].isin(class2))), first_col_met_table].values)
    
    tot_samples_class1 = samples_set1_class1 | samples_set2_class1
    tot_samples_class2 = samples_set1_class2 | samples_set2_class2
    tot_samples_classes = tot_samples_class1 | tot_samples_class2
    
    comp_to_consider = False
    if ((len(samples_set1_class1) >= l_thr) & (len(samples_set1_class2) >= l_thr) & (len(samples_set2_class1) >= l_thr) & (len(samples_set2_class2) >= l_thr)):
        comp_to_consider = True
        
    descr_table["comp_to_consider"] = comp_to_consider
    
    descr_table["n_samples_set1_class1"] = len(samples_set1_class1)
    descr_table["n_samples_set1_class2"] = len(samples_set1_class2)
    descr_table["n_samples_set2_class1"] = len(samples_set2_class1)
    descr_table["n_samples_set2_class2"] = len(samples_set2_class2)
        
        
    tmp_met.loc[tmp_met[first_col_met_table].isin(tot_samples_classes), f"task{taskn}_toconsider"] = "yes"
    
    tmp_met.loc[tmp_met[first_col_met_table].isin(tot_samples_class1), f"task{taskn}_classes"] = "class1"
    tmp_met.loc[tmp_met[first_col_met_table].isin(tot_samples_class2), f"task{taskn}_classes"] = "class2"
    
    if set1 != set2:
        tar = set2
    else:
        tar = ""
    
    descr_table_df = pd.DataFrame.from_dict({
        "taskn": [f"task{taskn}"],
        "n_samples_set1_class1": len(samples_set1_class1),
        "n_samples_set1_class2": len(samples_set1_class2),
        "n_samples_set2_class1": len(samples_set2_class1),
        "n_samples_set2_class2": len(samples_set2_class2),
        "comp_to_consider": comp_to_consider,
        "inp_table_path": os.path.join(out_path, prefix + ".tsv"),
        "tar": tar
    })
    
    tmp_met = copy(tmp_met[[first_col_met_table, f"task{taskn}_toconsider", f"task{taskn}_classes"]])
    tmp_met.set_index(first_col_met_table, inplace=True)
    tmp_met = copy(tmp_met.T)
    
    return [descr_table_df, tmp_met]

max_pool = 20

start = time.time()

with Pool(max_pool) as p:
    pool_outputs = list(p.imap(make_classes, list(tot_tasks_combs_codes.keys())))
    
tot_out_descrs = pd.concat([el[0] for el in pool_outputs])
tot_out_tasks = pd.concat([el[1] for el in pool_outputs])

end = time.time()
print(end - start)

descr_table_df = pd.DataFrame.from_dict(tot_tasks_combs_codes, orient="index")
descr_table_df["taskn"] = "task" + descr_table_df["taskn"].astype(str)

n_tot_out_descrs = copy(descr_table_df.merge(tot_out_descrs, on="taskn", how="inner"))

n_tot_out_descrs.to_csv(os.path.join(out_path, "tasks_descr.tsv"), sep="\t", index=False)

tot_tasks = copy(n_tot_out_descrs.loc[n_tot_out_descrs["comp_to_consider"] == True,
                                      ["task_descr",
                                       "taskn",
                                       "prefix",
                                       "inp_table_path",
                                       "tar"
                                      ]])

tot_tasks.to_csv(os.path.join(out_path, "tasks.tsv"), sep="\t", index=False)
tot_tasks.to_csv(os.path.join(out_path, "tasks.txt"), sep=":", index=False, header=False)

## Make the feature tables

tasks_profiles = [{
    "prefix": el,
    "inp_table": inp_tables[el]
}
    for el in inp_tables.keys()]
print(len(tasks_profiles))

tot_out_tasks["metID"] = tot_out_tasks.index
tot_out_tasks.rename(columns={"sample_id": "metID"}, inplace=True)
tot_out_tasks = copy(tot_out_tasks[["metID"] + [col for col in tot_out_tasks.columns if col != "metID"]])

for task in tasks_profiles:
    # print(task)
    
    prefix = task["prefix"]
    print("\n---", prefix)
    
    inp_table = task["inp_table"]
    
    path_inp_table = inp_table
    
    curr_table = pd.read_csv(path_inp_table, sep="\t").rename(columns={first_col_inp_tables: "metID"})
    curr_table["metID"] = prefix + curr_table["metID"]
        
    tmp_met = copy(metadata)[[first_col_met_table, study_col_met_table]].dropna(axis=0, how="any")            
    tmp_met.set_index(first_col_met_table, inplace=True)
    tmp_met = copy(tmp_met).T
    tmp_met["metID"] = tmp_met.index
    tmp_met.loc[tmp_met["metID"] != study_col_met_table, "metID"] = prefix + tmp_met.loc[
        tmp_met["metID"] != study_col_met_table, "metID"]    
    
    tot_cols = ["metID"] + [el for el in set(tot_out_tasks.columns) & set(tmp_met.columns) & set(curr_table.columns) if el != "metID"]
    
    tmp_met = copy(tmp_met[tot_cols])
    
    
    n_curr_table = pd.concat([tot_out_tasks[tot_cols], tmp_met[tot_cols], curr_table[tot_cols]])
    
    print( os.path.join(out_path, prefix + ".tsv") )
    n_curr_table.to_csv(os.path.join(out_path, prefix + ".tsv"), sep="\t", index=False)
    
print("Here")
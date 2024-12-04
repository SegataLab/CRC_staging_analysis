import os
from copy import copy
import pandas as pd

import numpy as np

import argparse as ap
import sys

import time
import datetime

# Example
# python step6_organize_AUCs_metaml.py -i curr_tmp_output_folder_03/inp_classification/ -l curr_tmp_output_folder_03/out_classification/ -m ../../data_March23/metadata_formatted_Feb24.tsv -t curr_tmp_output_folder_03/inp_classification/tasks_descr.tsv -o curr_tmp_output_folder_03

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
    
    p.add_argument('-i', '--inp_classification', type=str, default=None, help="")
    p.add_argument('-l', '--out_classification', type=str, default=None, help="")
    p.add_argument('-m', '--metadata', type=str, default=None, help="")
    p.add_argument('-t', '--tasks_descr', type=str, default=None, help="")
    p.add_argument('-o', '--output_folder', type=str, default=None, help=(""))
    
    return p.parse_args()


def main():
    args = read_params()
    
    printt(args)
    
    metadata = pd.read_csv(args.metadata, sep="\t")
    table_tasks = pd.read_csv(args.tasks_descr, sep="\t")

    table_tasks["comparison"].unique()

    out_prefix = "out_metaml_rf_"
    out_folders = {
        "lodo": os.path.join(args.out_classification, "out_ml_lodo/"),
        "cvperdat": os.path.join(args.out_classification, "out_ml_cvperdat/"),
        "crosspred": os.path.join(args.out_classification, "out_ml_crosspred/")
    }
    
    class_to_consider = [
        key for key in out_folders.keys() if len(out_folders[key]) > 0
    ]
    
    
    tables = {
        table: table_tasks.loc[
            ((table_tasks["task_descr"] == table) & (table_tasks["comp_to_consider"] == True))]
        for table in class_to_consider
    }

    out_aucs = {
        "taskn": [],
        "auc": [],
        "sds": []
    }

    missing = []
    empty = []

    for key in tables.keys():
        tasks_l = tables[key]["taskn"].values

        for task in tasks_l:

            tar_path = os.path.join(out_folders[key], f"{out_prefix}{task}.txt")

            if not os.path.exists(tar_path):

                # print("MISSING:", key, task, tar_path)
                missing.append((key, task, tar_path))
                continue
                # raise Exception("Here2")

            with open(tar_path, "r") as handle:
                tmp_lines_ = handle.readlines()

            if len(tmp_lines_) > 1:


                tot_samples = int(float(tmp_lines_[0].split("\t")[1].strip()))

                mean_auc, std_auc = [el for el in tmp_lines_ 
                                     if "auc" in el][0].strip().split("\t")[1:]

                out_aucs["taskn"].append(task)
                out_aucs["auc"].append(float(mean_auc))
                out_aucs["sds"].append(float(std_auc))

            else:

                # print("EMPTY:", key, task, tar_path)
                empty.append((key, task, tar_path))
                continue

    print("N missing classification tasks:", len(missing), "; N empty classification tasks", len(empty))
    df_out_aucs = pd.DataFrame.from_dict(out_aucs)

    n_tables = {}

    # print(" --->>>>>")
    # print(tables.keys())
    # print(" --->>>>>")
    
    for key in tables.keys():
        n_tables[key] = tables[key].merge(df_out_aucs, on="taskn", how="left")

        n_tables[key].to_csv(
            os.path.join(args.out_classification, 
                         f"table_of_aucs_rf_bench_{key}.tsv"),
                            sep="\t", index=False)

        # print(key)
        
        aucs = copy(n_tables[key])
        
        # print(" --->>>>>")
        # print(aucs["inp_table"].unique())
        # print(" --->>>>>")
        
        for target_col in aucs["target_col"].unique():
            for prefix in aucs["prefix"].unique():
            
                for comp in aucs["comparison"].unique():

                    tot_summaries = []

                    for inp_table in aucs["inp_table"].unique():

                        tmp_to_print = copy(aucs.loc[
                                            (aucs["inp_table"] == inp_table) & 
                                            (aucs["comparison"] == comp) & 
                                            (aucs["target_col"] == target_col) & 
                                            (aucs["prefix"] == prefix) & 
                                            (~aucs["auc"].isnull()),
                                            ["tar", "auc"]
                                        ].set_index("tar").T) #





                        if tmp_to_print.shape[1] > 0:

                            if "PRESCIENT" in tmp_to_print.columns:
                                tmp_to_print.drop(columns = ["PRESCIENT"], inplace = True)
                            tmp_to_print["average"] = tmp_to_print.apply(lambda x: np.mean(x), axis = 1)
                            tmp_to_print["inp_table"] = inp_table

                            print("\n --- --- --- --- --- --- --- --- --- --- --- --- ---\n", key, inp_table, comp)
                            print(tmp_to_print)

                            tot_summaries.append(tmp_to_print)


                    if len(tot_summaries) > 0:

                        tot_summaries_df = pd.concat(tot_summaries)

                        tot_summaries_df.to_csv(
                            os.path.join(args.out_classification, 
                                         f"table_of_aucs_rf_bench_{key}__{target_col}__{comp}__{prefix}.tsv"),
                            sep="\t", index=False
                        )
                
                


#         print("Here")

#         aucs = copy(n_tables[key])
#         print(aucs.head())


#         aucs["inp_table"].unique()


#         tot_summaries = []

#         for inp_table in aucs["inp_table"].unique():

#             for comp in aucs["comparison"].unique():


#                 tmp_df = copy(aucs.loc[
#                     (aucs["inp_table"] == inp_table) & 
#                     (aucs["comparison"] == comp) & 
#                     (aucs["tar"] != "PRESCIENT"), ['inp_table', 'comparison', 'tar', 'auc']])
#                 tmp_df = pd.concat([tmp_df, pd.DataFrame([tuple([inp_table, comp, 'Average', aucs.loc[
#                     (aucs["inp_table"] == inp_table) & 
#                     (aucs["comparison"] == comp) & 
#                     (aucs["tar"] != "PRESCIENT"), "auc"].astype(float).mean()])]).rename(columns = {0: "inp_table", 1: "comparison", 2: "tar", 3: "auc"})])

#                 # tmp_df['inp_table'] = inp_table

#                 tot_summaries.append(tmp_df)


#                 print('\n\n ----', inp_table, comp, aucs.loc[
#                     (aucs["inp_table"] == inp_table) & 
#                     (aucs["comparison"] == comp) & 
#                     (aucs["tar"] != "PRESCIENT"), "auc"].astype(float).mean())
#                 print(aucs.loc[
#                     (aucs["inp_table"] == inp_table) & 
#                     (aucs["comparison"] == comp) & 
#                     (aucs["tar"] != "PRESCIENT"), 
#                     ["comparison", "tar", "auc"]].sort_values(by = "auc"))


#         tot_summaries_df = pd.concat(tot_summaries)

    
if __name__ == '__main__':
    
    print(get_time(), " -- Starting")
    
    t0 = time.time()
    main()
    t1 = time.time()
    print(get_time(), " -- Total elapsed time {}s".format(int(t1 - t0)))
    
    
    print(get_time(), " -- Concluding")
    
    sys.exit(0)






#!/bin/bash

# parallel -j10 "bash run_metaml.sh {}" ::: `cat inp_classification/tasks.txt | grep lodo`
# parallel -j10 "bash run_metaml.sh {}" ::: `cat inp_classification/tasks.txt | grep crosspred`
# parallel -j10 "bash run_metaml.sh {}" ::: `cat inp_classification/tasks.txt | grep cvperdat`
# parallel -j10 "bash run_metaml.sh {}" ::: `cat inp_classification/tasks_filtered.txt`
# bash run_metaml.sh lodo:task240:mpa4__:inp_classification/mpa4__.tsv:VogtmannE_2016
# bash run_metaml.sh crosspred:task0:mpa4__:inp_classification/mpa4__.tsv:GuptaA_2019

# bash run_metaml.sh cvperdat:task256:mpa4__:inp_classification/mpa4__.tsv:

line=$1

class_script_metaml="classification_oncCRCstaging.py"
out_name="out_metaml_rf_";

type_of_task="$(echo ${line} | cut -d':' -f1)";
task="$(echo ${line} | cut -d':' -f2)";
profile="$(echo ${line} | cut -d':' -f3)";
curr_filename="$(echo ${line} | cut -d':' -f4)";
tar="$(echo ${line} | cut -d':' -f5)";

out_path="out_classification/out_ml_${type_of_task}/"
mkdir -p ${out_path}

echo "${task}, ${tar}, ${profile}, ${curr_filename}";

if [[ ${type_of_task} == "lodo" ]]; then

    python3 ${class_script_metaml} ${curr_filename} ${out_path}/${out_name}${task} \
            -d 1:${task}_classes:class2 \
            -sl ${task}_toconsider:yes \
            -l rf \
            -t study_name:${tar} \
            -nt 1000 \
            -nsl 5 \
            -mf auto \
            -nc 10 \
            -hv 0 \
            -r 1 \
            -p 1 \
            --no_norm \
            -z ${profile} -w 2>&1 | tee ${out_path}/${out_name}${task}.log;

fi

if [[ ${type_of_task} == "crosspred" ]]; then

    python3 ${class_script_metaml} ${curr_filename} ${out_path}/${out_name}${task} \
            -d 1:${task}_classes:class2 \
            -sl ${task}_toconsider:yes \
            -l rf \
            -t study_name:${tar} \
            -nt 1000 \
            -nsl 5 \
            -mf auto \
            -nc 10 \
            -hv 0 \
            -r 1 \
            -p 1 \
            -df \
            --no_norm \
            -z ${profile} -w 2>&1 | tee ${out_path}/${out_name}${task}.log;

fi

if [[ ${type_of_task} == "cvperdat" ]]; then

    python3 ${class_script_metaml} ${curr_filename} ${out_path}/${out_name}${task} \
            -d 1:${task}_classes:class2 \
            -sl ${task}_toconsider:yes \
            -l rf \
            -nt 1000 \
            -nsl 5 \
            -mf auto \
            -nc 10 \
            -hv 0 \
            -r 20 \
            -p 10 \
            -df \
            --no_norm \
            -z ${profile} -w 2>&1 | tee ${out_path}/${out_name}${task}.log;

fi
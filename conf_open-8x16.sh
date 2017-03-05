#!/bin/bash

#input deck for the analysis of the KPI-DBW2.

bc_name="open"
ens_name="open-8x16"
ens_dir="/home/gregm/DBW2/open-8x16/results/alg_wflow"

#
num_thermal="10000"

# time extent size
T="16"

# MD spacing between two consecutive measurements
meas_spacing_MD="10"

# block size in MD unit
block_size_MD="400"

# 
max_MD_sep=$[$block_size_MD/2]
#echo $max_MD_sep

#
err_max_MD_sep=$[$block_size_MD/4]
#echo $err_max_MD_sep

# number of replicas
num_reps="1"

# list of replicas
rep_names="open-8x16"

tau_guess="18"
D_guess="0.115"


#!/bin/bash

#input deck for the analysis of the KPI-DBW2.

# boundary condition
bc_name="periodic"

# ensemble name
ens_name="l2464f211b600m0102m0509m635a"

# directory of the wilson flow files
ens_dir="/Users/hummingtree/alice/DBW2-KPI/wflow_tcharge_cps/results/l2464f211b600m0102m0509m635a"

# number of thermalization steps, i.e. skip these configurations
num_thermal="284"

# time extent size
T="64"

# MD spacing between two consecutive measurements
meas_spacing_MD="5"

# block size in MD unit
block_size_MD="50"

# the maximum MD separation when computing correlations
#max_MD_sep=$[$block_size_MD]
max_MD_sep="25"
#echo $max_MD_sep

# the maximum MD separation when computing error of the correlations
# err_max_MD_sep=$[$block_size_MD/4]
err_max_MD_sep=15
#echo $err_max_MD_sep

# number of replicas
num_reps="1"

# list of replicas
rep_names="l2464f211b600m0102m0509m635a"

# initial values of the parameters for the fit
tau_guess="20"
D_guess="0.2"

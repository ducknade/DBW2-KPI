#!/bin/bash

#input deck for the analysis of the KPI-DBW2.

# boundary condition
bc_name="open"

# ensemble name
ens_name="open-10x20"

# directory of the wilson flow files
ens_dir="/home/gregm/DBW2/open-10x20/results/alg_wflow"

# number of thermalization steps, i.e. skip these configurations
num_thermal="10000"

# time extent size
T="20"

# MD spacing between two consecutive measurements
meas_spacing_MD="15"

# block size in MD unit
block_size_MD="1200"

# the maximum MD separation when computing correlations
max_MD_sep=$[$block_size_MD/2]
#echo $max_MD_sep

# the maximum MD separation when computing error of the correlations
err_max_MD_sep=$[$block_size_MD/4]
#echo $err_max_MD_sep

# number of replicas
num_reps="1"

# list of replicas
rep_names="open-10x20"

# initial values of the parameters for the fit
tau_guess="55"
D_guess="0.115"


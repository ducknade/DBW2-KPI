#!/bin/bash

#input deck for the analysis of the KPI-DBW2.

# boundary condition
bc_name="periodic"

# ensemble name
ens_name="L32_beta6.4_4.0_5"

# directory of the wilson flow files
ens_dir="/Users/hummingtree/alice/DBW2-KPI/wflow_tcharge_cps/results/$ens_name"

# number of thermalization steps, i.e. skip these configurations
num_thermal="3000"

# time extent size
T="32"

# MD spacing between two consecutive measurements
meas_spacing_MD="10"

# block size in MD unit
block_size_MD="1000"

# the maximum MD separation when computing correlations
#max_MD_sep=$[$block_size_MD]
max_MD_sep="500"
#echo $max_MD_sep

# the maximum MD separation when computing error of the correlations
# err_max_MD_sep=$[$block_size_MD/4]
err_max_MD_sep=250
#echo $err_max_MD_sep

# number of replicas
num_reps="1"

# list of replicas
rep_names="L32_beta6.4_4.0_5"

# initial values of the parameters for the fit
tau_guess="12.13"
D_guess="0.0635"

wfs_ref_min="0.0"
wfs_ref_max="0.30"
wfs_ref_add="0.05"


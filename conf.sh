#!/bin/bash

#input deck for the analysis of the KPI-DBW2.

# boundary condition
bc_name="periodic"

# ensemble name
ens_name="periodic-10x20"

# directory of the wilson flow files
ens_dir="/Users/hummingtree/alice/DBW2-KPI/wflow_tcharge_cps/results/periodic-10x20"

# number of thermalization steps, i.e. skip these configurations
num_thermal="300"

# time extent size
T="20"

# MD spacing between two consecutive measurements
meas_spacing_MD="15"

# block size in MD unit
block_size_MD="1500"

# the maximum MD separation when computing correlations
#max_MD_sep=$[$block_size_MD]
max_MD_sep="750"
#echo $max_MD_sep

# the maximum MD separation when computing error of the correlations
# err_max_MD_sep=$[$block_size_MD/4]
err_max_MD_sep=300
#echo $err_max_MD_sep

# number of replicas
num_reps="1"

# list of replicas
rep_names="periodic-10x20"

# initial values of the parameters for the fit
tau_guess="12.13"
D_guess="0.0635"

wfs_ref_min="0.0"
wfs_ref_max="0.25"
wfs_ref_add="0.05"


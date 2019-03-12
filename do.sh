#!/bin/bash

source ./conf.sh

cp conf.sh conf_${ens_name}.sh

# step 1: from raw cps data to Qslice, i.e. topological change density summed on every time slice
dir_1=$(pwd)"/cps_to_Qslice"

# step 2: from Qslice to correlation
dir_2=$(pwd)"/Qslice_to_correlation"

# step 3: from correlation to fit
dir_3=$(pwd)"/correlation_to_fit"

# execute step 1
cd $dir_1
./get_Qslice_histories.py "--"$bc_name $ens_name $num_thermal $ens_dir ${wfs_ref_min} ${wfs_ref_max} ${wfs_ref_add}

# execute step 2
cd $dir_2
make
wfs=${wfs_ref_min}
for wfs in `seq -f "%.4f" ${wfs_ref_min} ${wfs_ref_add} ${wfs_ref_max}`; do
  echo ${wfs}
./corr_jackknifes.x $bc_name ${ens_name}_wfs${wfs} $T $meas_spacing_MD $block_size_MD $max_MD_sep $err_max_MD_sep $num_reps ${rep_names}_wfs${wfs} | tee tmp.crap
Njack=$(grep "Num jackknife blocks per replica =" tmp.crap | cut -d "=" -f 2)
rm tmp.crap
done

cd ../fft
./fft.py ${ens_name} ${Njack} ${T} ${meas_spacing_MD}  
#./fft.py ${ens_name} 15 ${T} ${meas_spacing_MD}  

# execute step 3
#cd $dir_3
#matlab -nodesktop -r "DoFitJackknife($tau_guess, $D_guess, '$ens_name', '$bc_name', $T, $[$Njack+1]); quit;"
#
#cat $dir_3/results/$ens_name*

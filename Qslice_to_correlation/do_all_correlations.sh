#!/bin/bash

set -ex

# $ ./corr_jackknifes
# Usage: ./corr_jackknifes bc ens_name T meas_spacing_MD block_size_MD max_md_sep err_max_md_sep num_reps <rep names>

#./corr_jackknifes     open                open-8x16 16 10 400 200 100 1                open-8x16
#./corr_jackknifes periodic            periodic-8x16 16 10 400 200 100 1            periodic-8x16
#./corr_jackknifes periodic  shorttraj-periodic-8x16 16 10 400 200 100 1  shorttraj-periodic-8x16
#./corr_jackknifes periodic   longtraj-periodic-8x16 16 10 400 200 100 1   longtraj-periodic-8x16
./corr_jackknifes periodic  longtraj2-periodic-8x16 16 10 400 200 100 1  longtraj2-periodic-8x16

#./corr_jackknifes     open     open-10x20 20 15 1200 600 300 1     open-10x20
#./corr_jackknifes periodic periodic-10x20 20 15 1200 600 300 1 periodic-10x20
#
#./corr_jackknifes     open     open-12x24 24 21 4000 2000 1000 1     open-12x24
#./corr_jackknifes periodic periodic-12x24 24 21 4000 2000 1000 1 periodic-12x24
#
#./corr_jackknifes     open     open-14x28 28 28 12000 6000 3000 4     open-14x28-1     open-14x28-2     open-14x28-3     open-14x28-4
#./corr_jackknifes periodic periodic-14x28 28 28 12000 6000 3000 4 periodic-14x28-1 periodic-14x28-2 periodic-14x28-3 periodic-14x28-4
#
#./corr_jackknifes     open     open-16x32 32 40 40000 20000 10000 1     open-16x32
#./corr_jackknifes periodic periodic-16x32 32 40 40000 20000 10000 1 periodic-16x32
#
#
#./corr_jackknifes periodic iwasaki-closed-2 64 10 1000 1000 500 1 iwasaki-closed-2



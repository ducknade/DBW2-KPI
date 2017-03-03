#!/bin/bash

set -ev

echo "SKIPPING SOME ENSEMBLES"

#./get_Qslice_histories.py --periodic periodic-8x16 100
#./get_Qslice_histories.py --open         open-8x16 100

#./get_Qslice_histories.py --periodic periodic-10x20 300
#./get_Qslice_histories.py --open         open-10x20 300

#./get_Qslice_histories.py --periodic periodic-12x24 1000
#./get_Qslice_histories.py --open         open-12x24 1000

#./get_Qslice_histories.py --periodic periodic-14x28-1 3000
#./get_Qslice_histories.py --periodic periodic-14x28-2 3000
#./get_Qslice_histories.py --periodic periodic-14x28-3 3000
#./get_Qslice_histories.py --periodic periodic-14x28-4 3000
#./get_Qslice_histories.py --open         open-14x28-1 3000
#./get_Qslice_histories.py --open         open-14x28-2 3000
#./get_Qslice_histories.py --open         open-14x28-3 3000
#./get_Qslice_histories.py --open         open-14x28-4 3000

#./get_Qslice_histories.py --periodic periodic-16x32 10000
#./get_Qslice_histories.py --open         open-16x32 10000

#./get_Qslice_histories.py --periodic shorttraj-periodic-8x16 400
#./get_Qslice_histories.py --periodic  longtraj-periodic-8x16  20
./get_Qslice_histories.py --periodic  longtraj2-periodic-8x16  20

#./get_iwasaki_Qslice.py iwasaki-closed-2

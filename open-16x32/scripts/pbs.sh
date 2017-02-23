#!/bin/bash

#Instructions for the queueing system:
#PBS -q cuth
#PBS -r n
#PBS -l nodes=1,walltime=99:00:00
#PBS -N open-16x32

${PBS_O_WORKDIR}/../../../global-scripts/run.sh

#!/bin/bash

#SBATCH -p pdlabs
#SBATCH --time=5:00

declare -i n=$1
declare -i d=$2
declare -i k=$3

srun ./test_sequential $n $d $k

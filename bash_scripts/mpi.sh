#!/bin/bash

#SBATCH -p pdlabs
#SBATCH --time=0:05:00

declare -r exec=$1
declare -i n=$2
declare -i d=$3
declare -i k=$4

if [ "$exec" == "sync" ]
then
    srun ./test_synchronous $n $d $k $exec.txt
elif [ "$exec" == "async" ]
then
    srun ./test_asynchronous $n $d $k $exec.txt
else
    echo First argument must be "sync" or "async"
fi

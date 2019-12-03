#!/bin/bash

module load gcc openmpi netlib-lapack

declare -r exec=$1

if [ "$exec" = "sync" ]
then
	make test_synchronous
elif [ "$exec" == "async" ]
then
	make test_asynchronous
else
	echo First argument must be "sync" or "async"
fi

declare -i min_p=8
declare -i max_p=32
declare -i inc_p=$min_p

declare -i min_n=8000
declare -i max_n=40000
declare -i inc_n=$min_n

declare -i min_d=50
declare -i max_d=100
declare -i inc_d=$min_d

declare -i min_k=10
declare -i max_k=20
declare -i inc_k=$min_k

declare -i n
declare -i d
declare -i k
declare -i nodes
declare -i cores=8
declare -i p

for ((pc=$min_p; pc<=$max_p; pc += $inc_p))
do
	p=pc
	echo $p
	if (( min_n % p != 0 ))
	then
		continue
	fi
	nodes=p/$cores
    for ((nc=$min_n; nc<=$max_n; nc += $inc_n))
    do
    	n=nc
    	for ((dc=$min_d; dc<=$max_d; dc += $inc_d))
    	do
      		d=dc
        	for ((kc=$min_k; kc<=$max_k; kc += $inc_k))
        	do
        		k=kc
        		sbatch --nodes=$nodes --ntasks-per-node=$cores mpi.sh $exec $n $d $k
        	done
    	done
    done
done

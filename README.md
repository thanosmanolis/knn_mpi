# k-Nearest Neighbors

## Description
Sequential and parallel implementation (with MPI) of k-NN.

## How to run it

1. If your gcc version is previous than gcc-7, then change it in the Makefile
2. Run "make" to compile all implementations, or "make test_<__exec_name__>" to compile
only one implementation
3. Run the executable with:
    1. ./test_sequential "arg1" "arg2" "arg3"
    2. mpirun -np __number_of_processes__ ./test_synchronous __arg1__ __arg2__ __arg3__
    3. mpirun -np __number_of_processes__ ./test_asynchronous __arg1__ __arg2__ __arg3__

The three arguments, are the n (data), d (dimensions), k (kNN-neighbors) values accordingly. If no arguments are included at the run command, then the executable will run with default values.

# k-Nearest Neighbors

## Description
Sequential and parallel implementation (with MPI) of k-NN.

## How to run it

1. If your gcc version is previous than gcc-7, then change it in the Makefile
2. Type ``make`` to compile all implementations, or ``make test_<exec_name>`` to compile
only one implementation
3. Execution:
    1. Sequential: ``./test_sequential arg1 arg2 arg3``
    2. Synchronous: ``mpirun -np <number_of_processes> ./test_synchronous arg1 arg2 arg3``
    3. Asynchronous: ``mpirun -np <number_of_processes> ./test_asynchronous arg1 arg2 arg3``

The three arguments, are the n (number of data per-process), d (dimensions), k (kNN-neighbors) values accordingly. If no arguments are included at the run command, then the executable will run with default values.

# k-Nearest Neighbors

## Description
Sequential and parallel implementation (with MPI) of k-NN.

## How to run it

1. If your gcc version is previous than gcc-7, then change it in the Makefile
2. The main() function for testing is located into the tester files
3. Run "make" to compile all implementations, or "make test_*exec_name*" to compile
only one implementation
4. Run the executable with:
    1. ./test_sequential "arg1" "arg2" "arg3"
    2. ./test_synchronous "arg1" "arg2" "arg3"
    3. ./test_asynchronous "arg1" "arg2" "arg3"

The three arguments, are the n (data), d (dimensions), k (kNN-neighbors) values accordingly. If no arguments are included at the run command, then the executable will run with default values.

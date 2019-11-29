# k-Nearest Neighbors

## Description
Sequential and parallel implementation (with MPI) of k-NN.

## How to run it

1. If your gcc version is previous than gcc-7, then change it in the Makefile
2. The main function for testing is located into the tester files
    1. Run "make test_sequential" for sequential implementation
    2. Run "make test_synchronous" for MPI synchronous implementation
    3. Run "make test_asynchronous" for MPI asynchronous implementation

If you want to change the input values (size of Matrices, number of desired kNN Neighbors),
you have to modify the tester files.

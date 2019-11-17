# k-Nearest Neighbors

## Description
Sequential and parallel implementation of k-NN.

## How to run it

1. If your gcc version is previous than gcc-7, then change it in the Makefile
2. There are two ways of running it:
    1. Type "make lib" so that the Makefile produces four libraries, one for each implementation. Use those files along with inc/knnring.h at your own Makefile, with your main function, to produce executable files.
    2. Type "make" or "make exec_name" so that the Makefile produces the desired executables. That way, it uses the already existing main.c file located in the src folder.  

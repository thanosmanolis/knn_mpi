#############################################################
#    C/C++ Makefile											#
#    Author: Thanos Manolis <thanos.m14@gmail.com			#
#############################################################
#															#
#   'make lib'	  	  build the libraries .a				#
#   'make clean'  	  removes .o .a and executable files    #
#															#
#############################################################

# define the C/C++ compiler to use,default here is gcc-7
CC = gcc-7

# define the MPI compiler to use,default here is mpicc
MPICC = mpicc
MPIRUN = mpirun -np 4

# define flags
CFLAGS = -O3
LDFLAGS = -lm -lopenblas

# define command to remove files
RM = rm -rf

# always build those, even if "up-to-date"
.PHONY: lib

lib:
	cd src; $(CC) -c knnring_sequential.c $(CFLAGS); cd ..
	cd src; $(MPICC) -c knnring_synchronous.c $(CFLAGS); cd ..
	cd src; $(MPICC) -c knnring_asynchronous.c $(CFLAGS); cd ..
	cd src; ar rcs ../lib/knnring_sequential.a knnring_sequential.o; cd ..
	cd src; ar rcs ../lib/knnring_synchronous.a knnring_synchronous.o; cd ..
	cd src; ar rcs ../lib/knnring_asynchronous.a knnring_asynchronous.o; cd ..

clean:
	$(RM) src/*.o lib/*.a $(EXECS)

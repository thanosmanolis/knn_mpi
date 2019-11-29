#############################################################
#    C/C++ Makefile											#
#    Author: Thanos Manolis <thanos.m14@gmail.com			#
#############################################################
#															#
#   'make'  		  build all executable files			#
#   'make exec_name'  build executable file 'test_*'		#
#   'make lib'	  	  build the libraries .a				#
#   'make clean'  	  removes .o .a and executable files    #
#															#
#############################################################

# define the C/C++ compiler to use,default here is gcc-7
CC = gcc-7

# define the MPI compiler to use,default here is mpicc
MPICC = mpicc
MPIRUN = mpirun -np 4

# all the executables
EXECS = test_sequential test_synchronous

# define flags
CFLAGS = -O3
LDFLAGS = -lm -lopenblas

# define command to remove files
RM = rm -rf

# always build those, even if "up-to-date"
.PHONY: lib $(EXECS)

all: $(EXECS)

test_sequential:
	cd src; $(CC) main.c knnring_sequential.c -o ../$@ $(CFLAGS) $(LDFLAGS); cd ..
	./test_sequential

test_synchronous:
	cd src; $(MPICC) main.c knnring_synchronous.c -o ../$@ $(CFLAGS) $(LDFLAGS); cd ..
	$(MPIRUN) ./test_synchronous

lib:
	cd src; $(CC) -c knnring_sequential.c $(CFLAGS); cd ..
	cd src; $(MPICC) -c knnring_synchronous.c $(CFLAGS); cd ..
	cd src; ar rcs ../lib/knnring_sequential.a knnring_sequential.o; cd ..
	cd src; ar rcs ../lib/knnring_synchronous.a knnring_synchronous.o; cd ..

clean:
	$(RM) src/*.o lib/*.a $(EXECS)

#############################################################
#    C/C++ Makefile											#
#    Author: Thanos Manolis <thanos.m14@gmail.com			#
#############################################################
#															#
#   'make'  		  build all executable files			#
#   'make exec_name'  build executable file 'test_*'		#
#   'make clean'  	  removes .o .a and executable files    #
#															#
#############################################################

# define the C/C++ compiler to use,default here is gcc-7
CC = gcc

# define the MPI compiler to use, default here is mpicc
MPICC = mpicc
MPIRUN = mpirun -np 4

# all the executables
EXECS = test_sequential test_synchronous test_asynchronous

# define flags
CFLAGS = -O3
LDFLAGS = -lm -lopenblas

# define command to remove files
RM = rm -rf

# always build those, even if "up-to-date"
.PHONY: $(EXECS)

all: $(EXECS)

test_sequential:
	cd knnring; make lib; cd ..
	cd knnring; cp lib/*.a inc/*.h ../; cd ..
	$(CC) tester.c knnring_sequential.a -o $@ $(CFLAGS) $(LDFLAGS)

test_synchronous:
	cd knnring; make lib; cd ..
	cd knnring; cp lib/*.a inc/*.h ../; cd ..
	$(MPICC) tester_mpi.c knnring_synchronous.a -o $@ $(CFLAGS) $(LDFLAGS)

test_asynchronous:
	cd knnring; make lib; cd ..
	cd knnring; cp lib/*.a inc/*.h ../; cd ..
	$(MPICC) tester_mpi.c knnring_asynchronous.a -o $@ $(CFLAGS) $(LDFLAGS)

clean:
	$(RM) *.h *.a $(EXECS)

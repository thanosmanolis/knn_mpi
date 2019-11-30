/*!
  \file   tester.c
  \brief  Validate kNN ring implementation (MPI).

  \author Dimitris Floros
  \date   2019-11-25
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "mpi.h"

#include "knnring.h"
#include "tester_helper.h"

struct timeval startwtime, endwtime, start_excltime, end_excltime;
static clock_t st_time;
static clock_t en_time;
static struct tms st_cpu;
static struct tms en_cpu;
double p_time, excl_time;

/*
**********************************
*    RANDOM ALLOCATION HELPER    *
**********************************
*/

double * ralloc( int sz )
{
    double *X = (double *) malloc( sz *sizeof(double) );
    for (int i=0;i<sz;i++)
        X[i] = ( (double) (rand()) ) / (double) RAND_MAX;

    return X;
}

/*
**********************************
*    MPI TESTER MAIN FUNCTION    *
**********************************
*/

int testMPI( int const n, int const d, int const k, int const ap )
{
    int p, id;                    // MPI # processess and PID
    MPI_Status Stat;              // MPI status
    int dst, rcv, tag;            // MPI destination, receive, tag

    int isValid = 0;              // return value

    MPI_Comm_rank(MPI_COMM_WORLD, &id); // Task ID
    MPI_Comm_size(MPI_COMM_WORLD, &p);  // # tasks

    // allocate corpus for each process
    double * const corpus = (double * ) malloc( n*d * sizeof(double) );

    if (id == 0)
    {//! MASTER

        //! ======================= START POINT =======================
        gettimeofday (&start_excltime, NULL);

        //! Initialize data to begin with
        double const * const corpusAll = ralloc( n*d*p );

        //! ======================= END POINT =======================
        gettimeofday (&end_excltime, NULL);
        excl_time += (double)((end_excltime.tv_usec - start_excltime.tv_usec)/1.0e6
                  + end_excltime.tv_sec - start_excltime.tv_sec);

        //! Break to subprocesses
        for (int ip = 0; ip < p; ip++)
        {
            //! ======================= START POINT =======================
            gettimeofday (&start_excltime, NULL);

            for (int i=0; i<n; i++)
                for (int j=0; j<d; j++)
                    if (ap == COLMAJOR)
                        corpus_cm(i,j) = corpusAll_cm(i+ip*n,j);
                    else
                        corpus_rm(i,j) = corpusAll_rm(i+ip*n,j);

            //! ======================= END POINT =======================
            gettimeofday (&end_excltime, NULL);
            excl_time += (double)((end_excltime.tv_usec - start_excltime.tv_usec)/1.0e6
                      + end_excltime.tv_sec - start_excltime.tv_sec);


            //! Last chunk is mine
            if (ip == p-1)
                break;

            //! Which process to send? what tag?
            dst = ip+1;
            tag = 1;

            // send to correct process
            MPI_Send(corpus, n*d, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD);

        } // for (ip)

        //! Run distributed kNN
        knnresult const knnres = distrAllkNN( corpus, n, d, k);

        //! Prepare global kNN result object
        knnresult knnresall;
        knnresall.nidx  = (int *)   malloc( n*p*k*sizeof(int)    );
        knnresall.ndist = (double *)malloc( n*p*k*sizeof(double) );
        knnresall.m = n*p;
        knnresall.k = k;

        //! Put my results to correct spot
        for (int j = 0; j < k; j++)
            for (int i = 0; i < n; i++)
            {
                if (ap == COLMAJOR)
                {
                    knnresallnidx_cm(i+(p-1)*n,j)  = knnresnidx_cm(i,j);
                    knnresallndist_cm(i+(p-1)*n,j) = knnresndist_cm(i,j);
                }else
                {
                    knnresallnidx_rm(i+(p-1)*n,j)  = knnresnidx_rm(i,j);
                    knnresallndist_rm(i+(p-1)*n,j) = knnresndist_rm(i,j);
                }
            }

        //! Gather results back
        for (int ip = 0; ip < p-1; ip++)
        {
            rcv = ip+1;
            tag = 1;

            MPI_Recv( knnres.nidx, n*k, MPI_INT, rcv, tag, MPI_COMM_WORLD, &Stat);
            MPI_Recv( knnres.ndist, n*k, MPI_DOUBLE, rcv, tag, MPI_COMM_WORLD, &Stat);

            for (int j = 0; j < k; j++)
                for (int i = 0; i < n; i++){
                    if (ap == COLMAJOR)
                    {
                        knnresallnidx_cm(i+ip*n,j)  = knnresnidx_cm(i,j);
                        knnresallndist_cm(i+ip*n,j) = knnresndist_cm(i,j);
                    }else
                    {
                        knnresallnidx_rm(i+ip*n,j)  = knnresnidx_rm(i,j);
                        knnresallndist_rm(i+ip*n,j) = knnresndist_rm(i,j);
                    }
                }
        }

        //! ======================= START POINT =======================
        gettimeofday (&start_excltime, NULL);

        // ---------- Validate results
        isValid = validateResult( knnresall, corpusAll, corpusAll, n*p, n*p, d, k, ap );

        //! ======================= END POINT =======================
        gettimeofday (&end_excltime, NULL);
        excl_time += (double)((end_excltime.tv_usec - start_excltime.tv_usec)/1.0e6
                  + end_excltime.tv_sec - start_excltime.tv_sec);

    }else
    {//! SLAVE
        //! Get data from MASTER
        rcv = 0;
        tag = 1;

        MPI_Recv(corpus, n*d, MPI_DOUBLE, rcv, tag, MPI_COMM_WORLD, &Stat);

        //! Run distributed kNN
        knnresult const knnres = distrAllkNN( corpus, n, d, k);

        //! Send data back to MASTER
        dst = 0;
        tag = 1;

        MPI_Send(knnres.nidx, n*k, MPI_INT, dst, tag, MPI_COMM_WORLD);
        MPI_Send(knnres.ndist, n*k, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD);
    }

    //! Deallocate memory
    free( corpus );

    //! Return wheter validations passed or not
    return isValid;
}

int main(int argc, char *argv[])
{
    p_time = 0;
    excl_time = 0;

    int p, id;  // processes, PID
    int n = 0;  // # corpus elements per process
    int d;      // # dimensions
    int k;      // # neighbors

    MPI_Init(&argc, &argv);       // initialize MPI

    MPI_Comm_rank(MPI_COMM_WORLD, &id); // Task ID
    MPI_Comm_size(MPI_COMM_WORLD, &p); // # tasks

    //! By default, 100 experiments are done. You can change the
    //! nc, dc, kc values and also the n, d, k increment values.
    for(int nc=0; nc<5; nc++)
    {
        n += 2000/p;
        d = 0;

        for(int dc=0; dc<5; dc++)
        {
            d += 20;
            k = 0;

            for(int kc=0; kc<4; kc++)
            {
                k += 20;
                excl_time = 0;

                //! ======================= START POINT =======================
                gettimeofday (&startwtime, NULL);

                // ============================== RUN EXPERIMENTS
                int isValidC = testMPI( n, d, k, COLMAJOR );
                // int isValidR = testMPI( n, d, k, ROWMAJOR );

                // ============================== ONLY MASTER OUTPUTS
                //! MASTER gets result
                if (id == 0)
                {
                    printf("Tester validation: %s\n", STR_CORRECT_WRONG[isValidC]);

                    //! ======================= END POINT =======================
                    gettimeofday (&endwtime, NULL);
                    p_time = ( (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
                              + endwtime.tv_sec - startwtime.tv_sec) ) - excl_time;

                    printf(YELLOW "===== n*p: %d, d: %d, k: %d =====\n" RESET_COLOR, n*p, d, k);
                    printf(RED "%f sec\n" RESET_COLOR, p_time);
                }
            }
        }
    }

    //! Clean up
    MPI_Finalize();

    return 0;
}

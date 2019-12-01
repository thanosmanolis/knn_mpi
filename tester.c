/*!
  \file   tester.c
  \brief  Validate kNN ring implementation.

  \author Dimitris Floros
  \date   2019-11-13
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "knnring.h"
#include "tester_helper.h"

struct timeval startwtime, endwtime;
double p_time;

int main(int argc, char *argv[])
{
    int n, d, k;

    if(argc > 1)
    {
        n = atoi(argv[1]);  // # corpus
        d = atoi(argv[2]);  // # dimensions
        k = atoi(argv[3]);  // # neighbors
    }else
    {
        n = 1000;   // # corpus elements per process
        d = 10;     // # dimensions
        k = 5;      // # neighbors
    }

    int m = n;  // query

    double  * corpus = (double * ) malloc( n*d * sizeof(double) );
    double  * query  = (double * ) malloc( m*d * sizeof(double) );

    for (int i=0;i<n*d;i++)
        corpus[i] = ( (double) (rand()) ) / (double) RAND_MAX;

    for (int i=0;i<m*d;i++)
        query[i] = corpus[i];

    //! ========= START POINT =========
    gettimeofday (&startwtime, NULL);

    knnresult knnres = kNN( corpus, query, n, m, d, k );

    //! ========= END POINT =========
    gettimeofday (&endwtime, NULL);
    p_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
  		      + endwtime.tv_sec - startwtime.tv_sec);

    int isValidC = validateResult( knnres, corpus, query, n, m, d, k, COLMAJOR );
    // int isValidR = validateResult( knnres, corpus, query, n, m, d, k, ROWMAJOR );

    printf("Tester validation: %s\n", STR_CORRECT_WRONG[isValidC]);

    printf(YELLOW "===== n: %d, d: %d, k: %d =====\n" RESET_COLOR, n, d, k);
    printf(RED "Real Time: %f\n", p_time);

    free( corpus );
    free( query );

    return 0;
}

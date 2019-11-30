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
static clock_t st_time;
static clock_t en_time;
static struct tms st_cpu;
static struct tms en_cpu;
double p_time;

int main()
{
    int n = 0;  // corpus
    int m;      // query
    int d;      // dimensions
    int k;      // # neighbors

    //! By default, 100 experiments are done. You can change the
    //! nc, dc, kc values and also the n, d, k increment values.
    for(int nc=0; nc<5; nc++)
    {
        n += 2000;
        m = n;
        d = 0;

        for(int dc=0; dc<5; dc++)
        {
            d += 20;
            k = 0;

            for(int kc=0; kc<4; kc++)
            {
                k += 20;

                double  * corpus = (double * ) malloc( n*d * sizeof(double) );
                double  * query  = (double * ) malloc( m*d * sizeof(double) );

                for (int i=0;i<n*d;i++)
                    corpus[i] = ( (double) (rand()) ) / (double) RAND_MAX;

                for (int i=0;i<m*d;i++)
                    query[i] = corpus[i];

                //! ========= START POINT =========
                st_time = times(&st_cpu);
                gettimeofday (&startwtime, NULL);

                knnresult knnres = kNN( corpus, query, n, m, d, k );

                //! ========= END POINT =========
                gettimeofday (&endwtime, NULL);
                en_time = times(&en_cpu);
                p_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
              		      + endwtime.tv_sec - startwtime.tv_sec);

                int isValidC = validateResult( knnres, corpus, query, n, m, d, k, COLMAJOR );
                // int isValidR = validateResult( knnres, corpus, query, n, m, d, k, ROWMAJOR );

                printf("Tester validation: %s\n", STR_CORRECT_WRONG[isValidC]);

                printf(YELLOW "===== n: %d, d: %d, k: %d =====\n" RESET_COLOR, n, d, k);
                printf(RED "Real Time: %.3f, ", p_time);
                printf("User Time %jd, System Time %jd\n" RESET_COLOR,
                       (long)(en_cpu.tms_utime - st_cpu.tms_utime),
                       (long)(en_cpu.tms_stime - st_cpu.tms_stime));

                free( corpus );
                free( query );
            }
        }
    }

    return 0;
}

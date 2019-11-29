#include "../inc/knnring.h"

//! Clearing the shell using escape sequences
#define clear() printf("\033[H\033[J")
#define RED "\033[0;31m"
#define GREEN "\033[0;32m"
#define GREEN_BOLD "\033[1;32m"
#define YELLOW "\033[0;33m"
#define RESET_COLOR "\033[0m"

struct timeval startwtime, endwtime;
static clock_t st_time;
static clock_t en_time;
static struct tms st_cpu;
static struct tms en_cpu;

double p_time;

int main()
{
    int n=897;  // X
    int m=762;  // Y
    int d=7;    // dimensions
    int k=13;   // neighbors

    double *X = (double*)malloc(n*d * sizeof(double));
    double *Y = (double*)malloc(m*d * sizeof(double));

    for (int i=0;i<n*d;i++)
        X[i]= (double)rand()/(double)RAND_MAX;

    for (int i=0;i<m*d;i++)
        Y[i]= (double)rand()/(double)RAND_MAX;

    //! Set this timestamp as start
    st_time = times(&st_cpu);
    gettimeofday (&startwtime, NULL);

    knnresult knnres = kNN( X, Y, n, m, d, k );

    //! Set this timestamp as end
    en_time = times(&en_cpu);
    gettimeofday (&endwtime, NULL);
    p_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
  		      + endwtime.tv_sec - startwtime.tv_sec);

    printf(RED "%f sec\n" RESET_COLOR,p_time);

    printf(YELLOW "Real Time: %jd, User Time %jd, System Time %jd\n" RESET_COLOR,
           (long)(en_time - st_time),
           (long)(en_cpu.tms_utime - st_cpu.tms_utime),
           (long)(en_cpu.tms_stime - st_cpu.tms_stime));

    free(X);
    free(Y);
    free(knnres.ndist);
    free(knnres.nidx);

    return 0;
}

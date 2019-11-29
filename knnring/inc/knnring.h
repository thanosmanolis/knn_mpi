#ifndef KNN_H
#define KNN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <cblas.h>

#include <sys/time.h>
#include <sys/times.h>

struct timeval startwtime, endwtime;
static clock_t st_time;
static clock_t en_time;
static struct tms st_cpu;
static struct tms en_cpu;
double p_time;

/*
*************************************************
*    @file   knn.h                              *
*    @author amanolis <amanolis@ece.auth.gr>    *
*    @date   Wed Nov 13 16:03:30 2019           *
*    @brief  Main utility functions             *
*************************************************
*/

//! Definition of the kNN result struct
typedef struct knnresult{
    int    * nidx;     //! Indices (0-based) of nearest neighbors [m-by-k]
    double * ndist;    //! Distance of nearest neighbors          [m-by-k]
    int      m;        //! Number of query points                 [scalar]
    int      k;        //! Number of nearest neighbors            [scalar]
}knnresult;

/*
*****************************************************************
*    Compute k nearest neighbors of each point in X [n-by-d]    *
*                                                               *
*    - param X    Corpus data points      [n-by-d]              *
*    - param Y    Query data points       [m-by-d]              *
*    - param n    Number of data points   [scalar]              *
*    - param d    Number of dimensions    [scalar]              *
*    - param k    Number of neighbors     [scalar]              *
*                                                               *
*    Return The kNN result                                      *
*****************************************************************
*/

knnresult kNN(double * X, double * Y, int n, int m, int d, int k);

/*
*******************************************************
*    Compute distributed all-kNN of points in X       *
*                                                     *
*    - param X    Data points             [n-by-d]    *
*    - param n    Number of data points   [scalar]    *
*    - param d    Number of dimensions    [scalar]    *
*    - param k    Number of neighbors     [scalar]    *
*                                                     *
*    Return The kNN result                            *
*******************************************************
*/

knnresult distrAllkNN(double * X, int n, int d, int k);

/*
********************************************************
*    Functions that swap 2 numbers (double/integer)    *
********************************************************
*/

void swap_double(double* a, double* b);
void swap_int(int* a, int* b);

/*
********************************************************************
*    This function takes last element as pivot, places the pivot   *
*    element at its correct position in sorted array, and places   *
*    all smaller (smaller than pivot) to left of pivot and all     *
*    greater elements to right of pivot.                           *
********************************************************************
*/

int partition(double* arr, int* idx, int l, int r);

/*
************************************
*    Quickselect implementation    *
************************************
*/

double quickselect(double* arr, int* idx_list, int l, int r, int k, int *idx);

#endif

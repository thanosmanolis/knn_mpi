/*
***********************************************
*    k-nearest neighbors (k-NN) Sequential    *
***********************************************
*/
#include "../inc/knnring.h"
#include <mpi.h>
#include <stddef.h>

typedef struct knnresult knnresult;

knnresult distrAllkNN(double * X, int n, int d, int k)
{
    MPI_Status stat;
    MPI_Request	send_request,recv_request;
    int p, id;

    MPI_Comm_rank(MPI_COMM_WORLD, &id); // Task ID
    MPI_Comm_size(MPI_COMM_WORLD, &p); // # tasks

    knnresult knnres;

    double *Y    = (double *)malloc(n*d * sizeof(double));
    double *X_cp = (double *)malloc(n*d * sizeof(double));
    double *dist = (double *)malloc(k * sizeof(double));
    int    *idx  =    (int *)malloc(k * sizeof(int));

    memcpy(Y, X, n*d*sizeof(double));

    knnres = kNN(X, Y, n, n, d, k);

    int mul = id + p - 1;
    if(mul >= p)
        mul -= p;

    for(int i=0; i<n; i++)
        for(int z=0; z<k; z++)
            knnres.nidx[n*z + i] = knnres.nidx[n*z + i] + mul*n;

    int rcv;
    if(id == 0)
        rcv = p-1;
    else
        rcv = id-1;

    int dst;
    if(id == p-1)
        dst = 0;
    else
        dst = id+1;

    int tag = 1;

    for(int ip=0; ip<p; ip++)
    {
        if(ip > 0)
        {
            knnresult knn_temp;
            knn_temp = kNN(X, Y, n, n, d, k);

            mul--;
            if(mul < 0)
                mul = p-1;

            for(int i=0; i<n; i++)
                for(int z=0; z<k; z++)
                    knn_temp.nidx[n*z + i] = knn_temp.nidx[n*z + i] + mul*n;

            for(int i=0; i<n; i++)
            {
                int z1 = 0, z2 = 0, z3 = 0;

                // Traverse both arrays
                while (z1<k && z2<k && z3<k)
                {
                    if (knnres.ndist[n*z1 + i] < knn_temp.ndist[n*z2 + i])
                    {
                        dist[z3] = knnres.ndist[n*(z1) + i];
                        idx[z3++] = knnres.nidx[n*(z1++) + i];
                    }else
                    {
                        dist[z3] = knn_temp.ndist[n*(z2) + i];
                        idx[z3++] = knn_temp.nidx[n*(z2++) + i];
                    }
                }

                // Store remaining elements of first array
                while (z1 < k && z3<k)
                {
                    dist[z3] = knnres.ndist[n*(z1) + i];
                    idx[z3++] = knnres.nidx[n*(z1++) + i];
                }

                // Store remaining elements of second array
                while (z2 < k && z3<k)
                {
                    dist[z3] = knn_temp.ndist[n*(z2) + i];
                    idx[z3++] = knn_temp.nidx[n*(z2++) + i];
                }

                for(int z=0; z<k; z++)
                {
                    knnres.ndist[n*z + i] = dist[z];
                    knnres.nidx[n*z + i] = idx[z];
                }
            }
        }

        if(ip < p-1)
        {
            memcpy(X_cp, X, n*d*sizeof(double));

            MPI_Isend(X_cp, n*d, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, &send_request);
            MPI_Irecv(X, n*d, MPI_DOUBLE, rcv, tag, MPI_COMM_WORLD, &recv_request);

            MPI_Wait(&send_request,&stat);
            MPI_Wait(&recv_request,&stat);
        }
    }

    //! Free alocated memory
    free(Y);
    free(X_cp);
    free(dist);
    free(idx);

    return knnres;
}

knnresult kNN(double * X, double * Y, int n, int m, int d, int k)
{
    knnresult knnres;

    //! Fill the structure's m and k values
    knnres.m = m;
    knnres.k = k;

    //! Allocate the required memory
    double *dist     = malloc(m*n * sizeof(double));
    double *dist_cp  = malloc(n   * sizeof(double));
    int    *idx_list = malloc(n   * sizeof(int));
    knnres.nidx      = malloc(m*k * sizeof(int));
    knnres.ndist     = malloc(m*k * sizeof(double));

    //! dist = X^2 - 2*X*Y + Y^2
    //! Implementation of -2*X*Y with blas library
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, m, d, -2, X, n, Y, m, 0, dist, n);

    //! powX = X^2
    double *powX = malloc(n * sizeof(double));
    for(int i=0; i<n; i++)
    {
        powX[i] = 0.0;
        for(int z=0; z<d; z++)
            powX[i] += X[n*z + i] * X[n*z + i];
    }

    //! powY = Y^2
    double *powY = malloc(m * sizeof(double));
    for(int i=0; i<m; i++)
    {
        powY[i] = 0.0;
        for(int z=0; z<d; z++)
            powY[i] += Y[m*z + i] * Y[m*z + i];
    }

    //! dist += X^2 + Y^2
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<m; j++)
        {
            dist[n*j + i] += powX[i] + powY[j];

            //! If dist < 10^(-8), do it 0
            if( dist[n*j + i] < 1e-8 )
                dist[n*j + i] = 0.0;

            dist[n*j + i] = sqrt(dist[n*j + i]);
        }
    }

    //! Calculate k-nearest neighbors
    for(int i=0; i<m; i++)
    {
        //! Store the present point's distances and indices
        for(int j=0; j<n; j++)
        {
            dist_cp[j] = dist[n*i + j];
            idx_list[j] = j;
        }

        //! Perform quickselect for the whole array once, and then only for the elements
        //! on the left of the returned value. This is way more efficient
        int index;
        knnres.ndist[m*(k-1) + i] = quickselect(dist_cp, idx_list, 0, n-1, k, &index);
        knnres.nidx[m*(k-1) + i] = idx_list[index];

        for(int z=k-2; z>=0; z--)
        {
            knnres.ndist[m*z + i] = quickselect(dist_cp, idx_list, 0, z, z+1, &index);
            knnres.nidx[m*z + i] = idx_list[index];
        }
    }

    //! Free allocated memory
    free(dist);
    free(dist_cp);
    free(idx_list);

    return knnres;
}

void swap_double(double* a, double* b)
{
    double t = *a;
    *a = *b;
    *b = t;
}

void swap_int(int* a, int* b)
{
    int t = *a;
    *a = *b;
    *b = t;
}

int partition(double* arr, int* idx_list, int l, int r)
{
    double x = *(arr + r);
    int i = l;
    for (int j = l; j <= r - 1; j++)
    {
        if (*(arr + j) <= x)
        {
            swap_double(&*(arr + i), &*(arr + j));
            swap_int(&*(idx_list + i), &*(idx_list + j));
            i++;
        }
    }
    swap_double(&*(arr + i), &*(arr + r));
    swap_int(&*(idx_list + i), &*(idx_list + r));
    return i;
}

double quickselect(double* arr, int* idx_list, int l, int r, int k, int *idx)
{
    if (k > 0 && k <= r - l + 1)
    {
        int index = partition(arr, idx_list, l, r);

        if (index - l == k - 1)
        {
            *idx = index;
            return *(arr + index);
        }

        if (index - l > k - 1)
            return quickselect(arr, idx_list, l, index - 1, k, idx);

        return quickselect(arr, idx_list, index + 1, r, k - index + l - 1, idx);
    }

    return DBL_MAX;
}

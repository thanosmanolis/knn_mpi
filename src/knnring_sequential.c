/*
***********************************************
*    k-nearest neighbors (k-NN) Sequential    *
***********************************************
*/
#include "../inc/knnring.h"

typedef struct knnresult knnresult;

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

    //! dist += X^2 + Y^2
    for(int i=0; i<m; i++)
        for(int j=0; j<n; j++)
        {
            double sum = dist[n*i + j];
            for(int z=0; z<d; z++)
                sum += X[n*z + j] * X[n*z + j] + Y[m*z + i] * Y[m*z + i];

            dist[n*i + j] = sqrt(sum);

            //! If dist = NaN, do it 0
            if(dist[n*i + j] != dist[n*i + j])
                dist[n*i + j] = 0.0;
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

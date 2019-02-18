#include <omp.h>
#include <stdio.h>

int main()
{
    int a[5], i;
    #pragma omp parallel
    {
        // Perform some computation.
        #pragma omp for
        for (i = 0; i < 5; i++)
        {
            a[i] = i * i;
            printf(" A ---- Thread %d !!!\n", omp_get_thread_num());
        }
        // Print intermediate results.
        #pragma omp master
        for (i = 0; i < 5; i++)
        {
            printf("a[%d] = %d ", i, a[i]);
            printf(" ---- Thread %d !!!\n", omp_get_thread_num());
        }
        // Wait.
        #pragma omp barrier //barrier for the threads
        // Continue with the computation.
        #pragma omp for
        for (i = 0; i < 5; i++)
        {
            a[i] += i;
            printf(" B---- Thread %d !!!\n", omp_get_thread_num());
        }
    }
}

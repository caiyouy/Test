#include <stdio.h>
#include <omp.h>

int main()
{
		int coreNum=omp_get_num_procs();
		printf("Core Num is %d \n",coreNum);
#pragma omp parallel
		{
				int k=omp_get_thread_num();
				printf("ID: %d \n",k);
		}
		return 0;
}

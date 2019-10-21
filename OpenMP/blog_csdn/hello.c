#include <stdio.h>
#include <omp.h>

int main()
{
	//#pragma omp parallel
	//printf("Hello World.\n");
	
	#pragma omp parallel for
	for(int i=0;i<8;i++)
		for(int j=0;j<8;j++)
			printf("(%d,%d)\n",i,j);


	return 0;
}//omp stop when meeting bracket


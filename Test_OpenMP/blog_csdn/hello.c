#include <stdio.h>
#include <omp.h>

int main()
{
#pragma omp parallel
		printf("Hello World.\n");
		return 0;
}//omp stop when meeting bracket


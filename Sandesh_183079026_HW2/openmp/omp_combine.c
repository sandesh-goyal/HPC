#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <time.h>

//#define NUM_THREADS 2
#define N 1000
#define L 0
#define H 1000
#define RAND_MAX 2147483647		//0x7fffffff
//--------------------------------------------------GLOBAL_VARIABLE-------------------------------------------------
float a[N][N];
float b[N][N];
float c[N][N];

void print_mat();

//--------------------------------------------------------MAIN------------------------------------------------------
int main() 
{ 
	//srand(time(NULL));
	float t1, t2;
	int i,j,k;
	int num;

	//omp_set_num_threads(NUM_THREADS);
	//-----------------------------------------------------------------INITIALIZE (SERIAL)
	for(i=0; i<N; i++)
	{
		for(j=0; j<N; j++)
		{
			a[i][j] = (H-L)*((float) rand()/RAND_MAX);
			b[j][i] = (H-L)*((float) rand()/RAND_MAX);
		}
	}
	//-----------------------------------------------------------------MULTIPLY (PARALLEL)
	#pragma omp parallel for private(j,k)
	for(i=0; i<N; i++)
	{
		//num = omp_get_num_threads();
		//printf("number of threads = %d\n", num);
		for(j=0; j<N; j++)
		{
			c[i][j] = 0;
			for(k=0; k<N; k++)
			{
				c[i][j] = c[i][j] + (a[i][k] * b[j][k]);
			}
		}
	}
	//print_mat();
	//-----------------------------------------------------------------UPPER TRIANGLE TRANSFORMATION (PARALLEL)
	for(i=0; i<N-1; i++)
	{	
		#pragma omp parallel for private(k,t1,t2)
		for(j=i+1; j<N; j++)
		{	
			t1 = c[j][i];
			t2 = c[i][i];
			for(k=i; k<N; k++)
			{
				c[j][k] = c[j][k] - (float)t1/t2*c[i][k];
			}
		}
	}	
	//print_mat();
	
	return 0; 
}
//------------------------------------------------------PRINT_MAT----------------------------------------------------
void print_mat()
{
	int i,j;
	
	for(i=0; i<N; i++)
	{
		for(j=0; j<N; j++)
		{
			printf("%8.4f  ", a[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	for(i=0; i<N; i++)
	{
		for(j=0; j<N; j++)
		{
			printf("%8.4f  ", b[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	for(i=0; i<N; i++)
	{
		for(j=0; j<N; j++)
		{
			printf("%8.4f  ", c[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	printf("\n");
}

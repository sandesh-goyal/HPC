#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <mpi.h>

#define N 5000
#define TAG_MSL 1
#define TAG_MSH 2
#define TAG_MSM 3
#define TAG_SML 11
#define TAG_SMH 12
#define TAG_SMM 13
#define ROOT 0
#define L 0
#define H 1000
#define RAND_MAX 2147483647		//0x7fffffff

void initialize_mat();
void print_mat();
//--------------------------------------------------GLOBAL_VARIABLE-------------------------------------------------
int rank;
int size;
int min;
int max;
int part;
float t1;
float t2;

float a[N][N];
float b[N][N];
float c[N][N];

MPI_Status status;
MPI_Request request;

//--------------------------------------------------------MAIN------------------------------------------------------
int main(int argc, char *argv[])
{
	int i,j,k;
	double tic;
	double toc;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	tic = MPI_Wtime();
    //--------------------------------------------------------------------MATRIX MULTIPLICATION
    if(rank == ROOT)
    {	
    	initialize_mat();
        for(i=1; i<size; i++)
        {
            part = (N/(size-1));
            min = (i-1) * part;
            
            if(((i+1) == size) && ((N%(size-1)) != 0))
            {
                max = N;
            }
            else
            {
                max = min + part;
            }

            MPI_Isend(&min, 1, MPI_INT, i, TAG_MSL, MPI_COMM_WORLD, &request);
            MPI_Isend(&max, 1, MPI_INT, i, TAG_MSH, MPI_COMM_WORLD, &request);
            MPI_Isend(&a[min][0], (max-min)*N, MPI_FLOAT, i, TAG_MSM, MPI_COMM_WORLD, &request);
        }
    }

    MPI_Bcast(&b, N*N, MPI_FLOAT, ROOT, MPI_COMM_WORLD);

    if(rank > ROOT)
    {
        MPI_Recv(&min, 1, MPI_INT, 0, TAG_MSL, MPI_COMM_WORLD, &status);
        MPI_Recv(&max, 1, MPI_INT, 0, TAG_MSH, MPI_COMM_WORLD, &status);
        MPI_Recv(&a[min][0], (max-min)*N, MPI_FLOAT, 0, TAG_MSM, MPI_COMM_WORLD, &status);
        
        for (i=min; i<max; i++)
        {
       		for (j = 0; j < N; j++)
       		{
                for (k = 0; k < N; k++)
                {
                    c[i][j] = c[i][j] + (a[i][k] * b[j][k]);
                }
            }
        }

        MPI_Isend(&min, 1, MPI_INT, 0, TAG_SML, MPI_COMM_WORLD, &request);
        MPI_Isend(&max, 1, MPI_INT, 0, TAG_SMH, MPI_COMM_WORLD, &request);
        MPI_Isend(&c[min][0], (max-min)*N, MPI_FLOAT, 0, TAG_SMM, MPI_COMM_WORLD, &request);
    }

    if (rank == ROOT)
    {
        for (i = 1; i < size; i++)
        {
            MPI_Recv(&min, 1, MPI_INT, i, TAG_SML, MPI_COMM_WORLD, &status);
            MPI_Recv(&max, 1, MPI_INT, i, TAG_SMH, MPI_COMM_WORLD, &status);
            MPI_Recv(&c[min][0], (max-min)*N, MPI_FLOAT, i, TAG_SMM, MPI_COMM_WORLD, &status);
        }
        
    	//print_mat();
    }
    //--------------------------------------------------------------------UPPER TRIANGULAR TRANSFORM
    
	for(i=0; i<N-1; i++)
	{		
		if(rank == ROOT)
		{	
			part = ((N-i-1)/(size-1));
			for(j=1; j<size; j++)
		    {
		        min = ((j-1) * part) + i + 1;
		        
		        /*
		        //printf("i: %d part: %d\n", i, part);
		        //printf("i: %d min: %d\n", i, min);
		        //printf("check: 1 %d\n", i);
		        */
		        
		        if(((j+1) == size) && (((N-i-1)%(size-1)) != 0))
		        {
		            max = N;
		        }
		        else
		        {
		            max = min + part;
		        }
		        //printf("i: %d max: %d\n", i, max);
		        MPI_Isend(&min, 1, MPI_INT, j, TAG_MSL, MPI_COMM_WORLD, &request);
		        MPI_Isend(&max, 1, MPI_INT, j, TAG_MSH, MPI_COMM_WORLD, &request);
		        //MPI_Isend(&c[min][0], (max-min)*N, MPI_FLOAT, j, TAG_MSM, MPI_COMM_WORLD, &request);
		    }
		}
		//MPI_Bcast(&c[i][i], 1, MPI_FLOAT, ROOT, MPI_COMM_WORLD);
		MPI_Bcast(&c, N*N, MPI_FLOAT, ROOT, MPI_COMM_WORLD);		
		
		if(rank > ROOT)
		{
		    MPI_Recv(&min, 1, MPI_INT, 0, TAG_MSL, MPI_COMM_WORLD, &status);
		    MPI_Recv(&max, 1, MPI_INT, 0, TAG_MSH, MPI_COMM_WORLD, &status);
		    //MPI_Recv(&c[min][0], (max-min)*N, MPI_FLOAT, 0, TAG_MSM, MPI_COMM_WORLD, &status);
		    
		    t2 = c[i][i];
		    //printf("t2: %f\n", t2);
		    for(j=min; j<max; j++)
			{	
				t1 = c[j][i];
				//printf("t1: %f\n", t1);
				for(k=i; k<N; k++)
				{
					c[j][k] = c[j][k] - (float)t1/t2*c[i][k];
				}
			}

		    MPI_Isend(&min, 1, MPI_INT, 0, TAG_SML, MPI_COMM_WORLD, &request);
		    MPI_Isend(&max, 1, MPI_INT, 0, TAG_SMH, MPI_COMM_WORLD, &request);
		    MPI_Isend(&c[min][0], (max-min)*N, MPI_FLOAT, 0, TAG_SMM, MPI_COMM_WORLD, &request);
		}
		
		if (rank == ROOT)
		{
		    for (j = 1; j < size; j++)
		    {
		        MPI_Recv(&min, 1, MPI_INT, j, TAG_SML, MPI_COMM_WORLD, &status);
		        MPI_Recv(&max, 1, MPI_INT, j, TAG_SMH, MPI_COMM_WORLD, &status);
		        MPI_Recv(&c[min][0], (max-min)*N, MPI_FLOAT, j, TAG_SMM, MPI_COMM_WORLD, &status);
		    }

		}
	}
    
    /*
    if(rank == ROOT)
    {
        print_mat();
   	}
    */
    
    toc = MPI_Wtime();
    printf("Run Time for process %d = %lfs\n", rank, (double) (toc-tic));
    MPI_Finalize();
    return 0;
}
//---------------------------------------------------INITIALIZE_MAT-------------------------------------------------
void initialize_mat()
{	
	int i,j;
	
    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            a[i][j] = (H-L)*((float) rand()/RAND_MAX);
			b[j][i] = (H-L)*((float) rand()/RAND_MAX);
        }
    }
    
    //print_mat();
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

//v1.0
//**NO MEMORY OPTIMIZATION, NO SHARED MEMORY**
#include <stdio.h>

#define N 100								// Square Matrix dimension
#define BLOCK_DIM 25 						// max 32 for 1024 threads in a block, 25 is a good divisor for N
#define L 0									// LOW range for Random number
#define H 65535								// HIGH range for Random number
#define RAND_MAX 2147483647					// 0x7fffffff

//-----------------------------------------------------------------Function Declration
void random_ints(float a[N*N], float b[N*N]);
void print_mat(float arr[N*N]);
//-----------------------------------------------------------------GPU kernel for Multiplication
__global__ void mult(float *a, float *b, float *c)
{
	int rowid = (blockDim.y * blockIdx.y) + threadIdx.y;
	int columnid = (blockDim.x * blockIdx.x) + threadIdx.x;

	float sum = 0.0;
	for(int i = 0; i < N; i++)
	{
		sum += a[(rowid*N) + i] * b[(i*N) + columnid];
	}
	c[(rowid*N) + columnid] = sum;
}
//-----------------------------------------------------------------Main function
int main()
{	
	// Device copies of a, b, c
	float *dev_a, *dev_b, *dev_c; 	
	
	// 1D array size equivalent for 2D matrix		
	int size = N * N * sizeof(float); 
							
	// Allocate device copies of a, b, c			
	cudaMalloc((void**)&dev_a, size);		
	cudaMalloc((void**)&dev_b, size);
	cudaMalloc((void**)&dev_c, size);
	
	// Host copies of a, b, c, declared as 1D array instead of 2D
	// to avoid issue with pointer/address with cudaMemcpy
	float *a = (float *)malloc(size);		
	float *b = (float *)malloc(size);		
	float *c = (float *)malloc(size);
	
	// Initialize array a and b
	random_ints(a,b);						
	
	// Copy a and b to device
	cudaMemcpy(dev_a, a, size, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_b, b, size, cudaMemcpyHostToDevice);
	
	// 2D allocation of threads in block and blocks in grid, third dimesion is 1
	dim3 gridsize (N/BLOCK_DIM, N/BLOCK_DIM, 1);
  	dim3 blocksize(BLOCK_DIM, BLOCK_DIM, 1);
	
	// Launch mult() kernel with multiple blocks and threads
	mult<<<gridsize,blocksize >>>(dev_a, dev_b, dev_c);
	
	// Copy device result back to host copy of c
	cudaMemcpy(c, dev_c, size, cudaMemcpyDeviceToHost);
	
	// Print a, b, and result c=ab
	//print_mat(a);
	//print_mat(b);
	//print_mat(c);
	
	// Free host memory for a, b and c
	free(a);
	free(b);
	free(c);
	
	// Free device memory for a, b and c
	cudaFree(dev_a);
	cudaFree(dev_b);
	cudaFree(dev_c);
	return 0;
}
//-----------------------------------------------------------------Function definition for Random initialization
void random_ints(float a[N*N], float b[N*N])
{	
	int i;
	
    for(i=0; i<N*N; i++)
	{
		a[i] = (H-L)*((float) rand()/RAND_MAX);
		b[i] = (H-L)*((float) rand()/RAND_MAX);
	}
}
//-----------------------------------------------------------------Function definition for Print
void print_mat(float arr[N*N])
{	
	int i;
	for(i=0; i<N*N; i++)
	{
		printf("%f ", arr[i]);
		
		if((i+1)%N == 0)
			printf("\n");
	}
	printf("\n");
}

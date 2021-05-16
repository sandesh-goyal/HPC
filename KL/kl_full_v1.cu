/*
v1
for calculating maximum gain, square root method is used to get number of parallel threads calculating maximum gain
better for less number of MODULES
*/
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_TH 1024
//GPU KERNEL-----------------------------------------------------------------------CALCULATE INITIAL GAIN
__global__ void calc_init_gain(int *set1, int *set2, int *set1_int_gain, int *set1_ext_gain, int *set2_int_gain, int *set2_ext_gain, int *set1_d_gain, int *set2_d_gain, int *set1_id, int *set2_id, int set1_size, int set2_size, int MODULES)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j;
	
	if(i < set1_size)
	{
		set1_ext_gain[i] = 0;
    	set1_int_gain[i] = 0;
    	set2_ext_gain[i] = 0;
    	set2_int_gain[i] = 0;
		for(j=0; j<set1_size; j++)
		{
			set1_int_gain[i] += set1[(j*MODULES) + set1_id[i]];
			if (i< set2_size)
			{
				set2_ext_gain[i] += set1[(j*MODULES) + set2_id[i]];
				if (j < set2_size)
					set2_int_gain[i] += set2[(j*MODULES) + set2_id[i]];
			}
			if (j < set2_size)
			{
				set1_ext_gain[i] += set2[(j*MODULES) + set1_id[i]];
				//*initial_cutset_size += set2[(j*MODULES) + set1_id[i]];
			}
		}
		set1_d_gain[i] = set1_ext_gain[i] - set1_int_gain[i];
		if (i < set2_size)
			set2_d_gain[i] = set2_ext_gain[i] - set2_int_gain[i];
    }
}
//GPU KERNEL-----------------------------------------------------------------------CALCULATE GAIN BENEFIT
__global__ void cal_gain_benefit(int *set1, int *set2_id, int *set1_d_gain, int *set2_d_gain, int *gain_benefit, int set1_size, int set2_size, int k, int MODULES)
{
	int i = (blockDim.y * blockIdx.y) + threadIdx.y;
	int j = (blockDim.x * blockIdx.x) + threadIdx.x;
	
	if((i<set1_size-k) && (j<set2_size-k))
	{
		gain_benefit[((i*(set2_size-k))+j)*3 + 0] = i;
		gain_benefit[((i*(set2_size-k))+j)*3 + 1] = j;
		gain_benefit[((i*(set2_size-k))+j)*3 + 2] = set1_d_gain[i] + set2_d_gain[j] - (2*set1[(i*MODULES) + set2_id[j]]);
	}
}
//GPU KERNEL-----------------------------------------------------------------------CALCULATE MAX GAIN TEMP
__global__ void cal_max_g_t(int *gain_benefit, int *max_g, int *max_g_location, int t, int temp, int thread_num)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if(tid < thread_num)
	{
		int i, loc;
		loc = (tid+1)*temp;
		int my_max = gain_benefit[(tid*temp*3) + 2];
		int my_max_location = tid*temp;
		
		for(i=tid*temp; i<loc; i++)
		{
			if(i < t)
			{
				if(gain_benefit[i*3 + 2] > my_max)
				{
					my_max = gain_benefit[i*3 + 2];
					my_max_location = i;
				}
			}
		}
	
		max_g[tid] = my_max;
		max_g_location[tid] = my_max_location;
	}
}
//GPU KERNEL-----------------------------------------------------------------------CALCULATE MAX GAIN
__global__ void cal_max_g(int *max_g, int *max_g_location, int *maxt, int *max_location, int thread_num)
{
	int i;
	
	*maxt = max_g[0];
	*max_location = max_g_location[0];

	for(i=1; i<thread_num; i++)
	{
		if(max_g[i] > *maxt)
		{
			*maxt = max_g[i];
			*max_location = max_g_location[i];
		}
	}
}
//GPU KERNEL-----------------------------------------------------------------------UPDATE ITERATION
__global__ void update_iteration(int *iteration, int *set1_id, int *set2_id, int *gain_benefit, int *max_location, int *maxt, int k)
{
	iteration[(k*3) + 0] = set1_id[gain_benefit[(*max_location)*3 + 0]];
    iteration[(k*3) + 1] = set2_id[gain_benefit[(*max_location)*3 + 1]];
    iteration[(k*3) + 2] = *maxt;
    /*
    printf("----------------------Iteration %d\n", k+1);
	printf("max gain benefit: %d\n", iteration[(k*3) + 2]);
    printf("set1 swapped id: %d\n", iteration[(k*3) + 0]);
    printf("set2 swapped id: %d\n", iteration[(k*3) + 1]);
    */
}
//GPU KERNEL-----------------------------------------------------------------------GAIN SWAP
__global__ void gain_swap(int *set1, int *set1_id, int *set1_d_gain, int *set2, int *set2_id, int *set2_d_gain, int *gain_benefit, int *max_location, int MODULES, int set1_size, int set2_size, int k)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int temp;
	if(i < MODULES)
	{
		temp = set1[((gain_benefit[(*max_location)*3 + 0])*MODULES) + i];
		set1[((gain_benefit[(*max_location)*3 + 0])*MODULES) + i] = set1[((set1_size-1-k)*MODULES) + i];
		set1[((set1_size-1-k)*MODULES) + i] = temp;
	
		temp = set2[((gain_benefit[(*max_location)*3 + 1])*MODULES) + i];
		set2[((gain_benefit[(*max_location)*3 + 1])*MODULES) + i] = set2[((set2_size-1-k)*MODULES) + i];
		set2[((set2_size-1-k)*MODULES) + i] = temp;
	
		if(i == 0)
		{
				//------------------------------------------SWAP SET1 PARAMETERS
			temp = set1_id[gain_benefit[(*max_location)*3 + 0]];
			set1_id[gain_benefit[(*max_location)*3 + 0]] = set1_id[set1_size-1-k];
			set1_id[set1_size-1-k] = temp;
	
			temp = set1_d_gain[gain_benefit[(*max_location)*3 + 0]];
			set1_d_gain[gain_benefit[(*max_location)*3 + 0]] = set1_d_gain[set1_size-1-k];
			set1_d_gain[set1_size-1-k] = temp;
				//------------------------------------------SWAP SET2 PARAMETERS
			temp = set2_id[gain_benefit[(*max_location)*3 + 1]];
			set2_id[gain_benefit[(*max_location)*3 + 1]] = set2_id[set2_size-1-k];
			set2_id[set2_size-1-k] = temp;
	
			temp = set2_d_gain[gain_benefit[(*max_location)*3 + 1]];
			set2_d_gain[gain_benefit[(*max_location)*3 + 1]] = set2_d_gain[set2_size-1-k];
			set2_d_gain[set2_size-1-k] = temp;
		}
	}
}
//GPU KERNEL-----------------------------------------------------------------------GAIN UPDATE
__global__ void gain_update(int *set1, int *set2, int *set1_id, int *set2_id, int *set1_d_gain, int *set2_d_gain, int set1_size, int set2_size, int MODULES, int k)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if(tid < set1_size-k-1)
	{
		set1_d_gain[tid] = set1_d_gain[tid] + 2*(set1[((set1_size-1-k)*MODULES) + set1_id[tid]]) - 2*(set2[((set2_size-1-k)*MODULES) + set1_id[tid]]);
		if (tid < set2_size-k-1)
			set2_d_gain[tid] = set2_d_gain[tid] - 2*(set1[((set1_size-1-k)*MODULES) + set2_id[tid]]) + 2*(set2[((set2_size-1-k)*MODULES) + set2_id[tid]]);
	}
}
//GPU KERNEL-----------------------------------------------------------------------CALCULATE SET UPDATE LOCATION
__global__ void calc_set_update_location(int *iteration, int *set1_id, int *set2_id, int *set1t, int *set2t, int set1_size, int set2_size, int i)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x;
	
	if(j < set1_size)
	{
		if(iteration[(i*3) + 0] == set1_id[j])
		{
			*set1t = j;
		}
	}
	
	if(j < set2_size)
	{
		if(iteration[(i*3) + 1] == set2_id[j])
		{
			*set2t = j;
		}
	}
}
//GPU KERNEL-----------------------------------------------------------------------SET UPDATE
__global__ void set_update(int *set1, int *set2, int *set1_id, int *set2_id, int *set1t, int *set2t, int MODULES)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x;
	int temp;
		
	if(j < MODULES)
	{
		temp = set1[((*set1t)*MODULES) + j];
		set1[((*set1t)*MODULES) + j] = set2[((*set2t)*MODULES) + j];
		set2[((*set2t)*MODULES) + j] = temp;
	}
	
	if(j == 0)
	{
		temp = set1_id[*set1t];
		set1_id[*set1t] = set2_id[*set2t];
		set2_id[*set2t] = temp;
	}
}
//-------------------------------------------------------------------------------MAIN
int main(int argc, char *argv[]) 
{
	int *set1_id;
	int *set2_id;
	int *set1_int_gain;
	int *set1_ext_gain;
	int *set2_int_gain;
	int *set2_ext_gain;
	int *set1_d_gain;
	int *set2_d_gain;
	int *set1;
	int *set2;
	int *gain_benefit;
	int *iteration;

	int *dev_set1_id;
	int *dev_set2_id;
	int *dev_set1_int_gain;
	int *dev_set1_ext_gain;
	int *dev_set2_int_gain;
	int *dev_set2_ext_gain;
	int *dev_set1_d_gain;
	int *dev_set2_d_gain;
	int *dev_set1;
	int *dev_set2;
	int *dev_gain_benefit;
	int *dev_max_g;
	int *dev_max_g_location;
	int *dev_maxt;
	int *dev_max_location;
	int *dev_iteration;
	int *dev_set1t;
	int *dev_set2t;
	
	FILE *fptr;
	int i,j,k,t, temp, maxt, max_location, max_new;
	int initial_cutset_size = 0;
	int final_cutset_size = 0;
	int len = 0;
	int value = 0;
	int location = 0;
	int PINS = 0;
	int NETS = 0;
	int MODULES = 0;
	int PADS = 0;
	int CELLS = 0;
	int check_mat[1000];
	int g_max[1000];
	int g_count = 0;
	int check_len = 0;
	int set1_size = 0;
	int set2_size = 0;
	char c;
	char s[20];
	char IN_FILE[100];
    //----------------------------------------------------------Command line arguement check
    if(argc < 2)
    {
    	printf("INPUT FORMAT: executable filename.net\n");
    	exit(0);
    }
	strcpy(IN_FILE, argv[1]);
    //----------------------------------------------------------Open IN_FILE 
    fptr = fopen(IN_FILE, "r"); 
    if (fptr == NULL) 
    { 
        printf("Cannot open file \n"); 
        exit(0); 
    } 
    //----------------------------------------------------------Read #PINS #NETS #MODULES #PADS
    fgets(s, sizeof(s), fptr);			//read first line
    c = getc(fptr); 					//read next character
    
    while(c != '\n')					//read number of PINS
    {
    	PINS = PINS*10 + ((int)c - 48);
    	c = getc(fptr);
    }
    c = getc(fptr);						//read next character
    while(c != '\n')					//read number of NETS
    {
    	NETS = NETS*10 + ((int)c - 48);
    	c = getc(fptr);
    }
    c = getc(fptr); 					//read next character
    while(c != '\n')					//read number of MODULES
    {
    	MODULES = MODULES*10 + ((int)c - 48);
    	c = getc(fptr);
    }
    c = getc(fptr); 					//read next character
    while(c != '\n')    				//read number of CELLS
    {
    	CELLS = CELLS*10 + ((int)c - 48);
    	c = getc(fptr);
    }
    PADS = MODULES - CELLS;				//calculate number of PADS
    
    printf("**************DATASET DETAILS**************\n");
    //printf("PINS: \t\t%d\n", PINS);
    //printf("NETS: \t\t%d\n", NETS);
    printf("MODULES: \t%d\n", MODULES);
    printf("PADS: \t\t%d\n", PADS);
    printf("CELLS: \t\t%d\n", CELLS);
    
    //---------------------------------------------------------- HOST MEMORY ALLOCATION
    //calculate the size of set1 and set2
    set2_size = MODULES/2;				
    set1_size = MODULES - set2_size;	//set1_size >= set2_size
    
    //allocate memory for det id
    set1_id = (int *)calloc(set1_size, sizeof(int)); 
    set2_id = (int *)calloc(set2_size, sizeof(int)); 
    
    //allocate memory for internal and external gain for both sets
    set1_int_gain = (int *)calloc(set1_size, sizeof(int)); 
    set1_ext_gain = (int *)calloc(set1_size, sizeof(int)); 
    set2_int_gain = (int *)calloc(set2_size, sizeof(int)); 
    set2_ext_gain = (int *)calloc(set2_size, sizeof(int)); 
    set1_d_gain = (int *)calloc(set1_size, sizeof(int)); 
    set2_d_gain = (int *)calloc(set2_size, sizeof(int)); 
    
    //allocate memory for set1
    set1 = (int *)calloc(set1_size*MODULES, sizeof(int)); 
    for (i=0; i<set1_size; i++)
    {
         set1_id[i] = i;
    }
    //allocate memory for set2
    set2 = (int *)calloc(set2_size*MODULES, sizeof(int)); 
    for (i=0; i<set2_size; i++) 
    {
         set2_id[i] = i + set1_size;
    }    
    //allocate memory for gain benefit, 3 columns
    //SET1_ID(array location ref) SET2_ID(array location ref) GAIN_BENEFIT
    t = set1_size*set2_size;
    gain_benefit = (int *)calloc(t*3, sizeof(int)); 
    
    //allocate memory to store result of all iterations
    //SET1_ID(0)  SET2_ID(1)  SWAP_BENEFIT(2)
    iteration = (int *)calloc(set2_size*3, sizeof(int));
         
    printf("********HOST MEMORY ALLOCATION COMPLETED********\n");
    //START-----------------------------------------------------PARSER
    while (fgets(s, sizeof(s), fptr))
    {
		len = strlen(s);
		
		if(argv[1][0] == 'i')			//IBM DATASET HAS AN EXTRA SPACE
		{
			len -= 1;
		}
		
		if(s[len-2] == '1')
		{
			if(check_len > 1)			//INNER CONNECTION FOR OLD NET
			{
				for(i=0; i<check_len; i++)
				{
					for(j=0; j<check_len; j++)
					{
						if(i != j)
						{
							if(check_mat[i] > (set1_size - 1))
								set2[((check_mat[i] - set1_size)*MODULES) + check_mat[j]] += 1;
							else
								set1[(check_mat[i])*MODULES + check_mat[j]] += 1;
						}
					}
				}
			}
			
			//NEW NET CONNECTION
			value = 0;
			check_len = 0;
			for(i=1; i<len-5; i++)
			{
				value = value*10 + ((int)s[i] - 48);
			}
		
			if(s[0] == 'a')
			{
				value += PADS;
			}
			else
			{
				value -= 1;
			}
			
			continue;
		}
		
		location = 0;
		for(i=1; i<len-3; i++)
		{
			location = location*10 + ((int)s[i] - 48);
		}
		if(s[0] == 'a')
		{
			location += PADS;
		}
		else
		{
			location -= 1;
		}	

		if(value > (set1_size - 1))
			set2[((value-set1_size)*MODULES) + location] += 1;
		else	
			set1[(value)*MODULES + location] += 1;
		if(location > (set1_size - 1))
			set2[((location-set1_size)*MODULES) + value] += 1;
		else	
			set1[(location)*MODULES + value] += 1;
			
		check_mat[check_len] = location;
		check_len++;
	}
	fclose(fptr);
    //END-------------------------------------------------------PARSER
    //START-----------------------------------------------------INITIAL GAIN CALCULATE
    cudaMalloc((void**)&dev_set1_id, set1_size*sizeof(int));
    cudaMalloc((void**)&dev_set2_id, set2_size*sizeof(int));
    cudaMalloc((void**)&dev_set1_int_gain, set1_size*sizeof(int));
    cudaMalloc((void**)&dev_set1_ext_gain, set1_size*sizeof(int));
    cudaMalloc((void**)&dev_set2_int_gain, set2_size*sizeof(int));
    cudaMalloc((void**)&dev_set2_ext_gain, set2_size*sizeof(int));
    cudaMalloc((void**)&dev_set1_d_gain, set1_size*sizeof(int));
    cudaMalloc((void**)&dev_set2_d_gain, set2_size*sizeof(int));
    cudaMalloc((void**)&dev_set1, set1_size*MODULES*sizeof(int));
    cudaMalloc((void**)&dev_set2, set2_size*MODULES*sizeof(int));
    cudaMalloc((void**)&dev_gain_benefit, set1_size*set2_size*3*sizeof(int));
    
    temp = (int) sqrt(set1_size*set2_size);
    cudaMalloc((void**)&dev_max_g, temp*sizeof(int));
    cudaMalloc((void**)&dev_max_g_location, temp*sizeof(int));
    cudaMalloc((void**)&dev_maxt, sizeof(int));
    cudaMalloc((void**)&dev_max_location, sizeof(int));
    cudaMalloc((void**)&dev_iteration, set2_size*3*sizeof(int));
    cudaMalloc((void**)&dev_set1t, sizeof(int));
    cudaMalloc((void**)&dev_set2t, sizeof(int));
    printf("*******DEVICE MEMORY ALLOCATION COMPLETED*******\n");
    cudaMemcpy(dev_set1_id, set1_id, set1_size*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_set2_id, set2_id, set2_size*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_set1_int_gain, set1_int_gain, set1_size*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_set1_ext_gain, set1_ext_gain, set1_size*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_set2_int_gain, set2_int_gain, set2_size*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_set2_ext_gain, set2_ext_gain, set2_size*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_set1_d_gain, set1_d_gain, set1_size*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_set2_d_gain, set2_d_gain, set2_size*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_set1, set1, set1_size*MODULES*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_set2, set2, set2_size*MODULES*sizeof(int), cudaMemcpyHostToDevice);
    
    //START-----------------------------------------------------PRINT SET ID
    /*
    printf("-SET1--SET2-\n");
    for(i=0; i<set1_size; i++)
    {
    	if(i < set2_size)
    	{
    		printf("%d\t%d\n", set1_id[i], set2_id[i]);
    	}
    	else
    	{
    		printf("%d\n", set1_id[i]);
    	}
    }
    */
    //END-------------------------------------------------------PRINT SET ID
    //START-----------------------------------------------------PRINT SET
    /*
    printf("----------SET1----------\n");
    for(i=0; i<set1_size; i++)
    {
    	for(j=0; j<MODULES; j++)
    	{
    		printf("%d  ", set1[(i*MODULES)+j]);
    	}
    	printf("\n");
    }
    printf("----------SET2----------\n");
    for(i=0; i<set2_size; i++)
    {
    	for(j=0; j<MODULES; j++)
    	{
    		printf("%d  ", set2[(i*MODULES)+j]);
    	}
    	printf("\n");
    } 
    */
    //END-------------------------------------------------------PRINT SET
    int block_num, thread_num;
    //START-----------------------------------------------------INITIAL GAIN CALCULATE   
    block_num = ceil((double)set1_size/MAX_TH);

    calc_init_gain<<<block_num,MAX_TH>>>(dev_set1, dev_set2, dev_set1_int_gain, dev_set1_ext_gain, dev_set2_int_gain, dev_set2_ext_gain, dev_set1_d_gain, dev_set2_d_gain, dev_set1_id, dev_set2_id, set1_size, set2_size, MODULES);
    
    cudaMemcpy(set2_ext_gain, dev_set2_ext_gain, set2_size*sizeof(int), cudaMemcpyDeviceToHost);

	for(i=0; i<set2_size; i++)
    {
    	initial_cutset_size += set2_ext_gain[i];
    }
    //END-------------------------------------------------------INITIAL GAIN CALCULATE   
   	while(1)
    {
    	printf("------------------------------------RUN %d\n", g_count+1);
		for(k=0; k<set2_size; k++)
		{
			//START-----------------------------------------------------CALCULATE GAIN BENEFIT		
			dim3 blocksize(32, 32, 1);
			dim3 gridsize (ceil((double)(set2_size-k)/32), ceil((double)(set1_size-k)/32), 1);
	  		  		
			cal_gain_benefit<<<gridsize,blocksize>>>(dev_set1, dev_set2_id, dev_set1_d_gain, dev_set2_d_gain, dev_gain_benefit, set1_size, set2_size, k, MODULES);
			//END-------------------------------------------------------CALCULATE GAIN BENEFIT
			//START-----------------------------------------------------CALCULATE MAXIMUM GAIN SWAP
			t = (set1_size-k)*(set2_size-k);
			thread_num = (int) sqrt(t);
			block_num = ceil((double)thread_num/MAX_TH);
			temp = ceil((double)t/thread_num);
			
			//printf("thread_num: %d\n", thread_num);
			//printf("block_num: %d\n", block_num);
			//printf("temp: %d\n", temp);
		
			cal_max_g_t<<<block_num,MAX_TH>>>(dev_gain_benefit, dev_max_g, dev_max_g_location, t, temp, thread_num);
			cal_max_g <<<1,1>>>(dev_max_g, dev_max_g_location, dev_maxt, dev_max_location, thread_num);
			update_iteration <<<1,1>>>(dev_iteration, dev_set1_id, dev_set2_id, dev_gain_benefit, dev_max_location, dev_maxt, k);	
			//END-------------------------------------------------------CALCULATE MAXIMUM GAIN SWAP	
			//START-----------------------------------------------------SWAP
			block_num = ceil((double)MODULES/MAX_TH);
			gain_swap <<<block_num,MAX_TH>>>(dev_set1, dev_set1_id, dev_set1_d_gain, dev_set2, dev_set2_id, dev_set2_d_gain, dev_gain_benefit, dev_max_location, MODULES, set1_size, set2_size, k);	
			//END-------------------------------------------------------SWAP
			//START-----------------------------------------------------UPDATE GAIN (ONLY D)		
			if(k+1 < set2_size)
			{
				block_num = ceil((double)(set1_size-k-1)/MAX_TH);
				gain_update<<<block_num,MAX_TH>>>(dev_set1, dev_set2, dev_set1_id, dev_set2_id, dev_set1_d_gain, dev_set2_d_gain, set1_size, set2_size, MODULES, k);
			}
			//END-------------------------------------------------------UPDATE GAIN (ONLY D)
		}	
		//----------------------------------------------------------CALCULATE MAX CUMULATIVE GAIN
		cudaMemcpy(iteration, dev_iteration, set2_size*3*sizeof(int), cudaMemcpyDeviceToHost);
		maxt = 0;
		max_new = 0;
		max_location = 0;
		for(i=0; i<set2_size; i++)
		{
			max_new += iteration[(i*3) + 2];
			if(max_new > maxt)
			{
				maxt = max_new;
				max_location = i;
			}
		}
		//----------------------------------------------------------SWAP TO GET UPDATED SETS
		if(maxt > 0)
		{
			g_max[g_count] = maxt;
			g_count++;
			for(i=0; i<=max_location; i++)
			{
				block_num = ceil((double)set1_size/MAX_TH);
				calc_set_update_location<<<block_num,MAX_TH>>>(dev_iteration, dev_set1_id, dev_set2_id, dev_set1t, dev_set2t, set1_size, set2_size, i);
								
				block_num = ceil((double)MODULES/MAX_TH);
				set_update<<<block_num,MAX_TH>>>(dev_set1, dev_set2, dev_set1_id, dev_set2_id, dev_set1t, dev_set2t, MODULES);
			}
			//----------------------------------------------------------INITIAL GAIN FOR NEXT RUN
			block_num = ceil((double)set1_size/MAX_TH);
			calc_init_gain<<<block_num,MAX_TH>>>(dev_set1, dev_set2, dev_set1_int_gain, dev_set1_ext_gain, dev_set2_int_gain, dev_set2_ext_gain, dev_set1_d_gain, dev_set2_d_gain, dev_set1_id, dev_set2_id, set1_size, set2_size, MODULES);
				
		}
		else
		{
			/*
			cudaMemcpy(set1_id, dev_set1_id, set1_size*sizeof(int), cudaMemcpyDeviceToHost);
			cudaMemcpy(set2_id, dev_set2_id, set2_size*sizeof(int), cudaMemcpyDeviceToHost);
			cudaMemcpy(set1, dev_set1, set1_size*MODULES*sizeof(int), cudaMemcpyDeviceToHost);
			cudaMemcpy(set2, dev_set2, set2_size*MODULES*sizeof(int), cudaMemcpyDeviceToHost);
			*/
			break;
		}
	
	}
	//----------------------------------------------------------CALCULATE TOTAL GAIN OF ALL RUNS
	maxt = 0;
	for(i=0; i<g_count; i++)
	{
		maxt += g_max[i];
	}
	//----------------------------------------------------------CALCULATE FINAL CUT-SET SIZE
	final_cutset_size = initial_cutset_size - maxt;
	
	printf("****************FINAL RESULT***************\n");
	printf("Max Cumulative Gain: \t\t%d\n", maxt);
	//printf("Max Cumulative Gain Iteration: \t%d\n", max_location+1);
	printf("Initial Cutset Size: \t\t%d\n", initial_cutset_size);
    printf("Final Cutset Size: \t\t%d\n", final_cutset_size);
    printf("Number of Global Run: \t\t%d\n", g_count+1);
    printf("*******************************************\n");
    //START-----------------------------------------------------PRINT SET ID
    /*
    printf("------------\n");
    printf("-SET1--SET2-\n");
    for(i=0; i<set1_size; i++)
    {
    	if(i < set2_size)
    	{
    		printf("%d\t%d\n", set1_id[i], set2_id[i]);
    	}
    	else
    	{
    		printf("%d\n", set1_id[i]);
    	}
    }
    */
    //END-------------------------------------------------------PRINT SET ID
    //START-----------------------------------------------------PRINT SET
    /*
    printf("----------SET1----------\n");
    for(i=0; i<set1_size; i++)
    {
    	for(j=0; j<MODULES; j++)
    	{
    		printf("%d  ", set1[(i*MODULES)+j]);
    	}
    	printf("\n");
    }
    printf("----------SET2----------\n");
    for(i=0; i<set2_size; i++)
    {
    	for(j=0; j<MODULES; j++)
    	{
    		printf("%d  ", set2[(i*MODULES)+j]);
    	}
    	printf("\n");
    }   
    */
    //END-------------------------------------------------------PRINT SET
    //START-----------------------------------------------------FREE MEMORY
    free(set1);
    free(set2);
    free(set1_id);
    free(set2_id);
    free(set1_int_gain);
    free(set1_ext_gain);
    free(set2_int_gain);
    free(set2_ext_gain);
    free(set1_d_gain);
    free(set2_d_gain);
    free(gain_benefit);
    free(iteration);
    
    cudaFree(dev_set1_id);
	cudaFree(dev_set2_id);
	cudaFree(dev_set1_int_gain);
	cudaFree(dev_set1_ext_gain);
	cudaFree(dev_set2_int_gain);
	cudaFree(dev_set2_ext_gain);
	cudaFree(dev_set1_d_gain);
	cudaFree(dev_set2_d_gain);
	cudaFree(dev_set1);
	cudaFree(dev_set2);
	cudaFree(dev_gain_benefit);
	cudaFree(dev_max_g);
	cudaFree(dev_max_g_location);
	cudaFree(dev_maxt);
	cudaFree(dev_max_location);
	cudaFree(dev_iteration);
	cudaFree(dev_set1t);
	cudaFree(dev_set2t);
  
return 0; 
}

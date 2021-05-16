#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>

int *set1_id;
int *set2_id;
int *set1_int_gain;
int *set1_ext_gain;
int *set2_int_gain;
int *set2_ext_gain;
int *set1_d_gain;
int *set2_d_gain;
int *temp_swap;
int **set1;
int **set2;
int **gain_benefit;
int **iteration;

FILE *fptr;
int i,j,k, max, max_location, max_new;
int max_t[50];
int location_t[50];
int t, t1, t2, temp;
int th_num;
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

double tic, toc;

int main(int argc, char *argv[]) 
{	
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
    
    //----------------------------------------------------------MEMORY ALLOCATION
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
    set1 = (int **)calloc(set1_size, sizeof(int *)); 
    for (i=0; i<set1_size; i++)
    {
         set1[i] = (int *)calloc(MODULES, sizeof(int));
         set1_id[i] = i;
    }
    //allocate memory for set2
    set2 = (int **)calloc(set2_size, sizeof(int *)); 
    for (i=0; i<set2_size; i++) 
    {
         set2[i] = (int *)calloc(MODULES, sizeof(int)); 
         set2_id[i] = i + set1_size;
    }    
    //allocate memory for gain benefit, 3 columns
    //SET1_ID(array location ref) SET2_ID(array location ref) GAIN_BENEFIT
    t = set1_size*set2_size;
    gain_benefit = (int **)calloc(t, sizeof(int *)); 
    for (i=0; i<t; i++) 
         gain_benefit[i] = (int *)calloc(3, sizeof(int));
                 
    //allocate memory to temp_swap
    temp_swap = (int *)calloc(MODULES, sizeof(int));
    
    //allocate memory to store result of all iterations
    //SET1_ID(0)  SET2_ID(1)  SWAP_BENEFIT(2)
    iteration = (int **)calloc(set2_size, sizeof(int *));
    for (i=0; i<set2_size; i++) 
         iteration[i] = (int *)calloc(3, sizeof(int));
         
    printf("********MEMORY ALLOCATION COMPLETED********\n");   
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
				//#pragma omp parallel for private(j)
				for(i=0; i<check_len; i++)
				{
					for(j=0; j<check_len; j++)
					{
						if(i != j)
						{
							if(check_mat[i] > (set1_size - 1))
								set2[check_mat[i] - set1_size][check_mat[j]] += 1;
							else
								set1[check_mat[i]][check_mat[j]] += 1;
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
		//printf("location: %d\n", location);
		if(s[0] == 'a')
		{
			location += PADS;
		}
		else
		{
			location -= 1;
		}	

		if(value > (set1_size - 1))
			set2[value-set1_size][location] += 1;
		else	
			set1[value][location] += 1;
		if(location > (set1_size - 1))
			set2[location-set1_size][value] += 1;
		else	
			set1[location][value] += 1;
			
		check_mat[check_len] = location;
		check_len++;
	}
    //END-------------------------------------------------------PARSER
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
    		printf("%d  ", set1[i][j]);
    	}
    	printf("\n");
    }
    printf("----------SET2----------\n");
    for(i=0; i<set2_size; i++)
    {
    	for(j=0; j<MODULES; j++)
    	{
    		printf("%d  ", set2[i][j]);
    	}
    	printf("\n");
    } 
    */
    //END-------------------------------------------------------PRINT SET
    //START-----------------------------------------------------INITIAL GAIN CALCULATE
    //tic = omp_get_wtime();
    #pragma omp parallel for private(j) reduction (+:initial_cutset_size)
    for(i=0; i<set1_size; i++)
    {
    	for(j=0; j<set1_size; j++)
    	{
    		set1_int_gain[i] += set1[j][set1_id[i]];
    		if (i< set2_size)
    		{
    			set2_ext_gain[i] += set1[j][set2_id[i]];
    			if (j < set2_size)
    				set2_int_gain[i] += set2[j][set2_id[i]];
    		}
    		if (j < set2_size)
    		{
    			set1_ext_gain[i] += set2[j][set1_id[i]];
    			initial_cutset_size += set2[j][set1_id[i]];
    		}
    	}
    	set1_d_gain[i] = set1_ext_gain[i] - set1_int_gain[i];
    	if (i < set2_size)
    		set2_d_gain[i] = set2_ext_gain[i] - set2_int_gain[i];
    }
    //toc = omp_get_wtime();
    //printf("time for initial gain: %lf\n", (double)toc-tic);
    //END-------------------------------------------------------INITIAL GAIN CALCULATE
    while(1)
    {
    	printf("------------------------------------RUN %d\n", g_count+1);
		for(k=0; k<set2_size; k++)
		{
			t1 = set2_size-k;
			//START-----------------------------------------------------CALCULATE GAIN BENEFIT
			//tic = omp_get_wtime();
			#pragma omp parallel
			{
				int id;
				int i1, j1, lt, ti, my_max, my_location, n, start, end;
				th_num = omp_get_num_threads();
				n = (set1_size-k)/th_num;
				id = omp_get_thread_num();
			
				start = id*n;
				if(id != (th_num - 1))
				{
					end = start + n;
				}
				else
				{
					end = set1_size-k;
				}
				my_max = -9999;
				my_location = 0;
				
				for(i1=start; i1<end; i1++)
				{	
					lt = i1*t1;
					for(j1=0; j1<set2_size-k; j1++)
					{	
						gain_benefit[lt+j1][0] = i1;
						gain_benefit[lt+j1][1] = j1;
						ti = set1_d_gain[i1] + set2_d_gain[j1] - (2*set1[i1][set2_id[j1]]);
						gain_benefit[lt+j1][2] = ti;
						if(ti > my_max)
						{
							my_max = ti;
							my_location = lt+j1;
						}
					}
				}
				max_t[id] = my_max;
				location_t[id] = my_location;
			}
			//END-------------------------------------------------------CALCULATE GAIN BENEFIT
			//START-----------------------------------------------------CALCULATE MAXIMUM GAIN SWAP
			max = max_t[0];
			max_location = location_t[0];
			for(i=1; i<th_num; i++)
			{
				if(max_t[i] > max)
				{
					max = max_t[i];
					max_location = location_t[i];
				}
			}
			//toc = omp_get_wtime();
			//printf("time for max benefit: %lf\n", (double)toc-tic);
			//END-------------------------------------------------------CALCULATE MAXIMUM GAIN SWAP
			/*
			printf("----------------------Iteration %d\n", k+1);
			printf("max gain benefit: %d\n", max);
			printf("set1 swapped id: %d\n", set1_id[gain_benefit[max_location][0]]);
			printf("set2 swapped id: %d\n", set2_id[gain_benefit[max_location][1]]);
			*/
			iteration[k][0] = set1_id[gain_benefit[max_location][0]];
			iteration[k][1] = set2_id[gain_benefit[max_location][1]];
			iteration[k][2] = max;
			//START-----------------------------------------------------SWAP
				//------------------------------------------SWAP SET1 ALL PARAMETERS
			t = set1_size-1-k;
			temp_swap = set1[gain_benefit[max_location][0]];
			set1[gain_benefit[max_location][0]] = set1[t];
			set1[t] = temp_swap;
			
			temp = set1_id[gain_benefit[max_location][0]];
			set1_id[gain_benefit[max_location][0]] = set1_id[t];
			set1_id[t] = temp;
			
			temp = set1_d_gain[gain_benefit[max_location][0]];
			set1_d_gain[gain_benefit[max_location][0]] = set1_d_gain[t];
			set1_d_gain[t] = temp;
				//------------------------------------------SWAP SET2 ALL PARAMETERS
			t = set2_size-1-k;
			temp_swap = set2[gain_benefit[max_location][1]];
			set2[gain_benefit[max_location][1]] = set2[t];
			set2[t] = temp_swap;
			
			temp = set2_id[gain_benefit[max_location][1]];
			set2_id[gain_benefit[max_location][1]] = set2_id[t];
			set2_id[t] = temp;
			
			temp = set2_d_gain[gain_benefit[max_location][1]];
			set2_d_gain[gain_benefit[max_location][1]] = set2_d_gain[t];
			set2_d_gain[t] = temp;
			//END-------------------------------------------------------SWAP
			//START-----------------------------------------------------UPDATE GAIN (ONLY D)
			if(k+1 < set2_size)
			{
				t = set1_size-k-1;
				t1 = set2_size-k-1;
				//tic = omp_get_wtime();
				#pragma omp parallel for
				for(i=0; i<t; i++)
				{
					set1_d_gain[i] = set1_d_gain[i] + 2*(set1[t][set1_id[i]]) - 2*(set2[t1][set1_id[i]]);
					if (i < t1)
						set2_d_gain[i] = set2_d_gain[i] - 2*(set1[t][set2_id[i]]) + 2*(set2[t1][set2_id[i]]);
				}
				//toc = omp_get_wtime();
				//printf("time for gain update: %lf\n", (double)toc-tic);
			}
			//END-------------------------------------------------------UPDATE GAIN (ONLY D)
		}	
		//----------------------------------------------------------CALCULATE MAX CUMULATIVE GAIN
		max = 0;
		max_new = 0;
		max_location = 0;
		int set1t, set2t;
		for(i=0; i<set2_size; i++)
		{
			max_new += iteration[i][2];
			if(max_new > max)
			{
				max = max_new;
				max_location = i;
			}
		}
		//----------------------------------------------------------SWAP TO GET UPDATED SETS
		if(max > 0)
		{
			g_max[g_count] = max;
			g_count++;
			//#pragma omp parallel for private(j, set1t, set2t, temp, temp_swap)
			for(i=0; i<=max_location; i++)
			{
				#pragma omp parallel for
				for(j=0; j<set1_size; j++)
				{
					if(iteration[i][0] == set1_id[j])
						set1t = j;
					if(j < set2_size)
						if(iteration[i][1] == set2_id[j])
							set2t = j;
				}
		
				temp = set1_id[set1t];
				set1_id[set1t] = set2_id[set2t];
				set2_id[set2t] = temp;
		
				temp_swap = set1[set1t];
				set1[set1t] = set2[set2t];
				set2[set2t] = temp_swap;
			}
			//----------------------------------------------------------INITIAL GAIN FOR NEXT RUN
			#pragma omp parallel for private(j)
			for(i=0; i<set1_size; i++)
			{
				set1_ext_gain[i] = 0;
				set1_int_gain[i] = 0;
				set2_ext_gain[i] = 0;
				set2_int_gain[i] = 0;
				for(j=0; j<set1_size; j++)
				{
					set1_int_gain[i] += set1[j][set1_id[i]];
					if (i< set2_size)
					{
						set2_ext_gain[i] += set1[j][set2_id[i]];
						if (j < set2_size)
							set2_int_gain[i] += set2[j][set2_id[i]];
					}
					if (j < set2_size)
					{
						set1_ext_gain[i] += set2[j][set1_id[i]];
					}
				}
				set1_d_gain[i] = set1_ext_gain[i] - set1_int_gain[i];
				if (i < set2_size)
					set2_d_gain[i] = set2_ext_gain[i] - set2_int_gain[i];
			}
		}
		else
		{
			break;
		}
	}
	//----------------------------------------------------------CALCULATE TOTAL GAIN OF ALL RUNS
	max = 0;
	for(i=0; i<g_count; i++)
	{
		max += g_max[i];
	}
	//----------------------------------------------------------CALCULATE FINAL CUT-SET SIZE
	final_cutset_size = initial_cutset_size - max;
	
	printf("****************FINAL RESULT***************\n");
	printf("Max Cumulative Gain: \t\t%d\n", max);
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
    		printf("%d  ", set1[i][j]);
    	}
    	printf("\n");
    }
    printf("----------SET2----------\n");
    for(i=0; i<set2_size; i++)
    {
    	for(j=0; j<MODULES; j++)
    	{
    		printf("%d  ", set2[i][j]);
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
    //free(set1_int_gain);
    //free(set1_ext_gain);
    //free(set2_int_gain);
    //free(set2_ext_gain);
    //free(set1_d_gain);
    //free(set2_d_gain);
    free(gain_benefit);
    free(temp_swap);
    free(iteration);
    fclose(fptr);

return 0; 
}

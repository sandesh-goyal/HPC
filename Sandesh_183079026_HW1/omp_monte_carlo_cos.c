#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <time.h>

//#define RAND_MAX 0x7fffffff
#define steps 1000000
#define PI 3.14159
//------------------------------------------------FUNCTION_DECLARATION----------------------------------------------
double gen_random();
double calc_monte_carlo();
//--------------------------------------------------GLOBAL_VARIABLE-------------------------------------------------
double low_limit_x = -PI/2;
double high_limit_x = PI/2;
double low_limit_y = 0;
double high_limit_y = 1;
double r_num[steps][2];
//-----------------------------------------------------GEN_RANDOM---------------------------------------------------
double gen_random()
{
	int i;
	double lower_x = low_limit_x*100000;
	double upper_x = high_limit_x*100000;
	double lower_y = low_limit_y*100000;
	double upper_y = high_limit_y*100000;
	
	for(i=0; i<steps; i++)
	{
		r_num[i][0] = ((double)(rand() % (int)(upper_x - lower_x + 1)) + lower_x)/100000;
		r_num[i][1] = ((double)(rand() % (int)(upper_y - lower_y + 1)) + lower_y)/100000;
	}
}
//--------------------------------------------------CALC_MONTE_CARLO------------------------------------------------
double calc_monte_carlo()
{
	int i;
	double cosine_value;
	double check_pt_x;
	double check_pt_y;
	double check_count = 0;
	
	#pragma omp parallel for private(check_pt_x, check_pt_y) reduction(+:check_count)
	for(i=0; i<steps; i++)
	{
		check_pt_x = r_num[i][0];
		check_pt_y = r_num[i][1];
		
		if(check_pt_y <= cos(check_pt_x))
		{
			check_count = check_count + 1;
		}
		
	}
	
	cosine_value = (double)(high_limit_x - low_limit_x)*((double)(check_count)/(steps));
	return(cosine_value);
}
//--------------------------------------------------------MAIN------------------------------------------------------
int main() 
{ 
	double cosine_value;
	//srand(time(NULL));
	
	gen_random();
	cosine_value = calc_monte_carlo();
	printf("Integration of COS(X) from -PI/2 to PI/2 with %d steps = %lf\n", steps, cosine_value);

	return 0; 
}

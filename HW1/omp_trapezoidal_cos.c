#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <time.h>

#define steps 1000000
#define PI 3.14159
//------------------------------------------------FUNCTION_DECLARATION----------------------------------------------
double calc_trapezoidal();
//--------------------------------------------------GLOBAL_VARIABLE-------------------------------------------------
double low_limit_x = -PI/2;
double high_limit_x = PI/2;
//--------------------------------------------------CALC_TRAPEZOIDAL------------------------------------------------
double calc_trapezoidal()
{
	int i;
	double delta_x;
	double cosine_value;
	
	delta_x = (high_limit_x - low_limit_x)/steps;
	
	cosine_value = cos(low_limit_x);
	
	#pragma omp parallel for reduction (+:cosine_value)
	for(i = -steps/2; i < steps/2; i++)
	{
		cosine_value = cosine_value + 2*cos(i*delta_x);
	}
	
	cosine_value = cosine_value + cos(high_limit_x);
	cosine_value = (cosine_value*delta_x)/2;
	return(cosine_value);
}
//--------------------------------------------------------MAIN------------------------------------------------------
int main() 
{ 
	double cosine_value;

	cosine_value = calc_trapezoidal();
	printf("Integration of COS(X) from -PI/2 to PI/2 with %d steps = %lf\n", steps, cosine_value);
	
	return 0; 
}


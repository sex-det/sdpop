#include <stdio.h>
#include <float.h> 
#include <math.h>
 
int main(void)
{
    printf("LDBL_MIN      = %Le\n", 		LDBL_MIN);
    printf("LDBL_MAX      = %Le\n", 		LDBL_MAX);
    printf("LDBL_EPSILON  = %Le\n", 		LDBL_EPSILON);
    printf("LDBL_DIG      = %d\n", 		LDBL_DIG);
    printf("LDBL_MANT_DIG = %d\n", 		LDBL_MANT_DIG);
    printf("LDBL_MIN_EXP  = %d\n",  		LDBL_MIN_EXP);
    printf("LDBL_MIN_10_EXP  = %d\n",  	LDBL_MIN_10_EXP);
    printf("LDBL_MAX_EXP     = %d\n",  	LDBL_MAX_EXP);
    printf("LDBL_MAX_10_EXP  = %d\n",  	LDBL_MAX_10_EXP);
}
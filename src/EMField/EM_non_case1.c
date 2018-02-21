#include "APT_AllHeaders.h"

int GAPS_APT_Field_EM_non_case1(double *pTensor,double *pSpaceTime4,int Order,Gaps_IO_InputsContainer *pInputs)
{
	int MaxOrder = 3;
	double tt=pSpaceTime4[0],xx=pSpaceTime4[1],yy=pSpaceTime4[2],zz=pSpaceTime4[3];
	double a = 5;
	double b = 3;
	double beta = 1;
	double r = sqrt(xx*xx + yy*yy);
	double 	Ct = xx/r;
	double 	St = yy/r;
	if(-1 == Order)
	{
		if(pInputs->EMField_Cal_E)
		{	

			pTensor[0] = - ( ( -2 * beta ) / ( r*r*r * ( a + b * (2*Ct*Ct-1) ) ) * (Ct) + ( 2 * b * beta * 2*St*Ct ) / ( r*r*r * pow((a+b*(2*Ct*Ct-1)), 2) ) * (-St) );
			pTensor[1] = - ( ( -2 * beta ) / ( r*r*r * ( a + b * (2*Ct*Ct-1) ) ) * (St) + ( 2 * b * beta * 2*St*Ct ) / ( r*r*r * pow((a+b*(2*Ct*Ct-1)), 2) ) * (Ct)  );
			pTensor[2] = 0;
		}

		if(pInputs->EMField_Cal_B)
		{
			pTensor[3] = 0;
			pTensor[4] = 0;
			pTensor[5] =(b*b - a*a) / ( r*r*r*r * pow( (a + b * (2 * Ct*Ct -1)) ,1.5) );
		}
	}
	if(MaxOrder<Order)
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_Field_Uniform. This field function does NOT support tensors order larger than %d.\n",MaxOrder);
	}
	return 0;
}

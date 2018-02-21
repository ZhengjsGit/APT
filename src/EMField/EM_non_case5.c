#include "APT_AllHeaders.h"

int GAPS_APT_Field_EM_non_case5(double *pTensor,double *pSpaceTime4,int Order,Gaps_IO_InputsContainer *pInputs)
{
	int MaxOrder = 3;
	double tt=pSpaceTime4[0],xx=pSpaceTime4[1],yy=pSpaceTime4[2],zz=pSpaceTime4[3];
	double r = sqrt(xx*xx + yy*yy);
	double 	Ct = xx/r;
	double 	St = yy/r;
	if(-1 == Order)
	{

		double rsquare = xx*xx + yy*yy;
		if(pInputs->EMField_Cal_E)
		{	
			pTensor[0] = -( xx* (-2 - 2/pow( (xx*xx + yy*yy), 2) + 3*sqrt(xx*xx + yy*yy) ) );
			pTensor[1] = -( yy* (-2 - 2/pow( (xx*xx + yy*yy), 2) + 3*sqrt(xx*xx + yy*yy) ) );
			pTensor[2] = 0;
		}

		if(pInputs->EMField_Cal_B)
		{
			pTensor[3] = 0;
			pTensor[4] = 0;
			pTensor[5] = 2.5*r*r*r + xx*yy +xx*r - yy*yy*r;
		}
	}
	if(MaxOrder<Order)
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_Field_Uniform. This field function does NOT support tensors order larger than %d.\n",MaxOrder);
	}
	return 0;
}

#include "APT_AllHeaders.h"

int GAPS_APT_Field_EMhairerEquadBnc3(double *pTensor,double *pSpaceTime4,int Order,Gaps_IO_InputsContainer *pInputs)
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
			pTensor[0] = -(2*xx - 3);
			pTensor[1] = -(2*yy - 3);
			pTensor[2] = 0;
		}

		if(pInputs->EMField_Cal_B)
		{
			pTensor[3] = 0;
			pTensor[4] = 0;
			pTensor[5] = 3*(xx*yy+xx)/sqrt(rsquare);
		}
	}
	if(MaxOrder<Order)
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_Field_Uniform. This field function does NOT support tensors order larger than %d.\n",MaxOrder);
	}
	return 0;
}

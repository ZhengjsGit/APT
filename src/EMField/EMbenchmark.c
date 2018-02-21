#include "APT_AllHeaders.h"

int GAPS_APT_Field_EMbenchmark(double *pTensor,double *pSpaceTime4,int Order,Gaps_IO_InputsContainer *pInputs)
{
	int MaxOrder = 3;
	double tt=pSpaceTime4[0],xx=pSpaceTime4[1],yy=pSpaceTime4[2],zz=pSpaceTime4[3];
	if(-1 == Order)
	{
		if(pInputs->EMField_Cal_E)
		{	

			pTensor[0] = 0.01 * xx * pow( (xx * xx + yy * yy) , -1.5);
			pTensor[1] = 0.01 * yy * pow( (xx * xx + yy * yy) , -1.5);
			pTensor[2] = 0;
			//printf("E1=%f",pTensor[0]);
			//printf("E2=%f",pTensor[1]);
		}

		if(pInputs->EMField_Cal_B)
		{
			pTensor[3] = 0;
			pTensor[4] = 0;
			pTensor[5] = sqrt(xx*xx + yy*yy);
			//printf("B3=%f",pTensor[5]);
		}
	}
	if(MaxOrder<Order)
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_Field_Uniform. This field function does NOT support tensors order larger than %d.\n",MaxOrder);
	}
	return 0;
}

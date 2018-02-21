#include "APT_AllHeaders.h"

int GAPS_APT_Field_uncompact(double *pTensor,double *pSpaceTime4,int Order,Gaps_IO_InputsContainer *pInputs)
{
	int MaxOrder = 3;
	double tt=pSpaceTime4[0],xx=pSpaceTime4[1],yy=pSpaceTime4[2],zz=pSpaceTime4[3];
	if(-1 == Order)
	{
		if(pInputs->EMField_Cal_E)
		{	

			pTensor[0] =   10*xx;
			pTensor[1] =   10*yy;
			pTensor[2] = - 20*zz;
		}

		if(pInputs->EMField_Cal_B)
		{
			pTensor[3] = 100.0/3.0 + 50* (yy - zz);
			pTensor[4] = 0 + 50 * (xx + zz);
			pTensor[5] = 100 + 50 * (yy - xx);
			//printf("B3=%f",pTensor[5]);
		}
	}

	if(1 == Order)
	{
		//not true
		pTensor[0] = - yy * sqrt(xx*xx + yy*yy);
		pTensor[1] =  xx * sqrt(xx*xx + yy*yy);
		pTensor[2] = 0;
		pTensor[3] = 0;
	}
	if(MaxOrder<Order)
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_Field_Uniform. This field function does NOT support tensors order larger than %d.\n",MaxOrder);
	}
	return 0;
}

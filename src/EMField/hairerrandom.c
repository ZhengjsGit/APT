#include "APT_AllHeaders.h"

int GAPS_APT_Field_hairerrandom(double *pTensor,double *pSpaceTime4,int Order,Gaps_IO_InputsContainer *pInputs)
{
	int MaxOrder = 3;
	double tt=pSpaceTime4[0],xx=pSpaceTime4[1],yy=pSpaceTime4[2],zz=pSpaceTime4[3];
	if(-1 == Order)
	{
		if(pInputs->EMField_Cal_E)
		{	
			pTensor[0] = - (3*xx*xx + 4.0/5.0*xx*xx*xx);
			pTensor[1] = - (-3*yy*yy + 4*yy*yy*yy);
			pTensor[2] = - (4*zz*zz*zz);
		}

		if(pInputs->EMField_Cal_B)
		{
			pTensor[3] = 0; //0.5 * (yy-zz);
			pTensor[4] = 0;//0.5 * (xx+zz);
			pTensor[5] = 0.5 * (yy-xx);
			//printf("B3=%f",pTensor[5]);
		}
	}
	if(MaxOrder<Order)
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_Field_Uniform. This field function does NOT support tensors order larger than %d.\n",MaxOrder);
	}
	return 0;
}

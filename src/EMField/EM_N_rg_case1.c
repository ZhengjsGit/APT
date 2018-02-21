#include "APT_AllHeaders.h"

int GAPS_APT_Field_EM_N_rg_case1(double *pTensor,double *pSpaceTime4,int Order,Gaps_IO_InputsContainer *pInputs)
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
			pTensor[0] = r * 4*xx - 3 * rsquare * 4*xx;
			pTensor[1] = r * 4*yy - 3 * rsquare * 4*yy;
			pTensor[2] = 0;
		}

		if(pInputs->EMField_Cal_B)
		{
			pTensor[3] = 0;
			pTensor[4] = 0;
			pTensor[5] = xx/2/sqrt(rsquare);
		}
	}
	if(1 == Order)
	{
		pTensor[0] = - xx*yy /r/2;
		pTensor[1] = - yy*yy /r/2;
		pTensor[2] = 0;
		pTensor[3] = 0;
	}
	if(MaxOrder<Order)
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_Field_Uniform. This field function does NOT support tensors order larger than %d.\n",MaxOrder);
	}
	return 0;
}

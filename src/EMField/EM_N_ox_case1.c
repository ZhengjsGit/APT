#include "APT_AllHeaders.h"

int GAPS_APT_Field_EM_N_ox_case1(double *pTensor,double *pSpaceTime4,int Order,Gaps_IO_InputsContainer *pInputs)
{
	int MaxOrder = 3;
	double tt=pSpaceTime4[0],xx=pSpaceTime4[1],yy=pSpaceTime4[2],zz=pSpaceTime4[3];
	double r = sqrt(xx*xx + yy*yy);
	double 	Ct = xx/r;
	double 	St = yy/r;
	if(-1 == Order)
	{
		if(pInputs->EMField_Cal_E)
		{	

			pTensor[0] = (4*xx)/pow(4*pow(xx,2) + pow(yy,2),2);
			pTensor[1] = yy/pow(4*pow(xx,2) + pow(yy,2),2);
			pTensor[2] = 0;
		}

		if(pInputs->EMField_Cal_B)
		{
			pTensor[3] = 0;
			pTensor[4] = 0;
//			pTensor[5] =(5*Ct*Ct-1)*r*r;//(b*b - a*a) / ( r*r*r*r * pow( (a + b * (2 * Ct*Ct -1)) ,1.5) );
			pTensor[5] =4*xx*xx - yy*yy;//(b*b - a*a) / ( r*r*r*r * pow( (a + b * (2 * Ct*Ct -1)) ,1.5) );
		}
	}
	if(1 == Order)
	{
		pTensor[0] = -(1/4.0) *yy * r*r - 5 * xx * r*r * ( 1/2.0 * atan(xx/yy) - 1/4.0 * sin(2*atan(xx/yy)) );
		pTensor[1] = (1/4.0) *xx * r*r - 5 * yy * r*r * ( 1/2.0 * atan(xx/yy) - 1/4.0 * sin(2*atan(xx/yy)) );
		pTensor[2] = 0;
		pTensor[3] = 0;
	}
	if(MaxOrder<Order)
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_Field_Uniform. This field function does NOT support tensors order larger than %d.\n",MaxOrder);
	}
	return 0;
}

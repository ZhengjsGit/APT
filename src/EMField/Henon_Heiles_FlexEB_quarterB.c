#include "APT_AllHeaders.h"

int GAPS_APT_Field_Henon_Heiles_FlexEB_quarterB(double *pTensor,double *pSpaceTime4,int Order,Gaps_IO_InputsContainer *pInputs)
{
	int MaxOrder = 3;
	double tt=pSpaceTime4[0],xx=pSpaceTime4[1],yy=pSpaceTime4[2],zz=pSpaceTime4[3];
	double lamb = pInputs-> HenonLambda;
	double scaleE = pInputs->EMField_scaleE;
	double scaleB = pInputs->EMField_scaleB;
	double B0=pInputs->EMField_B0;
	double E0=pInputs->EMField_E0;
	if(-1 == Order)
	{
		if(pInputs->EMField_Cal_E)
		{
			pTensor[0] = -2.*scaleE*xx*(0.5 + lamb*yy);
			pTensor[1] = -1.*scaleE*(1.*yy + lamb*(pow(xx,2) - 1.*pow(yy,2)));
			pTensor[2] = 0;
		}

		if(pInputs->EMField_Cal_B)
		{
			pTensor[3] = 0;
			pTensor[4] = 0;
			pTensor[5] = (scaleB*(1.*pow(xx,2) + 1.*pow(yy,2)))/pow(pow(xx,2) + pow(yy,2),0.75);
		}
	}
	if(1 == Order)
	{
		pTensor[0] = scaleE*(pow(yy,2)*(0.5 - 0.3333333333333333*lamb*yy) + pow(xx,2)*(0.5 + lamb*yy));
		pTensor[1] = -0.4*scaleB*yy*pow(pow(xx,2) + pow(yy,2),0.25);
		pTensor[2] = 0.4*scaleB*xx*pow(pow(xx,2) + pow(yy,2),0.25);
		pTensor[3] = 0;
	}
	if(2 == Order)
	{
		pTensor[0] = 0;
		pTensor[1] = 0;
		pTensor[2] = 0;
		pTensor[3] = 0;
		pTensor[4] = 2*scaleE*xx*(0.5 + lamb*yy);
		pTensor[5] = (-0.2*scaleB*xx*yy)/pow(pow(xx,2) + pow(yy,2),0.75);
		pTensor[6] = (scaleB*(0.6000000000000001*pow(xx,2) + 0.4*pow(yy,2)))/pow(pow(xx,2) + pow(yy,2),0.75);
		pTensor[7] = 0;
		pTensor[8] = scaleE*(1.*yy + lamb*(pow(xx,2) - 1.*pow(yy,2)));
		pTensor[9] = (scaleB*(-0.4*pow(xx,2) - 0.6000000000000001*pow(yy,2)))/pow(pow(xx,2) + pow(yy,2),0.75);
		pTensor[10] = (0.2*scaleB*xx*yy)/pow(pow(xx,2) + pow(yy,2),0.75);
		pTensor[11] = 0;
		pTensor[12] = 0;
		pTensor[13] = 0;
		pTensor[14] = 0;
		pTensor[15] = 0;
	}
	if(3 == Order)
	{
		pTensor[0] = 0;
		pTensor[1] = 0;
		pTensor[2] = 0;
		pTensor[3] = 0;
		pTensor[4] = 0;
		pTensor[5] = 0;
		pTensor[6] = 0;
		pTensor[7] = 0;
		pTensor[8] = 0;
		pTensor[9] = 0;
		pTensor[10] = 0;
		pTensor[11] = 0;
		pTensor[12] = 0;
		pTensor[13] = 0;
		pTensor[14] = 0;
		pTensor[15] = 0;
		pTensor[16] = 0;
		pTensor[17] = 0;
		pTensor[18] = 0;
		pTensor[19] = 0;
		pTensor[20] = 2*scaleE*(0.5 + lamb*yy);
		pTensor[21] = (scaleB*(0.10000000000000005*pow(xx,2)*yy - 0.2*pow(yy,3)))/pow(pow(xx,2) + pow(yy,2),1.75);
		pTensor[22] = (scaleB*(0.30000000000000004*pow(xx,3) + 0.6000000000000001*xx*pow(yy,2)))/pow(pow(xx,2) + pow(yy,2),1.75);
		pTensor[23] = 0;
		pTensor[24] = 2*lamb*scaleE*xx;
		pTensor[25] = (scaleB*(-0.2*pow(xx,3) + 0.10000000000000005*xx*pow(yy,2)))/pow(pow(xx,2) + pow(yy,2),1.75);
		pTensor[26] = (scaleB*(-0.10000000000000005*pow(xx,2)*yy + 0.2*pow(yy,3)))/pow(pow(xx,2) + pow(yy,2),1.75);
		pTensor[27] = 0;
		pTensor[28] = 0;
		pTensor[29] = 0;
		pTensor[30] = 0;
		pTensor[31] = 0;
		pTensor[32] = 0;
		pTensor[33] = 0;
		pTensor[34] = 0;
		pTensor[35] = 0;
		pTensor[36] = 2*lamb*scaleE*xx;
		pTensor[37] = (scaleB*(-0.2*pow(xx,3) + 0.10000000000000005*xx*pow(yy,2)))/pow(pow(xx,2) + pow(yy,2),1.75);
		pTensor[38] = (scaleB*(-0.10000000000000005*pow(xx,2)*yy + 0.2*pow(yy,3)))/pow(pow(xx,2) + pow(yy,2),1.75);
		pTensor[39] = 0;
		pTensor[40] = scaleE*(1. - 2.*lamb*yy);
		pTensor[41] = (scaleB*(-0.6000000000000001*pow(xx,2)*yy - 0.30000000000000004*pow(yy,3)))/pow(pow(xx,2) + pow(yy,2),1.75);
		pTensor[42] = (scaleB*(0.2*pow(xx,3) - 0.10000000000000005*xx*pow(yy,2)))/pow(pow(xx,2) + pow(yy,2),1.75);
		pTensor[43] = 0;
		pTensor[44] = 0;
		pTensor[45] = 0;
		pTensor[46] = 0;
		pTensor[47] = 0;
		pTensor[48] = 0;
		pTensor[49] = 0;
		pTensor[50] = 0;
		pTensor[51] = 0;
		pTensor[52] = 0;
		pTensor[53] = 0;
		pTensor[54] = 0;
		pTensor[55] = 0;
		pTensor[56] = 0;
		pTensor[57] = 0;
		pTensor[58] = 0;
		pTensor[59] = 0;
		pTensor[60] = 0;
		pTensor[61] = 0;
		pTensor[62] = 0;
		pTensor[63] = 0;
	}
	if(MaxOrder<Order)
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_Field_Henon_Heiles_FlexEB_quarterB. This field function does NOT support tensors order larger than %d.\n",MaxOrder);
	}
	return 0;
}

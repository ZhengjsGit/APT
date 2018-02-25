#include "APT_AllHeaders.h"

int GAPS_APT_Field_Hairer51_FlexBE(double *pTensor,double *pSpaceTime4,int Order,Gaps_IO_InputsContainer *pInputs)
{
	int MaxOrder = 3;
	double tt=pSpaceTime4[0],xx=pSpaceTime4[1],yy=pSpaceTime4[2],zz=pSpaceTime4[3];
	double scaleE = pInputs->EMField_Hairer51_FlexBE_scaleE;
	double scaleB = pInputs->EMField_Hairer51_FlexBE_scaleB;
	double B0=pInputs->EMField_B0;
	double E0=pInputs->EMField_E0;
	if(-1 == Order)
	{
		if(pInputs->EMField_Cal_E)
		{
			pTensor[0] = scaleE*(-3 - 0.8*xx)*pow(xx,2);
			pTensor[1] = scaleE*(3 - 4*yy)*pow(yy,2);
			pTensor[2] = -4*scaleE*pow(zz,3);
		}

		if(pInputs->EMField_Cal_B)
		{
			pTensor[3] = 0;
			pTensor[4] = 0;
			pTensor[5] = scaleB*sqrt(pow(xx,2) + pow(yy,2));
		}
	}
	if(1 == Order)
	{
		pTensor[0] = scaleE*(pow(xx,3) + 0.2*pow(xx,4) - pow(yy,3) + pow(yy,4) + pow(zz,4));
		pTensor[1] = -(scaleB*yy*sqrt(pow(xx,2) + pow(yy,2)))/3.;
		pTensor[2] = (scaleB*xx*sqrt(pow(xx,2) + pow(yy,2)))/3.;
		pTensor[3] = 0;
	}
	if(2 == Order)
	{
		pTensor[0] = 0;
		pTensor[1] = 0;
		pTensor[2] = 0;
		pTensor[3] = 0;
		pTensor[4] = scaleE*(3 + 0.8*xx)*pow(xx,2);
		pTensor[5] = -(scaleB*xx*yy)/(3.*sqrt(pow(xx,2) + pow(yy,2)));
		pTensor[6] = (scaleB*(2*pow(xx,2) + pow(yy,2)))/(3.*sqrt(pow(xx,2) + pow(yy,2)));
		pTensor[7] = 0;
		pTensor[8] = scaleE*pow(yy,2)*(-3 + 4*yy);
		pTensor[9] = -(scaleB*(pow(xx,2) + 2*pow(yy,2)))/(3.*sqrt(pow(xx,2) + pow(yy,2)));
		pTensor[10] = (scaleB*xx*yy)/(3.*sqrt(pow(xx,2) + pow(yy,2)));
		pTensor[11] = 0;
		pTensor[12] = 4*scaleE*pow(zz,3);
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
		pTensor[20] = scaleE*xx*(6 + 2.4000000000000004*xx);
		pTensor[21] = -(scaleB*pow(yy,3))/(3.*pow(pow(xx,2) + pow(yy,2),1.5));
		pTensor[22] = (scaleB*(2*pow(xx,3) + 3*xx*pow(yy,2)))/(3.*pow(pow(xx,2) + pow(yy,2),1.5));
		pTensor[23] = 0;
		pTensor[24] = 0;
		pTensor[25] = -(scaleB*pow(xx,3))/(3.*pow(pow(xx,2) + pow(yy,2),1.5));
		pTensor[26] = (scaleB*pow(yy,3))/(3.*pow(pow(xx,2) + pow(yy,2),1.5));
		pTensor[27] = 0;
		pTensor[28] = 0;
		pTensor[29] = 0;
		pTensor[30] = 0;
		pTensor[31] = 0;
		pTensor[32] = 0;
		pTensor[33] = 0;
		pTensor[34] = 0;
		pTensor[35] = 0;
		pTensor[36] = 0;
		pTensor[37] = -(scaleB*pow(xx,3))/(3.*pow(pow(xx,2) + pow(yy,2),1.5));
		pTensor[38] = (scaleB*pow(yy,3))/(3.*pow(pow(xx,2) + pow(yy,2),1.5));
		pTensor[39] = 0;
		pTensor[40] = 6*scaleE*yy*(-1 + 2*yy);
		pTensor[41] = -(scaleB*(3*pow(xx,2)*yy + 2*pow(yy,3)))/(3.*pow(pow(xx,2) + pow(yy,2),1.5));
		pTensor[42] = (scaleB*pow(xx,3))/(3.*pow(pow(xx,2) + pow(yy,2),1.5));
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
		pTensor[60] = 12*scaleE*pow(zz,2);
		pTensor[61] = 0;
		pTensor[62] = 0;
		pTensor[63] = 0;
	}
	if(MaxOrder<Order)
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_Field_Hairer51_FlexBE. This field function does NOT support tensors order larger than %d.\n",MaxOrder);
	}
	return 0;
}

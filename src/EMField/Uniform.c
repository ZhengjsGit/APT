#include "APT_AllHeaders.h"

int GAPS_APT_Field_Uniform(double *pTensor,double *pSpaceTime4,int Order,Gaps_IO_InputsContainer *pInputs)
{
	int MaxOrder = 3;
	double tt=pSpaceTime4[0],xx=pSpaceTime4[1],yy=pSpaceTime4[2],zz=pSpaceTime4[3];
	double AngleEB = pInputs->EMField_Uniform_AngleEB;
	double B0=pInputs->EMField_B0;
	double E0=pInputs->EMField_E0;
	if(-1 == Order)
	{
		if(pInputs->EMField_Cal_E)
		{
			pTensor[0] = 0;
			pTensor[1] = E0*sin(AngleEB);
			pTensor[2] = E0*cos(AngleEB);
		}

		if(pInputs->EMField_Cal_B)
		{
			pTensor[3] = 0;
			pTensor[4] = 0;
			pTensor[5] = B0;
		}
	}
	if(1 == Order)
	{
		pTensor[0] = 0;
		pTensor[1] = -(B0*yy)/2.;
		pTensor[2] = (B0*xx)/2. - E0*tt*sin(AngleEB);
		pTensor[3] = -(E0*tt*cos(AngleEB));
	}
	if(2 == Order)
	{
		pTensor[0] = 0;
		pTensor[1] = 0;
		pTensor[2] = -(E0*sin(AngleEB));
		pTensor[3] = -(E0*cos(AngleEB));
		pTensor[4] = 0;
		pTensor[5] = 0;
		pTensor[6] = B0/2.;
		pTensor[7] = 0;
		pTensor[8] = 0;
		pTensor[9] = -B0/2.;
		pTensor[10] = 0;
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
		pTensor[20] = 0;
		pTensor[21] = 0;
		pTensor[22] = 0;
		pTensor[23] = 0;
		pTensor[24] = 0;
		pTensor[25] = 0;
		pTensor[26] = 0;
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
		pTensor[37] = 0;
		pTensor[38] = 0;
		pTensor[39] = 0;
		pTensor[40] = 0;
		pTensor[41] = 0;
		pTensor[42] = 0;
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
		fprintf(stderr,"ERROR: In function GAPS_APT_Field_Uniform. This field function does NOT support tensors order larger than %d.\n",MaxOrder);
	}
	return 0;
}
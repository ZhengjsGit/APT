#include "APT_AllHeaders.h"

int GAPS_APT_Field_Zjs_E_drift(double *pTensor,double *pSpaceTime4,int Order,Gaps_IO_InputsContainer *pInputs)
{
	int MaxOrder = 3;
	double tt=pSpaceTime4[0],xx=pSpaceTime4[1],yy=pSpaceTime4[2],zz=pSpaceTime4[3];
	double scaleE = pInputs->EMField_scaleE;
	double scaleB = pInputs->EMField_scaleB;
	if(-1 == Order)
	{
		if(pInputs->EMField_Cal_E)
		{
			pTensor[0] = (scaleE*xx*(-3*pow(xx,6) - 3*pow(yy,6) + 2*sqrt(pow(xx,2) + pow(yy,2)) + 2*pow(yy,4)*sqrt(pow(xx,2) + pow(yy,2)) + pow(xx,4)*(-9*pow(yy,2) + 2*sqrt(pow(xx,2) + pow(yy,2))) + pow(xx,2)*(-9*pow(yy,4) + 4*pow(yy,2)*sqrt(pow(xx,2) + pow(yy,2)))))/pow(pow(xx,2) + pow(yy,2),2.5);
			pTensor[1] = (scaleE*yy*(-3*pow(xx,6) - 3*pow(yy,6) + 2*sqrt(pow(xx,2) + pow(yy,2)) + 2*pow(yy,4)*sqrt(pow(xx,2) + pow(yy,2)) + pow(xx,4)*(-9*pow(yy,2) + 2*sqrt(pow(xx,2) + pow(yy,2))) + pow(xx,2)*(-9*pow(yy,4) + 4*pow(yy,2)*sqrt(pow(xx,2) + pow(yy,2)))))/pow(pow(xx,2) + pow(yy,2),2.5);
			pTensor[2] = 0;
		}

		if(pInputs->EMField_Cal_B)
		{
			pTensor[3] = 0;
			pTensor[4] = 0;
			pTensor[5] = (3*scaleB*xx*(1 + yy))/sqrt(pow(xx,2) + pow(yy,2));
		}
	}
	if(1 == Order)
	{
		pTensor[0] = (scaleE*(1 + pow(xx,4)*(-1 + sqrt(pow(xx,2) + pow(yy,2))) + 2*pow(xx,2)*pow(yy,2)*(-1 + sqrt(pow(xx,2) + pow(yy,2))) + pow(yy,4)*(-1 + sqrt(pow(xx,2) + pow(yy,2)))))/(pow(xx,2) + pow(yy,2));
		pTensor[1] = -((scaleB*xx*yy*(3 + yy))/sqrt(pow(xx,2) + pow(yy,2)));
		pTensor[2] = (scaleB*(pow(xx,2) - 3*yy)*yy)/sqrt(pow(xx,2) + pow(yy,2));
		pTensor[3] = 0;
	}
	if(2 == Order)
	{
		pTensor[0] = 0;
		pTensor[1] = 0;
		pTensor[2] = 0;
		pTensor[3] = 0;
		pTensor[4] = (scaleE*xx*(3*pow(xx,6) + 3*pow(yy,6) - 2*sqrt(pow(xx,2) + pow(yy,2)) - 2*pow(yy,4)*sqrt(pow(xx,2) + pow(yy,2)) + pow(xx,4)*(9*pow(yy,2) - 2*sqrt(pow(xx,2) + pow(yy,2))) + pow(xx,2)*(9*pow(yy,4) - 4*pow(yy,2)*sqrt(pow(xx,2) + pow(yy,2)))))/pow(pow(xx,2) + pow(yy,2),2.5);
		pTensor[5] = -((scaleB*pow(yy,3)*(3 + yy))/pow(pow(xx,2) + pow(yy,2),1.5));
		pTensor[6] = (scaleB*xx*yy*(pow(xx,2) + yy*(3 + 2*yy)))/pow(pow(xx,2) + pow(yy,2),1.5);
		pTensor[7] = 0;
		pTensor[8] = (scaleE*yy*(3*pow(xx,6) + 3*pow(yy,6) - 2*sqrt(pow(xx,2) + pow(yy,2)) - 2*pow(yy,4)*sqrt(pow(xx,2) + pow(yy,2)) + pow(xx,4)*(9*pow(yy,2) - 2*sqrt(pow(xx,2) + pow(yy,2))) + pow(xx,2)*(9*pow(yy,4) - 4*pow(yy,2)*sqrt(pow(xx,2) + pow(yy,2)))))/pow(pow(xx,2) + pow(yy,2),2.5);
		pTensor[9] = -((scaleB*xx*(pow(yy,3) + pow(xx,2)*(3 + 2*yy)))/pow(pow(xx,2) + pow(yy,2),1.5));
		pTensor[10] = (scaleB*(pow(xx,4) - 6*pow(xx,2)*yy - 3*pow(yy,3)))/pow(pow(xx,2) + pow(yy,2),1.5);
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
		pTensor[20] = -((scaleE*(-6*pow(xx,8) - 3*pow(yy,8) + 2*pow(yy,2)*sqrt(pow(xx,2) + pow(yy,2)) + 2*pow(yy,6)*sqrt(pow(xx,2) + pow(yy,2)) + pow(xx,6)*(-21*pow(yy,2) + 2*sqrt(pow(xx,2) + pow(yy,2))) + 3*pow(xx,4)*(-9*pow(yy,4) + 2*pow(yy,2)*sqrt(pow(xx,2) + pow(yy,2))) + 3*pow(xx,2)*(-5*pow(yy,6) - 2*sqrt(pow(xx,2) + pow(yy,2)) + 2*pow(yy,4)*sqrt(pow(xx,2) + pow(yy,2)))))/pow(pow(xx,2) + pow(yy,2),3.5));
		pTensor[21] = (3*scaleB*xx*pow(yy,3)*(3 + yy))/pow(pow(xx,2) + pow(yy,2),2.5);
		pTensor[22] = (scaleB*pow(yy,2)*(-(pow(xx,2)*(6 + yy)) + pow(yy,2)*(3 + 2*yy)))/pow(pow(xx,2) + pow(yy,2),2.5);
		pTensor[23] = 0;
		pTensor[24] = (scaleE*xx*yy*(3*pow(xx,6) + 9*pow(xx,4)*pow(yy,2) + 9*pow(xx,2)*pow(yy,4) + 3*pow(yy,6) + 8*sqrt(pow(xx,2) + pow(yy,2))))/pow(pow(xx,2) + pow(yy,2),3.5);
		pTensor[25] = -((scaleB*pow(yy,2)*(pow(yy,3) + pow(xx,2)*(9 + 4*yy)))/pow(pow(xx,2) + pow(yy,2),2.5));
		pTensor[26] = (scaleB*xx*(pow(xx,4) - 3*pow(yy,3) + 2*pow(xx,2)*yy*(3 + 2*yy)))/pow(pow(xx,2) + pow(yy,2),2.5);
		pTensor[27] = 0;
		pTensor[28] = 0;
		pTensor[29] = 0;
		pTensor[30] = 0;
		pTensor[31] = 0;
		pTensor[32] = 0;
		pTensor[33] = 0;
		pTensor[34] = 0;
		pTensor[35] = 0;
		pTensor[36] = (scaleE*xx*yy*(3*pow(xx,6) + 9*pow(xx,4)*pow(yy,2) + 9*pow(xx,2)*pow(yy,4) + 3*pow(yy,6) + 8*sqrt(pow(xx,2) + pow(yy,2))))/pow(pow(xx,2) + pow(yy,2),3.5);
		pTensor[37] = -((scaleB*pow(yy,2)*(pow(yy,3) + pow(xx,2)*(9 + 4*yy)))/pow(pow(xx,2) + pow(yy,2),2.5));
		pTensor[38] = (scaleB*xx*(pow(xx,4) - 3*pow(yy,3) + 2*pow(xx,2)*yy*(3 + 2*yy)))/pow(pow(xx,2) + pow(yy,2),2.5);
		pTensor[39] = 0;
		pTensor[40] = -((scaleE*(-3*pow(xx,8) - 6*pow(yy,8) - 6*pow(yy,2)*sqrt(pow(xx,2) + pow(yy,2)) + 2*pow(yy,6)*sqrt(pow(xx,2) + pow(yy,2)) + pow(xx,6)*(-15*pow(yy,2) + 2*sqrt(pow(xx,2) + pow(yy,2))) + 3*pow(xx,4)*(-9*pow(yy,4) + 2*pow(yy,2)*sqrt(pow(xx,2) + pow(yy,2))) + pow(xx,2)*(-21*pow(yy,6) + 2*sqrt(pow(xx,2) + pow(yy,2)) + 6*pow(yy,4)*sqrt(pow(xx,2) + pow(yy,2)))))/pow(pow(xx,2) + pow(yy,2),3.5));
		pTensor[41] = (scaleB*pow(xx,3)*(-2*pow(xx,2) + yy*(9 + yy)))/pow(pow(xx,2) + pow(yy,2),2.5);
		pTensor[42] = (3*scaleB*pow(xx,2)*(pow(yy,2) - pow(xx,2)*(2 + yy)))/pow(pow(xx,2) + pow(yy,2),2.5);
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
		fprintf(stderr,"ERROR: In function GAPS_APT_Field_Zjs_E_drift. This field function does NOT support tensors order larger than %d.\n",MaxOrder);
	}
	return 0;
}

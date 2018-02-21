#include "APT_AllHeaders.h"

int GAPS_APT_Field_EarthDipole(double *pTensor,double *pSpaceTime4,int Order,Gaps_IO_InputsContainer *pInputs)
{
	int MaxOrder = 3;
	double tt=pSpaceTime4[0],xx=pSpaceTime4[1],yy=pSpaceTime4[2],zz=pSpaceTime4[3];
	double R0 = pInputs->EMField_EarthDipole_R0;
	double B0=pInputs->EMField_B0;
	double E0=pInputs->EMField_E0;
	if(-1 == Order)
	{
		if(pInputs->EMField_Cal_E)
		{
			pTensor[0] = 0;
			pTensor[1] = 0;
			pTensor[2] = 0;
		}

		if(pInputs->EMField_Cal_B)
		{
			pTensor[3] = (B0*pow(R0,3)*xx*zz)/(sqrt(pow(xx,2) + pow(yy,2))*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),2));
			pTensor[4] = (B0*pow(R0,3)*yy*zz)/(sqrt(pow(xx,2) + pow(yy,2))*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),2));
			pTensor[5] = -(B0*pow(R0,3)*(pow(xx,2) + pow(yy,2) - pow(zz,2)))/(2.*sqrt(pow(xx,2) + pow(yy,2))*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),2));
		}
	}
	if(1 == Order)
	{
		pTensor[0] = 0;
		pTensor[1] = -(B0*pow(R0,3)*yy)/(2.*sqrt(pow(xx,2) + pow(yy,2))*(pow(xx,2) + pow(yy,2) + pow(zz,2)));
		pTensor[2] = (B0*pow(R0,3)*xx)/(2.*sqrt(pow(xx,2) + pow(yy,2))*(pow(xx,2) + pow(yy,2) + pow(zz,2)));
		pTensor[3] = 0;
	}
	if(2 == Order)
	{
		pTensor[0] = 0;
		pTensor[1] = 0;
		pTensor[2] = 0;
		pTensor[3] = 0;
		pTensor[4] = 0;
		pTensor[5] = (B0*pow(R0,3)*xx*yy*(3*pow(xx,2) + 3*pow(yy,2) + pow(zz,2)))/(2.*pow(pow(xx,2) + pow(yy,2),1.5)*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),2));
		pTensor[6] = -(B0*pow(R0,3)*(2*pow(xx,4) + pow(xx,2)*pow(yy,2) - pow(yy,2)*(pow(yy,2) + pow(zz,2))))/(2.*pow(pow(xx,2) + pow(yy,2),1.5)*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),2));
		pTensor[7] = 0;
		pTensor[8] = 0;
		pTensor[9] = -(B0*pow(R0,3)*(pow(xx,4) - 2*pow(yy,4) + pow(xx,2)*(-pow(yy,2) + pow(zz,2))))/(2.*pow(pow(xx,2) + pow(yy,2),1.5)*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),2));
		pTensor[10] = -(B0*pow(R0,3)*xx*yy*(3*pow(xx,2) + 3*pow(yy,2) + pow(zz,2)))/(2.*pow(pow(xx,2) + pow(yy,2),1.5)*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),2));
		pTensor[11] = 0;
		pTensor[12] = 0;
		pTensor[13] = (B0*pow(R0,3)*yy*zz)/(sqrt(pow(xx,2) + pow(yy,2))*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),2));
		pTensor[14] = -((B0*pow(R0,3)*xx*zz)/(sqrt(pow(xx,2) + pow(yy,2))*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),2)));
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
		pTensor[21] = (B0*pow(R0,3)*yy*(-12*pow(xx,6) - 3*pow(xx,4)*(7*pow(yy,2) + 2*pow(zz,2)) - 2*pow(xx,2)*(3*pow(yy,4) + pow(yy,2)*pow(zz,2) + pow(zz,4)) + pow(yy,2)*(3*pow(yy,4) + 4*pow(yy,2)*pow(zz,2) + pow(zz,4))))/(2.*pow(pow(xx,2) + pow(yy,2),2.5)*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),3));
		pTensor[22] = (B0*pow(R0,3)*xx*(6*pow(xx,6) + pow(xx,4)*(3*pow(yy,2) - 2*pow(zz,2)) - 2*pow(xx,2)*(6*pow(yy,4) + 7*pow(yy,2)*pow(zz,2)) - 3*pow(yy,2)*(3*pow(yy,4) + 4*pow(yy,2)*pow(zz,2) + pow(zz,4))))/(2.*pow(pow(xx,2) + pow(yy,2),2.5)*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),3));
		pTensor[23] = 0;
		pTensor[24] = 0;
		pTensor[25] = (B0*pow(R0,3)*xx*(3*pow(xx,6) + pow(xx,4)*(-6*pow(yy,2) + 4*pow(zz,2)) + pow(xx,2)*(-21*pow(yy,4) - 2*pow(yy,2)*pow(zz,2) + pow(zz,4)) - 2*pow(yy,2)*(6*pow(yy,4) + 3*pow(yy,2)*pow(zz,2) + pow(zz,4))))/(2.*pow(pow(xx,2) + pow(yy,2),2.5)*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),3));
		pTensor[26] = -(B0*pow(R0,3)*yy*(-12*pow(xx,6) - 3*pow(xx,4)*(7*pow(yy,2) + 2*pow(zz,2)) - 2*pow(xx,2)*(3*pow(yy,4) + pow(yy,2)*pow(zz,2) + pow(zz,4)) + pow(yy,2)*(3*pow(yy,4) + 4*pow(yy,2)*pow(zz,2) + pow(zz,4))))/(2.*pow(pow(xx,2) + pow(yy,2),2.5)*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),3));
		pTensor[27] = 0;
		pTensor[28] = 0;
		pTensor[29] = -((B0*pow(R0,3)*xx*yy*zz*(5*pow(xx,2) + 5*pow(yy,2) + pow(zz,2)))/(pow(pow(xx,2) + pow(yy,2),1.5)*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),3)));
		pTensor[30] = -((B0*pow(R0,3)*zz*(-4*pow(xx,4) - 3*pow(xx,2)*pow(yy,2) + pow(yy,4) + pow(yy,2)*pow(zz,2)))/(pow(pow(xx,2) + pow(yy,2),1.5)*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),3)));
		pTensor[31] = 0;
		pTensor[32] = 0;
		pTensor[33] = 0;
		pTensor[34] = 0;
		pTensor[35] = 0;
		pTensor[36] = 0;
		pTensor[37] = (B0*pow(R0,3)*xx*(3*pow(xx,6) + pow(xx,4)*(-6*pow(yy,2) + 4*pow(zz,2)) + pow(xx,2)*(-21*pow(yy,4) - 2*pow(yy,2)*pow(zz,2) + pow(zz,4)) - 2*pow(yy,2)*(6*pow(yy,4) + 3*pow(yy,2)*pow(zz,2) + pow(zz,4))))/(2.*pow(pow(xx,2) + pow(yy,2),2.5)*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),3));
		pTensor[38] = -(B0*pow(R0,3)*yy*(-12*pow(xx,6) - 3*pow(xx,4)*(7*pow(yy,2) + 2*pow(zz,2)) - 2*pow(xx,2)*(3*pow(yy,4) + pow(yy,2)*pow(zz,2) + pow(zz,4)) + pow(yy,2)*(3*pow(yy,4) + 4*pow(yy,2)*pow(zz,2) + pow(zz,4))))/(2.*pow(pow(xx,2) + pow(yy,2),2.5)*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),3));
		pTensor[39] = 0;
		pTensor[40] = 0;
		pTensor[41] = (B0*pow(R0,3)*yy*(9*pow(xx,6) - 6*pow(yy,6) + 2*pow(yy,4)*pow(zz,2) + 12*pow(xx,4)*(pow(yy,2) + pow(zz,2)) + pow(xx,2)*(-3*pow(yy,4) + 14*pow(yy,2)*pow(zz,2) + 3*pow(zz,4))))/(2.*pow(pow(xx,2) + pow(yy,2),2.5)*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),3));
		pTensor[42] = -(B0*pow(R0,3)*xx*(3*pow(xx,6) + pow(xx,4)*(-6*pow(yy,2) + 4*pow(zz,2)) + pow(xx,2)*(-21*pow(yy,4) - 2*pow(yy,2)*pow(zz,2) + pow(zz,4)) - 2*pow(yy,2)*(6*pow(yy,4) + 3*pow(yy,2)*pow(zz,2) + pow(zz,4))))/(2.*pow(pow(xx,2) + pow(yy,2),2.5)*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),3));
		pTensor[43] = 0;
		pTensor[44] = 0;
		pTensor[45] = (B0*pow(R0,3)*zz*(pow(xx,4) - 4*pow(yy,4) + pow(xx,2)*(-3*pow(yy,2) + pow(zz,2))))/(pow(pow(xx,2) + pow(yy,2),1.5)*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),3));
		pTensor[46] = (B0*pow(R0,3)*xx*yy*zz*(5*pow(xx,2) + 5*pow(yy,2) + pow(zz,2)))/(pow(pow(xx,2) + pow(yy,2),1.5)*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),3));
		pTensor[47] = 0;
		pTensor[48] = 0;
		pTensor[49] = 0;
		pTensor[50] = 0;
		pTensor[51] = 0;
		pTensor[52] = 0;
		pTensor[53] = -((B0*pow(R0,3)*xx*yy*zz*(5*pow(xx,2) + 5*pow(yy,2) + pow(zz,2)))/(pow(pow(xx,2) + pow(yy,2),1.5)*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),3)));
		pTensor[54] = -((B0*pow(R0,3)*zz*(-4*pow(xx,4) - 3*pow(xx,2)*pow(yy,2) + pow(yy,4) + pow(yy,2)*pow(zz,2)))/(pow(pow(xx,2) + pow(yy,2),1.5)*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),3)));
		pTensor[55] = 0;
		pTensor[56] = 0;
		pTensor[57] = (B0*pow(R0,3)*zz*(pow(xx,4) - 4*pow(yy,4) + pow(xx,2)*(-3*pow(yy,2) + pow(zz,2))))/(pow(pow(xx,2) + pow(yy,2),1.5)*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),3));
		pTensor[58] = (B0*pow(R0,3)*xx*yy*zz*(5*pow(xx,2) + 5*pow(yy,2) + pow(zz,2)))/(pow(pow(xx,2) + pow(yy,2),1.5)*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),3));
		pTensor[59] = 0;
		pTensor[60] = 0;
		pTensor[61] = (B0*pow(R0,3)*yy*(pow(xx,2) + pow(yy,2) - 3*pow(zz,2)))/(sqrt(pow(xx,2) + pow(yy,2))*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),3));
		pTensor[62] = -((B0*pow(R0,3)*xx*(pow(xx,2) + pow(yy,2) - 3*pow(zz,2)))/(sqrt(pow(xx,2) + pow(yy,2))*pow(pow(xx,2) + pow(yy,2) + pow(zz,2),3)));
		pTensor[63] = 0;
	}
	if(MaxOrder<Order)
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_Field_EarthDipole. This field function does NOT support tensors order larger than %d.\n",MaxOrder);
	}
	return 0;
}
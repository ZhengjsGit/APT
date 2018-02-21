#include "APT_AllHeaders.h"

int GAPS_APT_Field_Tokamak(double *pTensor,double *pSpaceTime4,int Order,Gaps_IO_InputsContainer *pInputs)
{
	int MaxOrder = 3;
	double tt=pSpaceTime4[0],xx=pSpaceTime4[1],yy=pSpaceTime4[2],zz=pSpaceTime4[3];
	double q = pInputs->EMField_Tokamak_q;
	double a = pInputs->EMField_Tokamak_a;
	double R0 = pInputs->EMField_Tokamak_R0;
	double B0=pInputs->EMField_B0;
	double E0=pInputs->EMField_E0;
	if(-1 == Order)
	{
		if(pInputs->EMField_Cal_E)
		{
			pTensor[0] = -((E0*R0*yy)/(pow(xx,2) + pow(yy,2)));
			pTensor[1] = (E0*R0*xx)/(pow(xx,2) + pow(yy,2));
			pTensor[2] = 0;
		}

		if(pInputs->EMField_Cal_B)
		{
			pTensor[3] = (B0*(-(q*R0*yy) + xx*zz))/(q*(pow(xx,2) + pow(yy,2)));
			pTensor[4] = (B0*(q*R0*xx + yy*zz))/(q*(pow(xx,2) + pow(yy,2)));
			pTensor[5] = (B0*(-1 + R0/sqrt(pow(xx,2) + pow(yy,2))))/q;
		}
	}
	if(1 == Order)
	{
		pTensor[0] = 0;
		pTensor[1] = (2*E0*q*R0*tt*yy + B0*(pow(R0,2)*yy + R0*(-2*yy*sqrt(pow(xx,2) + pow(yy,2)) + q*xx*zz) + yy*(pow(xx,2) + pow(yy,2) + pow(zz,2))))/(2.*q*(pow(xx,2) + pow(yy,2)));
		pTensor[2] = (-2*E0*q*R0*tt*xx + B0*(-(pow(R0,2)*xx) + R0*(2*xx*sqrt(pow(xx,2) + pow(yy,2)) + q*yy*zz) - xx*(pow(xx,2) + pow(yy,2) + pow(zz,2))))/(2.*q*(pow(xx,2) + pow(yy,2)));
		pTensor[3] = -(B0*R0*log(sqrt(pow(xx,2) + pow(yy,2))/R0))/2.;
	}
	if(2 == Order)
	{
		pTensor[0] = 0;
		pTensor[1] = (E0*R0*yy)/(pow(xx,2) + pow(yy,2));
		pTensor[2] = -((E0*R0*xx)/(pow(xx,2) + pow(yy,2)));
		pTensor[3] = 0;
		pTensor[4] = 0;
		pTensor[5] = (-4*E0*q*R0*tt*xx*yy + B0*(-2*pow(R0,2)*xx*yy - 2*xx*yy*pow(zz,2) + R0*(2*xx*yy*sqrt(pow(xx,2) + pow(yy,2)) - q*pow(xx,2)*zz + q*pow(yy,2)*zz)))/(2.*q*pow(pow(xx,2) + pow(yy,2),2));
		pTensor[6] = (2*E0*q*R0*tt*(pow(xx,2) - pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2)) - B0*(pow(R0,2)*(-pow(xx,2) + pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2)) - 2*R0*yy*(pow(xx,2)*yy + pow(yy,3) - q*xx*sqrt(pow(xx,2) + pow(yy,2))*zz) + sqrt(pow(xx,2) + pow(yy,2))*(pow(xx,4) + 2*pow(xx,2)*pow(yy,2) + pow(yy,4) - pow(xx,2)*pow(zz,2) + pow(yy,2)*pow(zz,2))))/(2.*q*pow(pow(xx,2) + pow(yy,2),2.5));
		pTensor[7] = -(B0*R0*xx)/(2.*(pow(xx,2) + pow(yy,2)));
		pTensor[8] = 0;
		pTensor[9] = (2*E0*q*R0*tt*(pow(xx,2) - pow(yy,2)) + B0*(pow(xx,4) + pow(yy,4) + pow(R0,2)*(pow(xx,2) - pow(yy,2)) - pow(yy,2)*pow(zz,2) - 2*R0*xx*(xx*sqrt(pow(xx,2) + pow(yy,2)) + q*yy*zz) + pow(xx,2)*(2*pow(yy,2) + pow(zz,2))))/(2.*q*pow(pow(xx,2) + pow(yy,2),2));
		pTensor[10] = (4*E0*q*R0*tt*xx*yy + B0*(2*pow(R0,2)*xx*yy + 2*xx*yy*pow(zz,2) + R0*(-2*xx*yy*sqrt(pow(xx,2) + pow(yy,2)) + q*pow(xx,2)*zz - q*pow(yy,2)*zz)))/(2.*q*pow(pow(xx,2) + pow(yy,2),2));
		pTensor[11] = -(B0*R0*yy)/(2.*(pow(xx,2) + pow(yy,2)));
		pTensor[12] = 0;
		pTensor[13] = (B0*q*R0*xx + 2*B0*yy*zz)/(2*q*pow(xx,2) + 2*q*pow(yy,2));
		pTensor[14] = (B0*q*R0*yy - 2*B0*xx*zz)/(2*q*pow(xx,2) + 2*q*pow(yy,2));
		pTensor[15] = 0;
	}
	if(3 == Order)
	{
		pTensor[0] = 0;
		pTensor[1] = 0;
		pTensor[2] = 0;
		pTensor[3] = 0;
		pTensor[4] = 0;
		pTensor[5] = (-2*E0*R0*xx*yy)/pow(pow(xx,2) + pow(yy,2),2);
		pTensor[6] = (E0*R0*(pow(xx,2) - pow(yy,2)))/pow(pow(xx,2) + pow(yy,2),2);
		pTensor[7] = 0;
		pTensor[8] = 0;
		pTensor[9] = (E0*R0*(pow(xx,2) - pow(yy,2)))/pow(pow(xx,2) + pow(yy,2),2);
		pTensor[10] = (2*E0*R0*xx*yy)/pow(pow(xx,2) + pow(yy,2),2);
		pTensor[11] = 0;
		pTensor[12] = 0;
		pTensor[13] = 0;
		pTensor[14] = 0;
		pTensor[15] = 0;
		pTensor[16] = 0;
		pTensor[17] = (-2*E0*R0*xx*yy)/pow(pow(xx,2) + pow(yy,2),2);
		pTensor[18] = (E0*R0*(pow(xx,2) - pow(yy,2)))/pow(pow(xx,2) + pow(yy,2),2);
		pTensor[19] = 0;
		pTensor[20] = 0;
		pTensor[21] = (-2*E0*q*R0*tt*yy*(-3*pow(xx,2) + pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2)) + B0*(-(pow(R0,2)*yy*(-3*pow(xx,2) + pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2))) - yy*(-3*pow(xx,2) + pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2))*pow(zz,2) + R0*(-2*pow(xx,4)*yy - pow(xx,2)*pow(yy,3) + pow(yy,5) + q*pow(xx,3)*sqrt(pow(xx,2) + pow(yy,2))*zz - 3*q*xx*pow(yy,2)*sqrt(pow(xx,2) + pow(yy,2))*zz)))/(q*pow(pow(xx,2) + pow(yy,2),3.5));
		pTensor[22] = (-2*E0*q*R0*tt*xx*(pow(xx,2) - 3*pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2)) - B0*(pow(R0,2)*xx*(pow(xx,2) - 3*pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2)) + xx*(pow(xx,2) - 3*pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2))*pow(zz,2) + R0*yy*(3*pow(xx,3)*yy + 3*xx*pow(yy,3) - 3*q*pow(xx,2)*sqrt(pow(xx,2) + pow(yy,2))*zz + q*pow(yy,2)*sqrt(pow(xx,2) + pow(yy,2))*zz)))/(q*pow(pow(xx,2) + pow(yy,2),3.5));
		pTensor[23] = (B0*R0*(pow(xx,2) - pow(yy,2)))/(2.*pow(pow(xx,2) + pow(yy,2),2));
		pTensor[24] = 0;
		pTensor[25] = (-2*E0*q*R0*tt*xx*(pow(xx,2) - 3*pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2)) - B0*(pow(R0,2)*xx*(pow(xx,2) - 3*pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2)) + xx*(pow(xx,2) - 3*pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2))*pow(zz,2) + R0*(-pow(xx,5) + pow(xx,3)*pow(yy,2) + 2*xx*pow(yy,4) - 3*q*pow(xx,2)*yy*sqrt(pow(xx,2) + pow(yy,2))*zz + q*pow(yy,3)*sqrt(pow(xx,2) + pow(yy,2))*zz)))/(q*pow(pow(xx,2) + pow(yy,2),3.5));
		pTensor[26] = (2*E0*q*R0*tt*yy*(-3*pow(xx,2) + pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2)) + B0*(pow(R0,2)*yy*(-3*pow(xx,2) + pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2)) + yy*(-3*pow(xx,2) + pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2))*pow(zz,2) + R0*(2*pow(xx,4)*yy + pow(xx,2)*pow(yy,3) - pow(yy,5) - q*pow(xx,3)*sqrt(pow(xx,2) + pow(yy,2))*zz + 3*q*xx*pow(yy,2)*sqrt(pow(xx,2) + pow(yy,2))*zz)))/(q*pow(pow(xx,2) + pow(yy,2),3.5));
		pTensor[27] = (B0*R0*xx*yy)/pow(pow(xx,2) + pow(yy,2),2);
		pTensor[28] = 0;
		pTensor[29] = -(B0*(q*R0*(pow(xx,2) - pow(yy,2)) + 4*xx*yy*zz))/(2.*q*pow(pow(xx,2) + pow(yy,2),2));
		pTensor[30] = -((B0*(q*R0*xx*yy + (-pow(xx,2) + pow(yy,2))*zz))/(q*pow(pow(xx,2) + pow(yy,2),2)));
		pTensor[31] = 0;
		pTensor[32] = 0;
		pTensor[33] = (E0*R0*(pow(xx,2) - pow(yy,2)))/pow(pow(xx,2) + pow(yy,2),2);
		pTensor[34] = (2*E0*R0*xx*yy)/pow(pow(xx,2) + pow(yy,2),2);
		pTensor[35] = 0;
		pTensor[36] = 0;
		pTensor[37] = (-2*E0*q*R0*tt*xx*(pow(xx,2) - 3*pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2)) - B0*(pow(R0,2)*xx*(pow(xx,2) - 3*pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2)) + xx*(pow(xx,2) - 3*pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2))*pow(zz,2) + R0*(-pow(xx,5) + pow(xx,3)*pow(yy,2) + 2*xx*pow(yy,4) - 3*q*pow(xx,2)*yy*sqrt(pow(xx,2) + pow(yy,2))*zz + q*pow(yy,3)*sqrt(pow(xx,2) + pow(yy,2))*zz)))/(q*pow(pow(xx,2) + pow(yy,2),3.5));
		pTensor[38] = (2*E0*q*R0*tt*yy*(-3*pow(xx,2) + pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2)) + B0*(pow(R0,2)*yy*(-3*pow(xx,2) + pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2)) + yy*(-3*pow(xx,2) + pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2))*pow(zz,2) + R0*(2*pow(xx,4)*yy + pow(xx,2)*pow(yy,3) - pow(yy,5) - q*pow(xx,3)*sqrt(pow(xx,2) + pow(yy,2))*zz + 3*q*xx*pow(yy,2)*sqrt(pow(xx,2) + pow(yy,2))*zz)))/(q*pow(pow(xx,2) + pow(yy,2),3.5));
		pTensor[39] = (B0*R0*xx*yy)/pow(pow(xx,2) + pow(yy,2),2);
		pTensor[40] = 0;
		pTensor[41] = (2*E0*q*R0*tt*yy*(-3*pow(xx,2) + pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2)) + B0*(pow(R0,2)*yy*(-3*pow(xx,2) + pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2)) + yy*(-3*pow(xx,2) + pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2))*pow(zz,2) + R0*xx*(3*pow(xx,3)*yy + 3*xx*pow(yy,3) - q*pow(xx,2)*sqrt(pow(xx,2) + pow(yy,2))*zz + 3*q*pow(yy,2)*sqrt(pow(xx,2) + pow(yy,2))*zz)))/(q*pow(pow(xx,2) + pow(yy,2),3.5));
		pTensor[42] = (2*E0*q*R0*tt*xx*(pow(xx,2) - 3*pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2)) + B0*(pow(R0,2)*xx*(pow(xx,2) - 3*pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2)) + xx*(pow(xx,2) - 3*pow(yy,2))*sqrt(pow(xx,2) + pow(yy,2))*pow(zz,2) + R0*(-pow(xx,5) + pow(xx,3)*pow(yy,2) + 2*xx*pow(yy,4) - 3*q*pow(xx,2)*yy*sqrt(pow(xx,2) + pow(yy,2))*zz + q*pow(yy,3)*sqrt(pow(xx,2) + pow(yy,2))*zz)))/(q*pow(pow(xx,2) + pow(yy,2),3.5));
		pTensor[43] = (B0*R0*(-pow(xx,2) + pow(yy,2)))/(2.*pow(pow(xx,2) + pow(yy,2),2));
		pTensor[44] = 0;
		pTensor[45] = -((B0*(q*R0*xx*yy + (-pow(xx,2) + pow(yy,2))*zz))/(q*pow(pow(xx,2) + pow(yy,2),2)));
		pTensor[46] = (B0*(q*R0*(pow(xx,2) - pow(yy,2)) + 4*xx*yy*zz))/(2.*q*pow(pow(xx,2) + pow(yy,2),2));
		pTensor[47] = 0;
		pTensor[48] = 0;
		pTensor[49] = 0;
		pTensor[50] = 0;
		pTensor[51] = 0;
		pTensor[52] = 0;
		pTensor[53] = -(B0*(q*R0*(pow(xx,2) - pow(yy,2)) + 4*xx*yy*zz))/(2.*q*pow(pow(xx,2) + pow(yy,2),2));
		pTensor[54] = -((B0*(q*R0*xx*yy + (-pow(xx,2) + pow(yy,2))*zz))/(q*pow(pow(xx,2) + pow(yy,2),2)));
		pTensor[55] = 0;
		pTensor[56] = 0;
		pTensor[57] = -((B0*(q*R0*xx*yy + (-pow(xx,2) + pow(yy,2))*zz))/(q*pow(pow(xx,2) + pow(yy,2),2)));
		pTensor[58] = (B0*(q*R0*(pow(xx,2) - pow(yy,2)) + 4*xx*yy*zz))/(2.*q*pow(pow(xx,2) + pow(yy,2),2));
		pTensor[59] = 0;
		pTensor[60] = 0;
		pTensor[61] = (B0*yy)/(q*(pow(xx,2) + pow(yy,2)));
		pTensor[62] = -((B0*xx)/(q*(pow(xx,2) + pow(yy,2))));
		pTensor[63] = 0;
	}
	if(MaxOrder<Order)
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_Field_Tokamak. This field function does NOT support tensors order larger than %d.\n",MaxOrder);
	}
	return 0;
}

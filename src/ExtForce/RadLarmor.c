#include "APT_AllHeaders.h"

int GAPS_APT_ExtForce_RadLarmor(double *pForce3,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	/**Step 1: Get pointers of particle data**/
	double *pP=GAPS_APT_GetP3(pPtc);
	double *pAclr=GAPS_APT_GetAclr3(pPtc);
	double *pGamma=GAPS_APT_GetGamma1(pPtc);
	/**Step 2: Get constant coefficients from pInputs**/
	double cons=pInputs->ExtForce_RadLarmor_Const;
	/**Step 3: Calculate and assign values to pForce3[0]~pForce3[2]**/
	double tempV[3];
	double pV[3];
	V3Smult(1/(*pGamma),pP,pV);
	double Ps;
	V3cross(pV,pAclr,tempV);
	Ps = cons * pow(*pGamma,6) * (V3INmult(pAclr,pAclr)-V3INmult(tempV,tempV));
	V3Smult(-1*Ps/(V3INmult(pV,pV)),pV,pForce3);
	return 0;
}

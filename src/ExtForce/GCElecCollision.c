#include "APT_AllHeaders.h"

int GAPS_APT_ExtForce_GCElecCollision(double *pForce3,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	/**Step 1: Get pointers of particle data**/
	double *pP=GAPS_APT_GetP3(pPtc);
	double *pGamma=GAPS_APT_GetGamma1(pPtc);
	/**Step 2: Get constant coefficients from pInputs**/
	double cons=pInputs->ExtForce_GCElecCollision_Const;
	double ne=pInputs->ExtForce_GCElecCollision_Ne;
	/**Step 3: Calculate and assign values to pForce3[0]~pForce3[2]**/
	double normP;
	normP = sqrt(pP[0]*pP[0]+pP[1]*pP[1]+pP[2]*pP[2]);

	V3Smult(-1* cons* pGamma[0]*pGamma[0]/(pow(normP,3)),pP,pForce3);
	return 0;
}

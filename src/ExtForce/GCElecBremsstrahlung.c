#include "APT_AllHeaders.h"

int GAPS_APT_ExtForce_GCElecBremsstrahlung(double *pForce3,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pP=GAPS_APT_GetP3(pPtc);
	double *pGamma=GAPS_APT_GetGamma1(pPtc);
	double cons=pInputs->ExtForce_GCElecBremsstrahlung_Const;
	double ne=pInputs->ExtForce_GCElecBremsstrahlung_Ne;
	double Zeff=pInputs->ExtForce_GCElecBremsstrahlung_Zeff;
	double normP;
	normP = sqrt(pP[0]*pP[0]+pP[1]*pP[1]+pP[2]*pP[2]);

	V3Smult(-1* cons* (*pGamma)*(log((*pGamma)*2)-1/3)/normP,pP,pForce3);
	return 0;
}


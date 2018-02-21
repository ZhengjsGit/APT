#include "APT_AllHeaders.h"

int GAPS_APT_SetParticlePosition_ParabolicTorus(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pX=GAPS_APT_GetX3(pPtc);
	double R0=pInputs->Init_X_ParabolicTorus_MajorRadius;
	double r_max=pInputs->Init_X_ParabolicTorus_rmax;
	double Unit_Space=pInputs->Unit_Space;

	double x0[3]={GenRandNum_Para(r_max),GenRandNum_Uniform(0,2*M_PI),GenRandNum_Uniform(0,2*M_PI)};
	TORD2CART(x0,pX,R0);
	return 0;
}


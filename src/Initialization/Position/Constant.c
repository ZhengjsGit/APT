#include "APT_AllHeaders.h"

int GAPS_APT_SetParticlePosition_Constant(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pX=GAPS_APT_GetX3(pPtc);
	pX[0]=pInputs->Init_X_Constant_X0[0];
	pX[1]=pInputs->Init_X_Constant_X0[1];
	pX[2]=pInputs->Init_X_Constant_X0[2];
	return 0;
}

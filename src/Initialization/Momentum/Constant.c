#include "APT_AllHeaders.h"

int GAPS_APT_SetParticleMomentum_Constant(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pP=GAPS_APT_GetP3(pPtc);
	pP[0]=pInputs->Init_P_Constant_P0[0];
	pP[1]=pInputs->Init_P_Constant_P0[1];
	pP[2]=pInputs->Init_P_Constant_P0[2];
	return 0;
}

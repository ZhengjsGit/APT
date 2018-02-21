#include "APT_AllHeaders.h"

int GAPS_APT_SetParticleAcceleration_Constant(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pAclr=GAPS_APT_GetAclr3(pPtc);
	pAclr[0]=pInputs->Init_Aclr_Constant_Aclr0[0];
	pAclr[1]=pInputs->Init_Aclr_Constant_Aclr0[1];
	pAclr[2]=pInputs->Init_Aclr_Constant_Aclr0[2];
	return 0;
}

#include "APT_AllHeaders.h"

int GAPS_APT_SetParticleMomentum_Maxwell(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pP=GAPS_APT_GetP3(pPtc);
	double v_disp[3];
	double sigma=sqrt(pInputs->Init_P_Maxwell_Temp);
	double mu=0.;

	pP[0]=GenRandNum_Norm(mu,sigma);
	pP[1]=GenRandNum_Norm(mu,sigma);
	pP[2]=GenRandNum_Norm(mu,sigma);
	return 0;
}

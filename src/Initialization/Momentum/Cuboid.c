#include "APT_AllHeaders.h"

int GAPS_APT_SetParticleMomentum_Cuboid(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	//fix vx and vz calc vy
	double *pP=GAPS_APT_GetP3(pPtc);
	double vx_min=pInputs->Init_P_Cuboid_Boundaries[0];
	double vx_max=pInputs->Init_P_Cuboid_Boundaries[1];
	double vy_min=pInputs->Init_P_Cuboid_Boundaries[2];
	double vy_max=pInputs->Init_P_Cuboid_Boundaries[3];
	double vz_min=pInputs->Init_P_Cuboid_Boundaries[4];
	double vz_max=pInputs->Init_P_Cuboid_Boundaries[5];
	double E_k=pInputs->Init_P_Cuboid_E_k;

	pP[0]=GenRandNum_Uniform(vx_min,vx_max);
	pP[2]=GenRandNum_Uniform(vz_min,vz_max);
	pP[1]=sqrt(2*E_k - pP[0]*pP[0] -pP[2]*pP[2]);
	return 0;
}

#include "APT_AllHeaders.h"

int GAPS_APT_SetParticlePosition_Cuboid(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pX=GAPS_APT_GetX3(pPtc);
	double x_min=pInputs->Init_X_Cuboid_Boundaries[0];
	double x_max=pInputs->Init_X_Cuboid_Boundaries[1];
	double y_min=pInputs->Init_X_Cuboid_Boundaries[2];
	double y_max=pInputs->Init_X_Cuboid_Boundaries[3];
	double z_min=pInputs->Init_X_Cuboid_Boundaries[4];
	double z_max=pInputs->Init_X_Cuboid_Boundaries[5];

	pX[0]=GenRandNum_Uniform(x_min,x_max);
	pX[1]=GenRandNum_Uniform(y_min,y_max);
	pX[2]=GenRandNum_Uniform(z_min,z_max);
	return 0;
}

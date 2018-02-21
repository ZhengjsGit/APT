#include "APT_AllHeaders.h"

int GAPS_APT_SetParticlePosition_Cylinder(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pX=GAPS_APT_GetX3(pPtc);
	double R_min	=pInputs->Init_X_Cylinder_Boundaries[0];
	double R_max	=pInputs->Init_X_Cylinder_Boundaries[1];
	double Theta_min=pInputs->Init_X_Cylinder_Boundaries[2];
	double Theta_max=pInputs->Init_X_Cylinder_Boundaries[3];
	double Z_min	=pInputs->Init_X_Cylinder_Boundaries[4];
	double Z_max	=pInputs->Init_X_Cylinder_Boundaries[5];

	double x0[3]={GenRandNum_Circ(R_min,R_max),GenRandNum_Uniform(Theta_min,Theta_max),GenRandNum_Uniform(Z_min,Z_max)};
	CYLD2CART(x0,pX);
	return 0;
}

#include "APT_AllHeaders.h"

int GAPS_APT_SetParticlePosition_Torus(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pX=GAPS_APT_GetX3(pPtc);
	double r_min	=pInputs->Init_X_Torus_Boundaries[0];
	double r_max	=pInputs->Init_X_Torus_Boundaries[1];
	double theta_min=pInputs->Init_X_Torus_Boundaries[2];
	double theta_max=pInputs->Init_X_Torus_Boundaries[3];
	double xi_min	=pInputs->Init_X_Torus_Boundaries[4];
	double xi_max	=pInputs->Init_X_Torus_Boundaries[5];
	double R0=pInputs->Init_X_Torus_MajorRadius;

	double x0[3]={GenRandNum_Circ(r_min,r_max),GenRandNum_Uniform(theta_min,theta_max),GenRandNum_Uniform(xi_min,xi_max)};
	TORD2CART(x0,pX,R0);
	return 0;
}

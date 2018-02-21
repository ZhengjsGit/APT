#include "APT_AllHeaders.h"

int GAPS_APT_SetParticleMomentum_Gyrocenter(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pST=GAPS_APT_GetX4(pPtc);
	double *pP=GAPS_APT_GetP3(pPtc);
	double *pGamma=GAPS_APT_GetGamma1(pPtc);

	double EB[6],*B=EB+3;
	(pPtc->FieldFunc)(EB,pST,-1,pInputs);

	double GammaMu			=pInputs->Init_P_Gyrocenter_SampleRegion[0];
	double GammaSigma		=pInputs->Init_P_Gyrocenter_SampleRegion[1];
	double PitchAngle_min	=pInputs->Init_P_Gyrocenter_SampleRegion[2];
	double PitchAngle_max	=pInputs->Init_P_Gyrocenter_SampleRegion[3];
	double Gyrophase_min	=pInputs->Init_P_Gyrocenter_SampleRegion[4];
	double Gyrophase_max	=pInputs->Init_P_Gyrocenter_SampleRegion[5];

	*pGamma = GenRandNum_Norm(GammaMu,GammaSigma);

	double p_perpOVERp=GenRandNum_Uniform(PitchAngle_min,PitchAngle_max);
	double normP=sqrt(pow((*pGamma),2)-1);
	double p0[3]={sqrt(pow(normP,2)-pow(p_perpOVERp*normP,2)) , p_perpOVERp*normP, GenRandNum_Uniform(Gyrophase_min,Gyrophase_max)};

	p_GC_TO_p_CART(p0,B,pP);
	return 0;
}

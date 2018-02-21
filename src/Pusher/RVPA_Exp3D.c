#include "APT_AllHeaders.h"
inline int p_minus2p_plus_exp(double dT,double *B,double gamma, double *Pp,double qOverM);

int GAPS_APT_Pusher_RVPA_Exp3D (Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	/**Step 1: Get pointers of particle data**/
	double *pT = GAPS_APT_GetT1(pPtc);
	double *pX = GAPS_APT_GetX3(pPtc);
	double *pGamma=GAPS_APT_GetGamma1(pPtc);
	double *pP = GAPS_APT_GetP3(pPtc);
	double *E = GAPS_APT_GetE3(pPtc);
	double *B = GAPS_APT_GetB3(pPtc);

	double *pMass = GAPS_APT_GetMass1(pPtc);
	double *pCharge= GAPS_APT_GetCharge1(pPtc);
	/**Step 2: Get parameters from pInputs**/
	double dT=pInputs->dT;
	
	/**Step 3: Update pT, pX, pP, E, B through an algorithm**/
	double F_ext[3],E_eff[3];
	//Update and Return E,B
	GAPS_APT_CalEB(E,B,pPtc,pInputs);
	
	//Update and Return sum of all extern forces
	GAPS_APT_MergeExtForce(F_ext,pPtc,pInputs);

	E_eff[0]=pCharge[0]*E[0]+F_ext[0];
	E_eff[1]=pCharge[0]*E[1]+F_ext[1];
	E_eff[2]=pCharge[0]*E[2]+F_ext[2];

	// Core of algorithm
	double dT_half=0.5*dT;	
	//p_minus
	pP[0]+=dT_half*E_eff[0];
	pP[1]+=dT_half*E_eff[1];
	pP[2]+=dT_half*E_eff[2];
	
	//p_plus
	*pGamma = GAPS_APT_CalGamma(pP,pMass);
	if(pInputs->EMField_Cal_B)
	{
		p_minus2p_plus_exp(dT, B,*pGamma, pP,pCharge[0]/pMass[0]);
	}

	//p_k+1
	pP[0]+=dT_half*E_eff[0];
	pP[1]+=dT_half*E_eff[1];
	pP[2]+=dT_half*E_eff[2];

	//x_k+1
	*pGamma = GAPS_APT_CalGamma(pP,pMass);
	double Inv_mass_gamma_dT=dT/((*pGamma)*pMass[0]);
	pX[0]+=Inv_mass_gamma_dT*pP[0];
	pX[1]+=Inv_mass_gamma_dT*pP[1];
	pX[2]+=Inv_mass_gamma_dT*pP[2];
	
	*pT +=dT;
	// End: Core of algorithm
	return 0;
}

inline int p_minus2p_plus_exp(double dT,double *B,double gamma, double *Pp,double qOverM)
{

	double b1,b2,b3,px,py,pz,Bnorm;
	double gamma_rev = 1/gamma;
	Bnorm = sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]);	
	double Bnorm_rev=qOverM/Bnorm;

	b1    = B[0]*Bnorm_rev;
	b2    = B[1]*Bnorm_rev;
	b3    = B[2]*Bnorm_rev;

	double C1=sin(dT*Bnorm*gamma_rev);
	double C2=1.-cos(dT*Bnorm*gamma_rev);

	px    = Pp[0];
	py    = Pp[1];
	pz    = Pp[2];
	
	Pp[0]=-((-1 + pow(b2,2)*C2 + pow(b3,2)*C2)*px) + b3*C1*py + b1*b2*C2*py - b2*C1*pz + b1*b3*C2*pz;
	Pp[1]=-(b3*C1*px) + b1*b2*C2*px + py - pow(b1,2)*C2*py - pow(b3,2)*C2*py + b1*C1*pz + b2*b3*C2*pz;
	Pp[2]=b2*C1*px + b1*b3*C2*px - b1*C1*py + b2*b3*C2*py + pz - pow(b1,2)*C2*pz - pow(b2,2)*C2*pz;

	return 0;
}



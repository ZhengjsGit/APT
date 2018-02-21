#include "APT_AllHeaders.h"
inline int p_minus2p_plus(double dT,double *pB, double Lgamma,double *pP,double qOverM);

int GAPS_APT_Pusher_RVPA_Cay3D (Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
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
		p_minus2p_plus(dT, B,*pGamma, pP,pCharge[0]/pMass[0]);
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

inline int p_minus2p_plus(double dT,double *pB, double Lgamma,double *pP,double qOverM)
{
	double Bx,By,Bz,px,py,pz,temp;
	Bx    = qOverM*pB[0];
	By    = qOverM*pB[1];
	Bz    = qOverM*pB[2];
	px    = pP[0];
	py    = pP[1];
	pz    = pP[2];
	double dT_2=dT*dT;
	
	temp  = pow(Bx,2)*dT_2 + pow(By,2)*dT_2 + pow(Bz,2)*dT_2 + 4*pow(Lgamma,2);
	pP[0] = (pow(Bx,2)*dT_2*px - pow(By,2)*dT_2*px - pow(Bz,2)*dT_2*px + 2*Bx*By*dT_2*py +  2*Bx*Bz*dT_2*pz + 4*Bz*dT*py*Lgamma - 4*By*dT*pz*Lgamma + 4*px*pow(Lgamma,2))/temp;

	pP[1]=(2*Bx*By*dT_2*px - pow(Bx,2)*dT_2*py + pow(By,2)*dT_2*py - pow(Bz,2)*dT_2*py +  2*By*Bz*dT_2*pz - 4*Bz*dT*px*Lgamma + 4*Bx*dT*pz*Lgamma + 4*py*pow(Lgamma,2))/temp;

	pP[2]=(2*Bx*Bz*dT_2*px + 2*By*Bz*dT_2*py - pow(Bx,2)*dT_2*pz - pow(By,2)*dT_2*pz +   pow(Bz,2)*dT_2*pz + 4*By*dT*px*Lgamma - 4*Bx*dT*py*Lgamma + 4*pz*pow(Lgamma,2))/temp;
	return 0;
}

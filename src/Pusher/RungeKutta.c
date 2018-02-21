#include "APT_AllHeaders.h"
int GAPS_APT_LorentzFlow(double *pFlow,Gaps_APT_Particle *pPtc,long Dim, Gaps_IO_InputsContainer *pInputs);

int GAPS_APT_Pusher_RungeKutta(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double dT=pInputs->dT;
	double *ptau= GAPS_APT_GetS1(pPtc);
	double *pX4 = GAPS_APT_GetX4(pPtc);
	double *pP4 = GAPS_APT_GetP4(pPtc);
	double *pMass = GAPS_APT_GetMass1(pPtc);
	long Dim=pInputs->Pusher_RungeKutta_Dim;
	long Order=pInputs->Pusher_RungeKutta_Order;
	double k1[8],k2[8],k3[8],k4[8];
	if(Dim != 3 && Dim!=4)
	{
		assert(Dim==3 || Dim== 4);
		fprintf(stderr,"ERROR: in funciton GAPS_APT_Pusher_RungeKutta: Wrong value of Pusher_RungeKutta_Dim: must be 3 or 4\n");
	}
	
	double XP8_init[8]={pX4[0],pX4[1],pX4[2],pX4[3],pP4[0],pP4[1],pP4[2],pP4[3]};
	int i;
	if (2==Order)
	{
		GAPS_APT_LorentzFlow(k1,pPtc,Dim,pInputs);
		double half_dT=0.5*dT;
		for(i=0;i<8;i++)
		{
			pX4[i]=XP8_init[i]+half_dT*k1[i];
		}
		GAPS_APT_LorentzFlow(k2,pPtc,Dim,pInputs);
		for(i=0;i<8;i++)
		{
			pX4[i]=XP8_init[i]+dT*k2[i];
		}
		if(Dim==3)
		{
			pX4[4] = GAPS_APT_CalGamma(&pX4[5],pMass);	// Update \gamma
			pX4[0]=XP8_init[0]+dT;
		}
		else
		{
			ptau[0]+=dT;
		}
	}
	else if (3==Order)
	{
		GAPS_APT_LorentzFlow(k1,pPtc,Dim,pInputs);
		double dTOver3=dT/3;
		double dT2Over3=2*dT/3;
		for(i=0;i<8;i++)
		{
			pX4[i]=XP8_init[i]+dTOver3*k1[i];
		}
		GAPS_APT_LorentzFlow(k2,pPtc,Dim,pInputs);
		for(i=0;i<8;i++)
		{
			pX4[i]=XP8_init[i]+dT2Over3*k2[i];
		}
		GAPS_APT_LorentzFlow(k3,pPtc,Dim,pInputs);
		for(i=0;i<8;i++)
		{
			pX4[i]=XP8_init[i]+0.25*dT*(k1[i]+3*k3[i]);
		}
		if(Dim==3)
		{
			pX4[4] = GAPS_APT_CalGamma(&pX4[5],pMass);	// Update \gamma
			pX4[0]=XP8_init[0]+dT;
		}
		else
		{
			ptau[0]+=dT;
		}
	}
	else if (4==Order)
	{
		double dTOver2=0.5*dT;
		double dTOver6=dT/6;
		GAPS_APT_LorentzFlow(k1,pPtc,Dim,pInputs);
		for(i=0;i<8;i++)
		{
			pX4[i]=XP8_init[i]+dTOver2*k1[i];
		}
		GAPS_APT_LorentzFlow(k2,pPtc,Dim,pInputs);
		for(i=0;i<8;i++)
		{
			pX4[i]=XP8_init[i]+dTOver2*k2[i];
		}
		GAPS_APT_LorentzFlow(k3,pPtc,Dim,pInputs);
		for(i=0;i<8;i++)
		{
			pX4[i]=XP8_init[i]+dT*k3[i];
		}
		GAPS_APT_LorentzFlow(k4,pPtc,Dim,pInputs);
		for(i=0;i<8;i++)
		{
			pX4[i]=XP8_init[i]+dTOver6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
		}
		if(Dim==3)
		{
			pX4[4] = GAPS_APT_CalGamma(&pX4[5],pMass);	// Update \gamma
			pX4[0]=XP8_init[0]+dT;
		}
		else
		{
			ptau[0]+=dT;
		}
	}
	else
	{
		assert(Order==2 || Order == 3 || Order == 4);
		fprintf(stderr,"ERROR: in funciton GAPS_APT_Pusher_RungeKutta: Wrong value of Pusher_RungeKutta_Order: must be 2, 3 or 4\n");
	}
	return 0;
}

int GAPS_APT_LorentzFlow(double *pFlow,Gaps_APT_Particle *pPtc,long Dim, Gaps_IO_InputsContainer *pInputs)
{
	double *pX4 = GAPS_APT_GetX4(pPtc);
	double *pP4 = GAPS_APT_GetP4(pPtc);
	double *E = GAPS_APT_GetE3(pPtc);
	double *B = GAPS_APT_GetB3(pPtc);
	double *pMass = GAPS_APT_GetMass1(pPtc);
	double *pCharge= GAPS_APT_GetCharge1(pPtc);

	double F_ext[3],E_eff[3];
	GAPS_APT_CalEB(E,B,pPtc,pInputs);
	GAPS_APT_MergeExtForce(F_ext,pPtc,pInputs);
	E_eff[0]=pCharge[0]*E[0]+F_ext[0];
	E_eff[1]=pCharge[0]*E[1]+F_ext[1];
	E_eff[2]=pCharge[0]*E[2]+F_ext[2];

	double gamma=GAPS_APT_CalGamma(pP4+1,pMass);
	double mass_gamma_rev=1./(pMass[0]*gamma);
	if (Dim==3)
	{
		pFlow[0] = 1;//dt/dt
		pFlow[1] = pP4[1]*mass_gamma_rev;
		pFlow[2] = pP4[2]*mass_gamma_rev;
		pFlow[3] = pP4[3]*mass_gamma_rev;

		pFlow[4] = 0.;//don't need for 3D case

		double PCrossB[3];
		V3cross(pP4+1,B,PCrossB);
		pFlow[5] = E_eff[0] + pCharge[0]*mass_gamma_rev*PCrossB[0];
		pFlow[6] = E_eff[1] + pCharge[0]*mass_gamma_rev*PCrossB[1];
		pFlow[7] = E_eff[2] + pCharge[0]*mass_gamma_rev*PCrossB[2];
	}
	else if (Dim==4)
	{
		double mass_rev=1./pMass[0];
		double chargeOverMass=pCharge[0]*mass_rev;
		pFlow[0] = gamma;//
		pFlow[1] = pP4[1]*mass_rev;
		pFlow[2] = pP4[2]*mass_rev;
		pFlow[3] = pP4[3]*mass_rev;

		pFlow[4] = (E_eff[0]*pP4[1]+E_eff[2]*pP4[2]+E_eff[3]*pP4[3])*mass_rev;

		double PCrossB[3];
		V3cross(pP4+1,B,PCrossB);
		pFlow[5] = gamma*E_eff[0] + chargeOverMass*PCrossB[0];
		pFlow[6] = gamma*E_eff[1] + chargeOverMass*PCrossB[1];
		pFlow[7] = gamma*E_eff[2] + chargeOverMass*PCrossB[2];
	}
	return 0;
}

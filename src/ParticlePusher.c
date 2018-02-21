/*System libraries*/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>

/*GAPS-IO libraries*/
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
#include "gapsio.h"
#include "IO_Tools.h"

/*Math and constants libraries*/
#include "WYLfunc.h"
#include "Constants.h"

/*GAPS-APT libraries*/
#include "APT_PubFunc.h"
#include "EM_Field.h"
#include "ParticleStruct.h"
#include "ParticlePusher.h"
#include "External_Forces.h"

extern Gaps_APT_Particle **Global_pCurrentPtc;

inline int GAPS_APT_Push_ParticleGroup(Gaps_APT_ParticleGroup *pPtcgrp,Gaps_APT_ParticlePusher pPusher,Gaps_IO_InputsContainer *pInputs)
{
	long idx_ptc;
	for(idx_ptc=0;idx_ptc<pPtcgrp->count_ptc;idx_ptc++)	
	{
		Gaps_APT_Particle *pPtc=pPtcgrp->ptcgrp+idx_ptc;
		double *die=GAPS_APT_GetDie1(pPtc);
		int Dead=(int)(*die);
		if(!Dead)
		{
			*Global_pCurrentPtc = pPtc;

			//Save phase space cache for calculation of work and acceleration
			GAPS_APT_SavePhaseSpace2Cache(pPtc,pInputs);
			// Push Particle
			(*pPusher)(pPtc,pInputs);

			//Calculate Work for all loaded forces
			GAPS_APT_CalForceWorks(pPtc,pInputs);

			//Calculate acceleration for Radiation
			GAPS_APT_CalAcceleration(pPtc,pInputs);
		}
	}
	return 0;
}

inline int GAPS_APT_SavePhaseSpace2Cache(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	if(pInputs->Open_Cal_Work || pInputs->Open_Cal_Acceleration)
	{
		double *pT = GAPS_APT_GetT1(pPtc);
		double *pX = GAPS_APT_GetX3(pPtc);
		double *pGamma=GAPS_APT_GetGamma1(pPtc);
		double *pP = GAPS_APT_GetP3(pPtc);
		pPtc->cache[0] = pT[0];
		pPtc->cache[1] = pX[0];
		pPtc->cache[2] = pX[1];
		pPtc->cache[3] = pX[2];
		pPtc->cache[4] = pGamma[0];
		pPtc->cache[5] = pP[0];
		pPtc->cache[6] = pP[1];
		pPtc->cache[7] = pP[2];
	}
	return 0;
}

inline int GAPS_APT_CalAcceleration(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	if(pInputs->Open_Cal_Acceleration)
	{
		double *pAclr = GAPS_APT_GetAclr3(pPtc);
		double *pGamma=GAPS_APT_GetGamma1(pPtc);
		double *pP = GAPS_APT_GetP3(pPtc);
		double *P0=pPtc->cache+5;
		double gamma0=pPtc->cache[4];

		double dT_inv = 1/(pInputs->dT);
		double gamma0_inv=1/gamma0;
		double gamma1_inv=1/(*pGamma);

		pAclr[0]=(pP[0]*gamma1_inv-P0[0]*gamma0_inv)*dT_inv;
		pAclr[1]=(pP[1]*gamma1_inv-P0[1]*gamma0_inv)*dT_inv;
		pAclr[2]=(pP[2]*gamma1_inv-P0[2]*gamma0_inv)*dT_inv;
		//	printf("Aclr=%e,%e,%e\n",pAclr[0],pAclr[1],pAclr[2]);
	}
	return 0;
}

inline int GAPS_APT_CalEB(double *E,double *B,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *ST=GAPS_APT_GetX4(pPtc);
	double EB[6];
	GAPS_APT_CalField(EB,pPtc,ST,-1,pInputs);
	E[0]=EB[0];
	E[1]=EB[1];
	E[2]=EB[2];
	B[0]=EB[3];
	B[1]=EB[4];
	B[2]=EB[5];
	return 0;
}

inline int GAPS_APT_CalA(double *A,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *ST=GAPS_APT_GetX4(pPtc);
	GAPS_APT_CalField(A,pPtc,ST,1,pInputs);
	return 0;
}


inline int GAPS_APT_UpdatePtcData_EB(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *ST=GAPS_APT_GetX4(pPtc);
	double EB[6];
	double *E=GAPS_APT_GetE3(pPtc);
	double *B=GAPS_APT_GetB3(pPtc);
	GAPS_APT_CalField(EB,pPtc,ST,-1,pInputs);
	E[0]=EB[0];
	E[1]=EB[1];
	E[2]=EB[2];
	B[0]=EB[3];
	B[1]=EB[4];
	B[2]=EB[5];
	return 0;
}

inline int GAPS_APT_UpdatePtcData_A4(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *ST=GAPS_APT_GetX4(pPtc);
	double *pA4=GAPS_APT_GetA4(pPtc);
	GAPS_APT_CalField(pA4,pPtc,ST,1,pInputs);
	return 0;
}


//Generated functions
#include "Pusher_Set.h"

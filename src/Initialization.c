/*System libraries*/
#include "APT_AllHeaders.h"

int GAPS_APT_ParticleInitialCondition(Gaps_APT_ParticleGroup *pPtcgrp,Gaps_IO_InputsContainer *pInputs)
{
	long i;
	long num_ptc_init=(pInputs->Init_Num_Particles)/(pPtcgrp->proc_size);
	//could be zero danger
	for(i=0;i<num_ptc_init;i++)
	{
		Gaps_APT_Particle *pPtc=pPtcgrp->ptcgrp + i;
		GAPS_APT_SetParticleStatus(pPtc,pInputs);
		GAPS_APT_SetParticleType(pPtc,pInputs);
		GAPS_APT_SetParticlePosition(pPtc,pInputs);
		GAPS_APT_SetParticleMomentum(pPtc,pInputs);
		GAPS_APT_SetParticleAcceleration(pPtc,pInputs);
		(pPtcgrp->count_ptc)++;
	}

	MPI_Barrier(pPtcgrp->mpi_comm);
	//Initialize EM field and Force
	for(i=0;i<num_ptc_init;i++)
	{
		Gaps_APT_Particle *pPtc=pPtcgrp->ptcgrp + i;
		GAPS_APT_UpdatePtcData_EB(pPtc,pInputs);	
		GAPS_APT_UpdatePtcData_A4(pPtc,pInputs);	
		GAPS_APT_UpdatePtcData_ExtForce(pPtc,pInputs);	
	}
	return 0;
}

int GAPS_APT_SetParticleStatus(Gaps_APT_Particle *pPtc, Gaps_IO_InputsContainer *pInputs)
{
	if(!(strcmp(pInputs->Init_Status_Type,"Constant")))
	{
		GAPS_APT_SetParticleStatus_Constant(pPtc,pInputs);
	}
	else
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_SetParticleType: Initial particle status is wrong! You input Init_Ptc_Status = %s\n",pInputs->Init_Status_Type);
	}
	return 0;
}

int GAPS_APT_SetParticleType(Gaps_APT_Particle *pPtc, Gaps_IO_InputsContainer *pInputs)
{
	if(!(strcmp(pInputs->Init_Ptc_Type,"Constant")))
	{
		GAPS_APT_SetParticleType_Constant(pPtc,pInputs);
	}
	else
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_SetParticleType: Initial particle type is wrong! You input Init_Aclr_Type = %s\n",pInputs->Init_Ptc_Type);
	}
	return 0;
}

int GAPS_APT_SetParticleStatus_Constant(Gaps_APT_Particle *pPtc, Gaps_IO_InputsContainer *pInputs)
{
	double *die=GAPS_APT_GetDie1(pPtc);
	*die = pInputs->Init_Status_Constant_IsDead;
	return 0;
}

int GAPS_APT_SetParticleType_Constant(Gaps_APT_Particle *pPtc, Gaps_IO_InputsContainer *pInputs)
{
	double *charge=GAPS_APT_GetCharge1(pPtc);
	double *mass=GAPS_APT_GetMass1(pPtc);

	*charge = pInputs->Init_Ptc_Constant_ChargeMass[0];
	*mass= pInputs->Init_Ptc_Constant_ChargeMass[1];
	return 0;
}
//Generated functions
#include "Init_Set.h"

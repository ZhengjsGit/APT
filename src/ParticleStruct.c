#include "APT_AllHeaders.h"

int GAPS_APT_InitParticleGroup(Gaps_APT_ParticleGroup *pPtcgrp,Gaps_IO_InputsContainer *pInputs,int Rank,int Mpi_size,MPI_Datatype Datatype,MPI_Comm Comm)
{
	pPtcgrp->proc_rank = Rank;
	pPtcgrp->proc_size = Mpi_size;
	pPtcgrp->mpi_datatype= Datatype;
	pPtcgrp->mpi_comm= Comm;
	
	long total_num=pInputs->num_total_particles;
	pPtcgrp->num_total_particles = total_num-total_num%Mpi_size;
	pPtcgrp->num_local_ptc=pInputs->num_total_particles/Mpi_size;
//	printf("num_total=%d,num_local=%d\n",pPtcgrp->num_total_particles,pPtcgrp->num_local_ptc);

	pPtcgrp->ptcgrp = (Gaps_APT_Particle *)calloc(pPtcgrp->num_local_ptc,sizeof(Gaps_APT_Particle));
	pPtcgrp->count_ptc= 0;
	long i;
	for(i=0;i<pPtcgrp->num_local_ptc;i++)
	{
		pPtcgrp->ptcgrp[i].Pusher_Type = pInputs->Pusher_Type;
		pPtcgrp->ptcgrp[i].cache=calloc(GAPS_APT_CONST_PTC_PHASESPACE_CACHE_SIZE,sizeof(double));
		GAPS_APT_SetExtForceInfo(pPtcgrp->ptcgrp+i,pInputs);
		GAPS_APT_SetFieldFunction(pPtcgrp->ptcgrp+i,pInputs);
	}
	int ptc_data_total_len=GAPS_APT_GetTotalPtcDataLen(pPtcgrp);
	pPtcgrp->ptc_data_total_len=ptc_data_total_len;
	for(i=0;i<pPtcgrp->num_local_ptc;i++)
	{
		pPtcgrp->ptcgrp[i].data=calloc(ptc_data_total_len,sizeof(double));
	}
	return 0;
}

int GAPS_APT_EraseParticleGroup(Gaps_APT_ParticleGroup *pPtcgrp)
{
	long i;
	for(i=0;i<pPtcgrp->num_local_ptc;i++)
	{
		free(pPtcgrp->ptcgrp[i].cache);
		free(pPtcgrp->ptcgrp[i].data);
		free(pPtcgrp->ptcgrp[i].ForceType);
		free(pPtcgrp->ptcgrp[i].ForceName);
	}
	free(pPtcgrp->ptcgrp);
	return 0;
}

inline int GAPS_APT_GetTotalPtcDataLen(Gaps_APT_ParticleGroup *pPtcgrp)
{
	return GAPS_APT_CONST_PTC_DATA_LEN+((pPtcgrp->ptcgrp)->num_forces) * (GAPS_APT_CONST_EXTERN_FORCE_DIM+1);
}

inline double *GAPS_APT_GetForceData_FromData(double *data,int force_idx)
{
	return data+GAPS_APT_CONST_PTC_DATA_LEN+force_idx*GAPS_APT_CONST_EXTERN_FORCE_DIM;
}

inline double *GAPS_APT_GetForceData(Gaps_APT_Particle *pPtc,int force_idx)
{
	return pPtc->data+GAPS_APT_CONST_PTC_DATA_LEN+force_idx*GAPS_APT_CONST_EXTERN_FORCE_DIM;
}

inline double *GAPS_APT_GetForceWorkData_FromData(double *data,int force_idx,int num_forces)
{
	return data+GAPS_APT_CONST_PTC_DATA_LEN+num_forces*GAPS_APT_CONST_EXTERN_FORCE_DIM+force_idx;
}

inline double *GAPS_APT_GetForceWorkData(Gaps_APT_Particle *pPtc,int force_idx)
{
	return pPtc->data+GAPS_APT_CONST_PTC_DATA_LEN+(pPtc->num_forces)*GAPS_APT_CONST_EXTERN_FORCE_DIM + force_idx;
}

double *GAPS_APT_GetForcePointer(Gaps_APT_Particle *pPtc,int Force_Type)
{
	if(0==(pPtc->num_forces))
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_GetForcePointer: no force is loaded!\n");		
		return NULL;
	}
	else
	{
		int index=0,i;	
		for(i=0;i<pPtc->num_forces;i++)
		{
			if(pPtc->ForceType[i] == Force_Type)
			{
				index=i;
			}
		}
		return GAPS_APT_GetForceData(pPtc,index);
	}
}

inline int GAPS_APT_CalField(double *pTensor,Gaps_APT_Particle *pPtc,double *pSpaceTime4,int Order, Gaps_IO_InputsContainer *pInputs)
{
	GAPS_APT_ResetTensor(pTensor,Order);
	(pPtc->FieldFunc)(pTensor,pSpaceTime4,Order,pInputs);
	return 0;
}

inline double GAPS_APT_CalField_Element(double *pTensor,Gaps_APT_Particle *pPtc,double *pSpaceTime4,int Order,int *pIndex, Gaps_IO_InputsContainer *pInputs)//Need optimization
{
	GAPS_APT_ResetTensor(pTensor,Order);
	(pPtc->FieldFunc)(pTensor,pSpaceTime4,Order,pInputs);
	return GAPS_APT_TensorValue(pTensor,pIndex,4,Order);
}

double GAPS_APT_IntegralField(int n,Gaps_APT_Particle *pPtc,int Order,int *pIndex,int idxInt,double *pST0,double *pST1,Gaps_IO_InputsContainer *pInputs)
{
	int len;
	if(Order>=1)
	{
		len=GAPS_APT_PowerInt(4,Order);
	}
	else if(Order==-1)
	{
		len=6;
	}
	else
	{
		assert(Order>=-1);
		fprintf(stderr,"ERROR: In function GAPS_APT_IntegralField: Order must be integers: -1, 1, 2,... You input: Order = %d\n ",Order);
	}
	double *Tensor=(double *)calloc(len,sizeof(double));
	
	double a=pST0[idxInt];
	double b=pST1[idxInt];
	double h=(b-a)/n;

	//printf("len=%d,n=%d,pIndex=%d,psT0={%e,%e,%e,%e},pST1={%e,%e,%e,%e}\n",len,n,pIndex[0],pST0[0],pST0[1],pST0[2],pST0[3],pST1[0],pST1[1],pST1[2],pST1[2]);
	double f_a = GAPS_APT_CalField_Element(Tensor,pPtc,pST0,Order,pIndex,pInputs);
//	printf("f_a=%e\n",f_a);
	double f_b = GAPS_APT_CalField_Element(Tensor,pPtc,pST1,Order,pIndex,pInputs);
	double XI0=f_a+f_b;

	double XI1=0.;
	double XI2=0.;
	int i;
	for(i=1;i<n;i+=2)
	{
        pST0[idxInt] = a+i*h;
        XI1 += GAPS_APT_CalField_Element(Tensor,pPtc,pST0,Order,pIndex,pInputs);
	}
	for(i=2;i<n-1;i+=2)
	{
        pST0[idxInt] = a+i*h;
        XI2 += GAPS_APT_CalField_Element(Tensor,pPtc,pST0,Order,pIndex,pInputs);
	}
	free(Tensor);
	pST0[idxInt]=a;//Reset value of pST0
	return h*(XI0+2*XI1+4*XI2)/3;
}

//Generated functions
#include "EMField_Set.h"
#include "ExtForce_Set.h"

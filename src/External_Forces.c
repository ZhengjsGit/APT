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
#include "External_Forces.h"

//Generated
#include "ExtForce_FuncContainer.h"

inline int GAPS_APT_CalExtForce(double *pForce3,int Type,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	int count_force=pPtc->num_forces;
	if(0!=count_force)
	{
		(*(Global_ExtForce_Container[Type]))(pForce3,pPtc,pInputs);
	}
	else
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_CalExtForce: No external force is loaded!\n");
	}
	return 0;
}

inline int GAPS_APT_UpdatePtcData_ExtForce(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	int count_force=pPtc->num_forces;
	if(0!=count_force)
	{
		int i;
		int *Types=pPtc->ForceType;
		double *ForceData=GAPS_APT_GetForceData(pPtc,0);
		for(i=0;i<count_force;i++)
		{
			int DataOffset=i*GAPS_APT_CONST_EXTERN_FORCE_DIM;
			(*(Global_ExtForce_Container[Types[i]]))(ForceData+DataOffset,pPtc,pInputs);
		}
	}
	return 0;
}

inline int GAPS_APT_MergeExtForce(double *F_ext3,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	int count_force=pPtc->num_forces;
	if(0!=count_force)
	{
		F_ext3[0]=0.;	
		F_ext3[1]=0.;	
		F_ext3[2]=0.;	

		int i;
		int *Types=pPtc->ForceType;
		double *ForceData=GAPS_APT_GetForceData(pPtc,0);
		for(i=0;i<count_force;i++)
		{
			int DataOffset=i*GAPS_APT_CONST_EXTERN_FORCE_DIM;
			(*(Global_ExtForce_Container[Types[i]]))(ForceData+DataOffset,pPtc,pInputs);
			F_ext3[0]+=ForceData[DataOffset+0];	
			F_ext3[1]+=ForceData[DataOffset+1];	
			F_ext3[2]+=ForceData[DataOffset+2];	
		}
	}
	else
	{
		F_ext3[0]=0.;	
		F_ext3[1]=0.;	
		F_ext3[2]=0.;	
	}
	return 0;
}

inline int GAPS_APT_CalForceWorks(Gaps_APT_Particle *pPtc, Gaps_IO_InputsContainer *pInputs)
{
	if(pInputs->Open_Cal_Work)
	{
		double *pX = GAPS_APT_GetX3(pPtc);
		double *X0=pPtc->cache+1;
		double dX[3]={pX[0]-X0[0],pX[1]-X0[1],pX[2]-X0[2]};
		double gamma0=pPtc->cache[4];

		if(pInputs->EMField_Cal_E)
		{
			//Calculate work of electric field
			double *pE=GAPS_APT_GetE3(pPtc);
			double *pEwork=GAPS_APT_GetEwork1(pPtc);
			*pEwork +=pE[0]*dX[0]+pE[1]*dX[1]+pE[2]*dX[2];
		//printf("Ework=%e\n",pEwork[0]);
		}

		int count_force=pPtc->num_forces;
		if(0!=count_force)
		{
			int i;
			double *ForceData=GAPS_APT_GetForceData(pPtc,0);
			double *WorkData=GAPS_APT_GetForceWorkData(pPtc,0);
			for(i=0;i<count_force;i++)
			{
				int DataOffset=i*GAPS_APT_CONST_EXTERN_FORCE_DIM;
				double *force=ForceData+DataOffset;
				WorkData[i] += force[0]*dX[0]+force[1]*dX[1]+force[2]*dX[2];
			}
		}
	}
	return 0;
}

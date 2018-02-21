#include "ParticleStructMacro.h"
#define GAPS_APT_CONST_EXTERN_FORCE_DIM 3
#define GAPS_APT_CONST_PTC_PHASESPACE_CACHE_SIZE 8

/*! Double direction chain list for particles */
typedef struct {
		double *data;/*! Array of particle data */		
		int Pusher_Type;
		double *cache;
		Gaps_APT_Field FieldFunc;

		//Extern forces
		int num_forces;
		char **ForceName;
		int *ForceType;
} Gaps_APT_Particle;

typedef struct {
		//Global parameters
		long num_total_particles;//total number of particle
		int ptc_data_total_len;
		
		//Local parameters
		long num_local_ptc;
		long count_ptc;

		/*MPI parameters*/
		int proc_rank;
		int proc_size;
		MPI_Datatype mpi_datatype;
		MPI_Comm mpi_comm;

		Gaps_APT_Particle *ptcgrp;
} Gaps_APT_ParticleGroup;

int GAPS_APT_InitParticleGroup(Gaps_APT_ParticleGroup *pPtcgrp,Gaps_IO_InputsContainer *pInputs,int Rank,int Mpi_size,MPI_Datatype Datatype,MPI_Comm Comm);
int GAPS_APT_EraseParticleGroup(Gaps_APT_ParticleGroup *pPtcgrp);
int GAPS_APT_SetFieldFunction(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
int GAPS_APT_SetExtForceInfo(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
inline int GAPS_APT_GetTotalPtcDataLen(Gaps_APT_ParticleGroup *pPtcgrp);
inline double *GAPS_APT_GetForceData_FromData(double *data,int force_idx);
inline double *GAPS_APT_GetForceWorkData_FromData(double *data,int force_idx,int num_forces);
inline double *GAPS_APT_GetForceData(Gaps_APT_Particle *pPtc,int force_idx);
inline double *GAPS_APT_GetForceWorkData(Gaps_APT_Particle *pPtc,int force_idx);

double *GAPS_APT_GetForcePointer(Gaps_APT_Particle *pPtc,int Force_Type);

inline int GAPS_APT_CalField(double *pTensor,Gaps_APT_Particle *pPtc,double *pSpaceTime4,int Order, Gaps_IO_InputsContainer *pInputs);
inline double GAPS_APT_CalField_Element(double *pTensor,Gaps_APT_Particle *pPtc,double *pSpaceTime4,int Order,int *pIndex, Gaps_IO_InputsContainer *pInputs);
double GAPS_APT_IntegralField(int n,Gaps_APT_Particle *pPtc,int Order,int *pIndex,int idxInt,double *pST0,double *pST1,Gaps_IO_InputsContainer *pInputs);

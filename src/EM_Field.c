#include "APT_AllHeaders.h"

#define GAPS_APT_CONST_EB_TENSOR_DIM 6

extern Gaps_APT_DiscreteTensor *Global_pDisTensor;
extern Gaps_APT_Particle **Global_pCurrentPtc;
extern Gaps_APT_ParticleGroup *Global_pPtcgrp;

inline int GAPS_APT_PowerInt(int x,int n)
{
	int y=1;
	int i;
	if(n<=0)
	{
		y=1;
	}
	else
	{
		for(i=0;i<n;i++)
		{
			y*=x;
		}
	}
	return y;
}

inline int GAPS_APT_ResetTensor(double *pTensor,int Order)
{
	int i=0;
	if(-1 == Order)
	{
		for(i=0;i<GAPS_APT_CONST_EB_TENSOR_DIM;i++)
		{
			pTensor[i]=0.;	
		}
	}
	else if(Order > 0)
	{
		int len=GAPS_APT_PowerInt(4,Order);
		for(i=0;i<len;i++)
		{
			pTensor[i]=0.;	
		}
	}
	else
	{
		fprintf(stderr,"ERROR: In GAPS_APT_ResetTensor: Order of tensors in APT cannot be 0 or less than -1");	
	}
	return 0;
}

inline int GAPS_APT_TensorLen(int Order)
{
	if(Order == -1)
	{
		return 6;
	}
	else if(Order > 0)
	{
		return GAPS_APT_PowerInt(4,Order);
	}
	else
	{
		return 0;
	}
}

double * GAPS_APT_GetFieldVector(double *data,long *pWT_idx,Gaps_APT_DiscreteTensor *pDisTensor)
{
	long i=pWT_idx[0];	
	long j=pWT_idx[1];	
	long k=pWT_idx[2];	
	long nx=pDisTensor->xnum;
	long ny=pDisTensor->ynum;
	long dim=pDisTensor->dim;

	return data+dim*(i+nx*(j+ny*k));
}

//Need optimization
int GAPS_APT_Field_Discrete(double *pTensor,double *pSpaceTime4,int Order, Gaps_IO_InputsContainer *pInputs)
{
//Global variables: Global_pDisTensor, Global_pCurrentPtc
	long Idx[4],Rel_Idx[4];
	int InBoundary;
	InBoundary=GAPS_APT_Field_Discrete_GetGridIdx(Idx,Rel_Idx,Global_pDisTensor,pSpaceTime4);
	//printf("Idx=%ld,%ld,%ld,%ld, IN?=%d\n",Idx[0],Idx[1],Idx[2],Idx[3],InBoundary);
	if(InBoundary)
	{

		//Weighting grid: (i,j), (i,j+1), (i+1,j), (i+1,j+1)
		long lenT = GAPS_APT_TensorLen(Order);
		double *Tij=(double *)calloc(lenT,sizeof(double));
		int i,j,k,l;
		for(i=0;i<2;i++)
		{
			for(j=0;j<2;j++)	
			{
				for(k=0;k<2;k++)	
				{	
					long WT_idx[3]={Idx[1]+i,Idx[2]+j,Idx[3]+k};
					double *Field=GAPS_APT_GetFieldVector(Global_pDisTensor->data[Idx[0]],WT_idx,Global_pDisTensor);
					for(l=0;l<lenT;l++)
					{
						Tij[l]+=Field[l];
					}
				}
			}
		}
		for(l=0;l<lenT;l++)
		{
			pTensor[l]=Tij[l]*0.125;
		}
	}
	else
	{
	//printf("Fuck1\n");
		double *die=GAPS_APT_GetDie1(*Global_pCurrentPtc);	
		*die=1;
	}
	return 0;
}

inline int GAPS_APT_Field_Discrete_GetGridIdx(long *pIdx4,long *pRelIdx4,Gaps_APT_DiscreteTensor *pDisTensor,double *pST4)
{
	int InX,InY,InZ,InT;
	int StaticField;

	double tmin=pDisTensor->boundary[0];
	double tmax=pDisTensor->boundary[1];
	if(tmin == tmax)
	{
		StaticField=1;	
	}
	else
	{
		StaticField=0;	
	}

	double xmin=pDisTensor->boundary[2];
	double xmax=pDisTensor->boundary[3];
	double ymin=pDisTensor->boundary[4];
	double ymax=pDisTensor->boundary[5];
	double zmin=pDisTensor->boundary[6];
	double zmax=pDisTensor->boundary[7];

	double DT=(pST4[0]-0.)/(pDisTensor->dt);
	double DX=(pST4[1]-xmin)/(pDisTensor->dx);
	double DY=(pST4[2]-ymin)/(pDisTensor->dy);
	double DZ=(pST4[3]-zmin)/(pDisTensor->dz);

	//printf("t=%e,x=%e,y=%e,z=%e\n",pST4[0],pST4[1],pST4[2],pST4[3]);
	//printf("xmin=%e,xmax=%e,ymin=%e,ymax=%e,zmin=%e,zmax=%e,DT=%e,DX=%e,DY=%e,DZ=%e\n",xmin,xmax,ymin,ymax,zmin,zmax,DT,DX,DY,DZ);

	if(StaticField)
	{
		InT=1;	
		pIdx4[0]=0;
		pRelIdx4[0] = 0;
	}
	else
	{
		if(DT<0. || DT >= (pDisTensor->tnum))
		{
			InT=0;
		}
		else
		{
			InT=1;
			pIdx4[0]=(long)DT;	
			pRelIdx4[0] = DT-pIdx4[0];
		}
	}

	if(DX<0. || DX >= (pDisTensor->xnum))
	{
		InX=0;	
	}
	else
	{
		InX=1;
		pIdx4[1]=(long)DX;
		pRelIdx4[1] = DX-pIdx4[1];
	}

	if(DY<0. || DY >= (pDisTensor->ynum))
	{
		InY=0;	
	}
	else
	{
		InY=1;
		pIdx4[2]=(long)DY;
		pRelIdx4[2] = DY-pIdx4[2];
	}

	if(DZ<0. || DZ >= (pDisTensor->znum))
	{
		InZ=0;	
	}
	else
	{
		InZ=1;
		pIdx4[3]=(long)DZ;
		pRelIdx4[3] = DZ-pIdx4[3];
	}
	//printf("InT=%d,InX=%d,InY=%d,InZ=%d\n",InT,InX,InY,InZ);

	return (InT && InX && InY && InZ);
}

int GAPS_APT_Field_Discrete_LoadData(Gaps_APT_DiscreteTensor *pDisTensor,char *pFilename)
{
#ifdef DATA_FORMAT_HDF5
	hid_t file_id, dataset_id;
	herr_t status;

	file_id = H5Fopen(pFilename, H5F_ACC_RDONLY, H5P_DEFAULT);

	double N_grid[4];
	dataset_id = H5Dopen2(file_id, "/N_grid", H5P_DEFAULT);
	status = H5Dread(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_grid);
	status = H5Dclose(dataset_id);

	double Order;
	dataset_id = H5Dopen2(file_id, "/Order", H5P_DEFAULT);
	status = H5Dread(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Order);
	status = H5Dclose(dataset_id);

	double DX;
	dataset_id = H5Dopen2(file_id, "/DX", H5P_DEFAULT);
	status = H5Dread(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &DX);
	status = H5Dclose(dataset_id);

	double DT;
	dataset_id = H5Dopen2(file_id, "/DT", H5P_DEFAULT);
	status = H5Dread(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &DT);
	status = H5Dclose(dataset_id);

	double OriginPoint[3];
	dataset_id = H5Dopen2(file_id, "/OriginPoint", H5P_DEFAULT);
	status = H5Dread(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, OriginPoint);
	status = H5Dclose(dataset_id);

	long num_data=1;
	pDisTensor->xnum=(long)N_grid[0];
	num_data*=(pDisTensor->xnum);
	pDisTensor->ynum=(long)N_grid[1];
	num_data*=(pDisTensor->ynum);
	pDisTensor->znum=(long)N_grid[2];
	num_data*=(pDisTensor->znum);
	pDisTensor->tnum=(long)N_grid[3];
	num_data*=(pDisTensor->tnum);
	pDisTensor->order=(long)Order;

	if(-1 == pDisTensor->order)
	{
		pDisTensor->dim=6;
	}
	else if(0 < pDisTensor->order)
	{
		pDisTensor->dim = GAPS_APT_PowerInt(4,pDisTensor->order);
	}
	num_data*=(pDisTensor->dim);
	pDisTensor->dt=DT;
	pDisTensor->dx=DX;
	pDisTensor->dy=pDisTensor->dx;
	pDisTensor->dz=pDisTensor->dx;

	pDisTensor->boundary[0] = 0.;
	pDisTensor->boundary[1] = 0.+(pDisTensor->dt)*((pDisTensor->tnum)-1);
	pDisTensor->boundary[2] = OriginPoint[0];
	pDisTensor->boundary[3] = OriginPoint[0]+(pDisTensor->dx)*((pDisTensor->xnum)-1);
	pDisTensor->boundary[4] = OriginPoint[1];
	pDisTensor->boundary[5] = OriginPoint[1]+(pDisTensor->dy)*((pDisTensor->ynum)-1);
	pDisTensor->boundary[6] = OriginPoint[2];
	pDisTensor->boundary[7] = OriginPoint[2]+(pDisTensor->dz)*((pDisTensor->znum)-1);

	pDisTensor->data = (double **)calloc(pDisTensor->tnum,sizeof(double *));
	pDisTensor->data1D = (double *)calloc(num_data,sizeof(double));
	long num_per_time = num_data/(pDisTensor->tnum);
	long i;
	pDisTensor->num_data_per_moment = num_per_time;

	dataset_id = H5Dopen2(file_id, "/Data", H5P_DEFAULT);
	for(i=0;i<(pDisTensor->tnum);i++)
	{
		long offset = i*num_per_time;
		pDisTensor->data[i] = &(pDisTensor->data1D[offset]);	
	}
	double data[pDisTensor->tnum][num_per_time];
	status = H5Dread(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pDisTensor->data1D);
	//status = H5Dread(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	status = H5Dclose(dataset_id);
	status = H5Fclose(file_id);
#endif

	/*	printf("DisF:\n");
		printf("N_grid=%ld,%ld,%ld,%ld\n",(long)N_grid[0],(long)N_grid[1],(long)N_grid[2],(long)N_grid[3]);
		printf("Dim=%ld\n",pDisTensor->dim);
		printf("Dx=%f,Dt=%f\n",DX,DT);
		printf("OP=%f,%f,%f\n",OriginPoint[0],OriginPoint[1],OriginPoint[2]);
		long j;
		for(i=0;i<num_data;i++)
		{
		printf("%f ",pDisTensor->data1D[i]);	
		}
		printf("DisF:\n");
		for(i=0;i<num_per_time;i++)
		{
		printf("%f ",pDisTensor->data[0][i]);	
		}
		*/
#ifdef DATA_FORMAT_GAPSIO
#endif
	return 0;
}

int GAPS_APT_EraseDiscreteTensor(Gaps_APT_DiscreteTensor *pDisTensor)
{
	free(pDisTensor->data1D);
	free(pDisTensor->data);
	return 0;	
}

int GAPS_APT_GatherPtcInfo(double *pData)
{
	int num_local_ptc=Global_pPtcgrp->num_local_ptc;
	int len_local_data=6*num_local_ptc;
	double *local_data_cache = (double *)calloc(len_local_data,sizeof(double));
	
	//Copy data to cache
	long i;
	for(i=0;i<num_local_ptc;i++)
	{
		Gaps_APT_Particle *pPtc=Global_pPtcgrp->ptcgrp+i;
		double XP[6];
		double *X=GAPS_APT_GetX3(pPtc),*P=GAPS_APT_GetP3(pPtc);
		memcpy(XP+0,X,sizeof(double)*3);
		memcpy(XP+3,P,sizeof(double)*3);
		memcpy(local_data_cache+i*6,XP,sizeof(double)*6);
	}

	int num_total_ptc=Global_pPtcgrp->num_total_particles;
	int len_global_data=6*num_total_ptc;
	double *global_data_cache = (double *)calloc(len_global_data,sizeof(double));
	if(NULL == global_data_cache)
	{
		assert(NULL != global_data_cache);
		fprintf(stderr,"ERROR: During particle information gathering: fail to allocate RAM, particle number might be too large\n");
	}

//	printf("n_total=%d,n_local=%d\n",num_total_ptc,num_local_ptc);
	//Gather data by Rank 0
	int ToRank=0,tag=Global_pPtcgrp->proc_rank;
	if(0==tag)
	{
		memcpy(global_data_cache,local_data_cache,sizeof(double)*len_local_data);
		int rank_j;
		MPI_Status status;
		for(rank_j=1;rank_j<Global_pPtcgrp->proc_size;rank_j++)
		{
			MPI_Recv(global_data_cache+rank_j*len_local_data,len_local_data,Global_pPtcgrp->mpi_datatype,rank_j,rank_j,Global_pPtcgrp->mpi_comm,&status);
		}
	}
	else
	{
		MPI_Send(local_data_cache,len_local_data,Global_pPtcgrp->mpi_datatype,ToRank,tag,Global_pPtcgrp->mpi_comm);
	}
	//MPI_Barrier(Global_pPtcgrp->mpi_comm);

	//Rank 0 broadcast data to all processors
	MPI_Bcast(global_data_cache,len_global_data,Global_pPtcgrp->mpi_datatype,0,Global_pPtcgrp->mpi_comm);

	//Copy data to the output
	memcpy(pData,global_data_cache,len_global_data*sizeof(double));

	free(global_data_cache);
	free(local_data_cache);
	return 0;
}

/*System libraries*/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <sys/stat.h>

/*GAPS-IO libraries*/
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
#include "gapsio.h"
#include "IO_Tools.h"

/*HDF5 header*/
#ifdef DATA_FORMAT_HDF5
#include "hdf5.h"
#endif

/*Math and constants libraries*/
#include "WYLfunc.h"
#include "Constants.h"

/*GAPS-APT libraries*/
#include "APT_PubFunc.h"
#include "EM_Field.h"
#include "ParticleStruct.h"
#include "External_Forces.h"
#include "ParticlePusher.h"
#include "Initialization.h"
#include "APT_Output.h"

int GAPS_APT_GenOutputFileName(char *pFname,char *pDirIN,char *filename)
{
	//Generate output filename
	char pDirOut[GAPS_APT_CONST_DIRSTR_LEN];
	if(!strcmp(pDirIN,"c"))//Deflaut
	{
		strcpy(pDirOut,"./");
	}
	else
	{
		GAPS_APT_ProducePath(pDirOut,pDirIN);
	}

	strcat(pDirOut,filename);
	strcpy(pFname,pDirOut);
	return 0;
}

int GAPS_APT_ProducePath(char *pDirOut,char *pDirIN)
{
	char * position;
	int length,pp;

	if(pDirIN!=NULL)
	{
		position=strrchr(pDirIN, '/');

		if (position)
		{
			pp = (int)(position-pDirIN);
			length=strlen(pDirIN);

			if(pp+1<length)
			{
				strcpy(pDirOut, pDirIN);
				strcat(pDirOut,"/");
			}
			else
			{
				strcpy(pDirOut,pDirIN);
			}
		}
		else
		{
			strcpy(pDirOut,pDirIN);
			strcat(pDirOut,"/");
		}
	}
	else
	{
		char WorkingDir[GAPS_APT_CONST_DIRSTR_LEN];
		assert(NULL != getcwd(WorkingDir,sizeof(WorkingDir)));
		strcpy(pDirOut,WorkingDir);
		strcat(pDirOut,"/");
	}
	return 0;
}

int GAPS_APT_InitOutputFiles(Gaps_APT_OutputFiles *pOutFiles,Gaps_APT_ParticleGroup *pPtcgrp,char *DIR_Output,Gaps_IO_InputsContainer *pInputs)
{
	long len_data_per_ptc=GAPS_APT_CONST_PTC_SAVEDATA_LEN;

#ifdef DATA_FORMAT_GAPSIO
	Gaps_IO_Type TYPE=GAPS_IO_FLOAT64;
	//PTC.dat
	long PTC_Dim=2, PTC_Dimarray[2]={pInputs->num_total_particles,len_data_per_ptc};
	GAPS_IO_InitDataInfo(&(pOutFiles->PTC),TYPE,PTC_Dim,PTC_Dimarray);

	char fnamePTC[GAPS_APT_CONST_DIRSTR_LEN];
	GAPS_APT_GenOutputFileName(fnamePTC,DIR_Output,"PTC.dat");
	GAPS_IO_InitOFile(&(pOutFiles->PTC),fnamePTC);

	long EB_Dimarray[2]={pInputs->num_total_particles,3};
	long Work_Dimarray[2]={pInputs->num_total_particles,1};

	//E.dat
	if(pInputs->EMField_Cal_E)
	{
		GAPS_IO_InitDataInfo(&(pOutFiles->E),TYPE,PTC_Dim,EB_Dimarray);
		char fnameF[GAPS_APT_CONST_DIRSTR_LEN];
		GAPS_APT_GenOutputFileName(fnameF,DIR_Output,"E.dat");
		GAPS_IO_InitOFile(&(pOutFiles->E),fnameF);

		if(pInputs->Open_Cal_Work)
		{
			//Ework.dat
			GAPS_IO_InitDataInfo(&(pOutFiles->Work_E),TYPE,PTC_Dim,Work_Dimarray);

			char fnameWorkF[GAPS_APT_CONST_DIRSTR_LEN];
			GAPS_APT_GenOutputFileName(fnameWorkF,DIR_Output,"Work_E.dat");
			GAPS_IO_InitOFile(&(pOutFiles->Work_E),fnameWorkF);
		}
	}

	//B.dat
	if(pInputs->EMField_Cal_B)
	{
		GAPS_IO_InitDataInfo(&(pOutFiles->B),TYPE,PTC_Dim,EB_Dimarray);

		char fnameF[GAPS_APT_CONST_DIRSTR_LEN];
		GAPS_APT_GenOutputFileName(fnameF,DIR_Output,"B.dat");
		GAPS_IO_InitOFile(&(pOutFiles->B),fnameF);
	}

	//Extern forces
	Gaps_APT_Particle *ptc0=pPtcgrp->ptcgrp;
	pOutFiles->num_forces=ptc0->num_forces;
	pOutFiles->Force=(Gaps_IO_DataFile *)calloc(pOutFiles->num_forces,sizeof(Gaps_IO_DataFile));
	pOutFiles->Work_Force=(Gaps_IO_DataFile *)calloc(pOutFiles->num_forces,sizeof(Gaps_IO_DataFile));
	pOutFiles->ForceName=(char **)calloc(pOutFiles->num_forces,sizeof(char *));
	pOutFiles->ForceType=(int *)calloc(pOutFiles->num_forces,sizeof(int));
	int i;
	long F_Dimarray[2]={pInputs->num_total_particles,GAPS_APT_CONST_EXTERN_FORCE_DIM};
	for(i=0;i<pOutFiles->num_forces;i++)
	{
		pOutFiles->ForceName[i] = ptc0->ForceName[i];
		pOutFiles->ForceType[i] = ptc0->ForceType[i];

		//Forces
		GAPS_IO_InitDataInfo(pOutFiles->Force+i,TYPE,2,F_Dimarray);
		char F_filename[100]="";
		strcat(F_filename,pOutFiles->ForceName[i]);
		strcat(F_filename,".dat");

		char fnameF[GAPS_APT_CONST_DIRSTR_LEN];
		GAPS_APT_GenOutputFileName(fnameF,DIR_Output,F_filename);

		GAPS_IO_InitOFile(pOutFiles->Force+i,fnameF);

		//Force Work
		if(pInputs->Open_Cal_Work)
		{
			GAPS_IO_InitDataInfo(pOutFiles->Work_Force+i,TYPE,2,Work_Dimarray);
			char Fwork_filename[100]="Work_";
			strcat(Fwork_filename,pOutFiles->ForceName[i]);
			strcat(Fwork_filename,".dat");

			char fnameWorkF[GAPS_APT_CONST_DIRSTR_LEN];
			GAPS_APT_GenOutputFileName(fnameWorkF,DIR_Output,Fwork_filename);

			GAPS_IO_InitOFile(pOutFiles->Work_Force+i,fnameWorkF);
		}
	}

#elif DATA_FORMAT_HDF5
	//Create Hdf5 file
	pOutFiles->offset_time=0;

	char filename[GAPS_APT_CONST_DIRSTR_LEN];
	GAPS_APT_GenOutputFileName(filename,DIR_Output,"Data.h5");

	hid_t DataType = H5T_IEEE_F64LE;
	pOutFiles->DataType=DataType;
	int Rank=2;
	pOutFiles->Rank=Rank;
	long num_steps_saved= (pInputs->num_steps)/(pInputs->SavePerNSteps);
	int DimInfo_PTC[2]={len_data_per_ptc*(pInputs->num_total_particles),num_steps_saved};
	int DimInfo_1D[2]={(pInputs->num_total_particles),num_steps_saved};
	int DimInfo_3D[2]={3*(pInputs->num_total_particles),num_steps_saved};

	pOutFiles->File = GAPS_IO_HDF5_CreateFile(filename);

	//PTC_Data
	pOutFiles->PTC=GAPS_IO_HDF5_CreateDataset(pOutFiles->File,"PTC",DataType,Rank,DimInfo_PTC);

	//E.dat
	if(pInputs->EMField_Cal_E)
	{
		pOutFiles->E=GAPS_IO_HDF5_CreateDataset(pOutFiles->File,"E",DataType,Rank,DimInfo_3D);

		if(pInputs->Open_Cal_Work)
		{
			//Ework.dat
			pOutFiles->Work_E=GAPS_IO_HDF5_CreateDataset(pOutFiles->File,"Work_E",DataType,Rank,DimInfo_1D);
		}
	}

	//B.dat
	if(pInputs->EMField_Cal_B)
	{
		pOutFiles->B=GAPS_IO_HDF5_CreateDataset(pOutFiles->File,"B",DataType,Rank,DimInfo_3D);
	}

	//Extern forces
	Gaps_APT_Particle *ptc0=pPtcgrp->ptcgrp;
	pOutFiles->num_forces=ptc0->num_forces;
	pOutFiles->Force=(hid_t *)calloc(pOutFiles->num_forces,sizeof(hid_t));
	pOutFiles->Work_Force=(hid_t *)calloc(pOutFiles->num_forces,sizeof(hid_t));
	pOutFiles->ForceName=(char **)calloc(pOutFiles->num_forces,sizeof(char *));
	pOutFiles->ForceType=(int *)calloc(pOutFiles->num_forces,sizeof(int));
	int i;
	int DimInfo_F[2]={GAPS_APT_CONST_EXTERN_FORCE_DIM*(pInputs->num_total_particles),num_steps_saved};
	for(i=0;i<pOutFiles->num_forces;i++)
	{
		pOutFiles->ForceName[i] = ptc0->ForceName[i];
		pOutFiles->ForceType[i] = ptc0->ForceType[i];

		//Forces
		pOutFiles->Force[i]=GAPS_IO_HDF5_CreateDataset(pOutFiles->File,pOutFiles->ForceName[i],DataType,Rank,DimInfo_F);

		//Force Work
		if(pInputs->Open_Cal_Work)
		{
			char Fwork_filename[100]="Work_";
			strcat(Fwork_filename,pOutFiles->ForceName[i]);
			pOutFiles->Work_Force[i]=GAPS_IO_HDF5_CreateDataset(pOutFiles->File,Fwork_filename,DataType,Rank,DimInfo_1D);
		}
	}

#endif
	return 0;
}

int GAPS_APT_EraseOutputFiles(Gaps_APT_OutputFiles *pOutFiles,Gaps_IO_InputsContainer *pInputs)
{
#ifdef DATA_FORMAT_GAPSIO
	GAPS_IO_DeleteDataInfo(&(pOutFiles->PTC));	
	if(pInputs->EMField_Cal_E)
	{
		GAPS_IO_DeleteDataInfo(&(pOutFiles->E));	
		if(pInputs->Open_Cal_Work)
		{
			GAPS_IO_DeleteDataInfo(&(pOutFiles->Work_E));	
		}
	}
	if(pInputs->EMField_Cal_B)
	{
		GAPS_IO_DeleteDataInfo(&(pOutFiles->B));	
	}

	int i;
	for(i=0;i<pOutFiles->num_forces;i++)
	{
		GAPS_IO_DeleteDataInfo(pOutFiles->Force+i);	
		if(pInputs->Open_Cal_Work)
		{
			GAPS_IO_DeleteDataInfo(pOutFiles->Work_Force+i);	
		}
	}
	free(pOutFiles->Force);
	free(pOutFiles->Work_Force);
	free(pOutFiles->ForceName);
	free(pOutFiles->ForceType);

#elif DATA_FORMAT_HDF5

	H5Dclose(pOutFiles->PTC);	
	if(pInputs->EMField_Cal_E)
	{
		H5Dclose(pOutFiles->E);	
		if(pInputs->Open_Cal_Work)
		{
			H5Dclose(pOutFiles->Work_E);	
		}
	}
	if(pInputs->EMField_Cal_B)
	{
		H5Dclose(pOutFiles->B);	
	}

	int i;
	for(i=0;i<pOutFiles->num_forces;i++)
	{
		H5Dclose(pOutFiles->Force[i]);	
		if(pInputs->Open_Cal_Work)
		{
			H5Dclose(pOutFiles->Work_Force[i]);	
		}
	}
	H5Fclose(pOutFiles->File);

	free(pOutFiles->Force);
	free(pOutFiles->Work_Force);
	free(pOutFiles->ForceName);
	free(pOutFiles->ForceType);
#endif
	return 0;
}
inline int GAPS_APT_SaveData(Gaps_APT_OutputFiles *pOutFiles,Gaps_APT_ParticleGroup *pPtcgrp,long Step,Gaps_IO_InputsContainer *pInputs)
{
	if(pInputs->OpenDataSaving==1)
	{
/*
//		if(pInputs->SaveMode == 1)
//		{
			Gaps_APT_Particle *pPtc = pPtcgrp -> ptcgrp;
			double *pX=GAPS_APT_GetX3(pPtc);
			double *pP=GAPS_APT_GetP3(pPtc);
			double A[4];
			GAPS_APT_CalA(A,pPtc,pInputs);
//			printf("%e",A[0]);
//			printf("\n");
//			printf("%e",A[1]);
//			printf("\n");
//			printf("%e",A[2]);
//			printf("\n");
			double temp[2];
			double ITV;
			temp[0] = pInputs->PoincareDelta[0];
			temp[1] = pInputs->PoincareDelta[1];
			ITV = pInputs->PoincareInterval;
			//if(fabs(pX[1] - temp[0])< ITV && (pP[1]+A[1]-temp[1]) < ITV  && fabs(pX[2])< ITV && fabs(pP[2]+A[2]) < ITV)
			if(fabs(pX[1] - temp[0])< ITV && fabs(pX[2])< ITV && fabs(pP[2]+A[2]) < ITV)
			{
				GAPS_APT_OutPutParticle(pOutFiles,pPtcgrp,pInputs);
			}
		}
*/
		{
			if(Step%(pInputs->SavePerNSteps)==0)
			{
				GAPS_APT_OutPutParticle(pOutFiles,pPtcgrp,pInputs);
			}	
		
		}

	}
	return 0;
}

int GAPS_APT_OutPutParticle(Gaps_APT_OutputFiles *pOutFiles,Gaps_APT_ParticleGroup *pPtcgrp,Gaps_IO_InputsContainer *pInputs)
{
	//Send Particle number
	int flag=0;
	GAPS_APT_SendParticleData2SaveRank(pPtcgrp,0,flag);	

	//Recv Particle number
	long *num_ptc=(long *)calloc(pPtcgrp->proc_size,sizeof(long));
	if((pPtcgrp->proc_rank)==0)
	{	
		MPI_Status status;	
		int rank_i;
		int tag;
		for(rank_i=0;rank_i<pPtcgrp->proc_size;rank_i++)
		{
			tag=rank_i;
			MPI_Recv(&num_ptc[rank_i],1,MPI_LONG,rank_i,tag,pPtcgrp->mpi_comm,&status);
		}
	}

	MPI_Barrier(pPtcgrp->mpi_comm);

	//Send Particle data to rank 0 
	flag=1;
	if((pPtcgrp->proc_rank)!=0)
	{
		GAPS_APT_SendParticleData2SaveRank(pPtcgrp,0,flag);	
	}

	//Recv Particle data and write
	if((pPtcgrp->proc_rank)==0)
	{
		long offset_ptc=0;
		MPI_Status status;	
		int rank_i,tag;
		for(rank_i=0;rank_i<(pPtcgrp->proc_size);rank_i++)
		{
			tag=rank_i;
			long data_len=num_ptc[rank_i]*(pPtcgrp->ptc_data_total_len);

			double *recv_buffer=(double *)calloc(data_len,sizeof(double));
			if(rank_i==0)
			{
				Gaps_APT_Particle *pPtc;
				double *tmp_pointer;
				long i;
				for(i=0;i<num_ptc[rank_i];i++)
				{
					pPtc = pPtcgrp->ptcgrp + i;
					tmp_pointer=recv_buffer+i*(pPtcgrp->ptc_data_total_len);
					memcpy(tmp_pointer,pPtc->data,sizeof(double)*(pPtcgrp->ptc_data_total_len));
				}
			}
			else
			{
				MPI_Recv(recv_buffer,data_len,pPtcgrp->mpi_datatype,rank_i,tag,pPtcgrp->mpi_comm,&status);
				//int tt;
				//for(tt=0;tt<data_len;tt++)
				//		printf("rb=%e\n",recv_buffer[tt]);
			}

			GAPS_APT_WritePtcgrpData2Disk(pOutFiles,recv_buffer,num_ptc[rank_i],pPtcgrp->ptc_data_total_len,offset_ptc,pInputs);
			free(recv_buffer);
			offset_ptc+=num_ptc[rank_i];
		}
	}

#ifdef DATA_FORMAT_HDF5
	(pOutFiles->offset_time)++;
#endif

	free(num_ptc);
	return 0;
}

int GAPS_APT_SendParticleData2SaveRank(Gaps_APT_ParticleGroup *pPtcgrp,int SaveRank,int flag)
{
	int Rank=pPtcgrp->proc_rank;
	int tag=Rank;
	if(0==flag)
	{
		MPI_Send(&(pPtcgrp->num_local_ptc),1,MPI_LONG,SaveRank,tag,pPtcgrp->mpi_comm);
	}
	else
	{
		long num_ptc=pPtcgrp->num_local_ptc;
		Gaps_APT_Particle *pPtc;

		long data_len=(pPtcgrp->ptc_data_total_len)*num_ptc;
		double *send_buffer=(double *)calloc(data_len,sizeof(double));
		double *tmp_pointer;
		long i;
		for(i=0;i<num_ptc;i++)
		{
			pPtc = pPtcgrp->ptcgrp + i;
			tmp_pointer=send_buffer+i*(pPtcgrp->ptc_data_total_len);
			memcpy(tmp_pointer,pPtc->data,sizeof(double)*(pPtcgrp->ptc_data_total_len));
		}
		MPI_Send(send_buffer,data_len,pPtcgrp->mpi_datatype,SaveRank,tag,pPtcgrp->mpi_comm);
		free(send_buffer);
	}
	return 0;
}

int GAPS_APT_WritePtcgrpData2Disk(Gaps_APT_OutputFiles *pOutFiles,double *data,long num_ptc,int data_per_ptc,long offset_ptc,Gaps_IO_InputsContainer *pInputs)
{
	long len_data_per_ptc=GAPS_APT_CONST_PTC_SAVEDATA_LEN;
	long i;
	for(i=0;i<num_ptc;i++)
	{
		double *AllPtcData= data+i*data_per_ptc;	
		double *die=GAPS_APT_GetDie_FromData(AllPtcData);
		double *pT=GAPS_APT_GetT1_FromData(AllPtcData);
		double *pX=GAPS_APT_GetX3_FromData(AllPtcData);
		double *pP=GAPS_APT_GetP3_FromData(AllPtcData);
		double *pAclr=GAPS_APT_GetAclr3_FromData(AllPtcData);

#ifdef DATA_FORMAT_GAPSIO
		GAPS_IO_FWrite(&(pOutFiles->PTC),die,1);
		GAPS_IO_FWrite(&(pOutFiles->PTC),pT,1);
		GAPS_IO_FWrite(&(pOutFiles->PTC),pX,3);
		GAPS_IO_FWrite(&(pOutFiles->PTC),pP,3);
		GAPS_IO_FWrite(&(pOutFiles->PTC),pAclr,3);

		if(pInputs->EMField_Cal_E)
		{
			double *pE=GAPS_APT_GetE3_FromData(AllPtcData);
			GAPS_IO_FWrite(&(pOutFiles->E),pE,3);
			if(pInputs->Open_Cal_Work)
			{
				double *pEwork=GAPS_APT_GetEwork1_FromData(AllPtcData);
				GAPS_IO_FWrite(&(pOutFiles->Work_E),pEwork,1);
			}
		}

		if(pInputs->EMField_Cal_B)
		{
			double *pB=GAPS_APT_GetB3_FromData(AllPtcData);
			GAPS_IO_FWrite(&(pOutFiles->B),pB,3);
		}

		//Output forces
		int j;
		for(j=0;j<pOutFiles->num_forces;j++)
		{
			//Forces
			double *pF=GAPS_APT_GetForceData_FromData(AllPtcData,j);
			GAPS_IO_FWrite(&(pOutFiles->Force[j]),pF,GAPS_APT_CONST_EXTERN_FORCE_DIM);

			//Force Work
			if(pInputs->Open_Cal_Work)
			{
				double *pFwork=GAPS_APT_GetForceWorkData_FromData(AllPtcData,j,pOutFiles->num_forces);
				GAPS_IO_FWrite(&(pOutFiles->Work_Force[j]),pFwork,1);
			}
		}

#elif DATA_FORMAT_HDF5
		int Offset_PTC[2]={(i+offset_ptc)*len_data_per_ptc,pOutFiles->offset_time};
		int Count_PTC[2]={len_data_per_ptc,1};
		int Offset_3D[2]={(i+offset_ptc)*3,pOutFiles->offset_time};
		int Count_3D[2]={3,1};
		int Offset_1D[2]={(i+offset_ptc),pOutFiles->offset_time};
		int Count_1D[2]={1,1};

		double ptc_cache[len_data_per_ptc];
		memcpy(ptc_cache+0,die,sizeof(double)*1);
		memcpy(ptc_cache+1,pT,sizeof(double)*1);
		memcpy(ptc_cache+2,pX,sizeof(double)*3);
		memcpy(ptc_cache+5,pP,sizeof(double)*3);
		memcpy(ptc_cache+8,pAclr,sizeof(double)*3);

		GAPS_IO_HDF5_WriteData(pOutFiles->PTC, pOutFiles->DataType, pOutFiles->Rank,Count_PTC, Offset_PTC, ptc_cache);

		if(pInputs->EMField_Cal_E)
		{
			double *pE=GAPS_APT_GetE3_FromData(AllPtcData);
			GAPS_IO_HDF5_WriteData(pOutFiles->E, pOutFiles->DataType, pOutFiles->Rank,Count_3D, Offset_3D, pE);
			if(pInputs->Open_Cal_Work)
			{
				double *pEwork=GAPS_APT_GetEwork1_FromData(AllPtcData);
				GAPS_IO_HDF5_WriteData(pOutFiles->Work_E, pOutFiles->DataType, pOutFiles->Rank,Count_1D, Offset_1D, pEwork);
			}
		}

		if(pInputs->EMField_Cal_B)
		{
			double *pB=GAPS_APT_GetB3_FromData(AllPtcData);
			GAPS_IO_HDF5_WriteData(pOutFiles->B, pOutFiles->DataType, pOutFiles->Rank,Count_3D, Offset_3D, pB);
		}

		//Output forces
		int j;
		for(j=0;j<pOutFiles->num_forces;j++)
		{
			//Forces
			double *pF=GAPS_APT_GetForceData_FromData(AllPtcData,j);
			GAPS_IO_HDF5_WriteData(pOutFiles->Force[j], pOutFiles->DataType, pOutFiles->Rank,Count_3D, Offset_3D, pF);

			//Force Work
			if(pInputs->Open_Cal_Work)
			{
				double *pFwork=GAPS_APT_GetForceWorkData_FromData(AllPtcData,j,pOutFiles->num_forces);
				GAPS_IO_HDF5_WriteData(pOutFiles->Work_Force[j], pOutFiles->DataType, pOutFiles->Rank,Count_1D, Offset_1D, pFwork);
			}
		}
#endif
	}

	return 0;
}

int GAPS_APT_CommPtcNum_Disp(Gaps_APT_ParticleGroup *pPtcgrp)
{
	if(pPtcgrp->proc_rank==0)
	{	
		long tmp;
		long num_ptc = pPtcgrp->count_ptc;
		MPI_Status status;	
		int rank_i;
		printf("0: %ld\t",num_ptc);
		for(rank_i=1;rank_i<pPtcgrp->proc_size;rank_i++)
		{
			long tag=rank_i;
			MPI_Recv(&tmp,1,MPI_LONG,rank_i,tag,pPtcgrp->mpi_comm,&status);
			num_ptc+=tmp;
			printf("%d: %ld\t",rank_i,tmp);
		}
		printf("\n");
		printf("PTC_NUM = %ld\n",num_ptc);
	}
	else
	{
		MPI_Send(&(pPtcgrp->count_ptc),1,MPI_LONG,0,pPtcgrp->proc_rank,pPtcgrp->mpi_comm);
	}
	return 0;
}

//Hdf5 functions
#ifdef DATA_FORMAT_HDF5

hid_t GAPS_IO_HDF5_CreateFile(char *FileName)
{
	return H5Fcreate (FileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
}

hid_t GAPS_IO_HDF5_OpenFile(char *FileName)
{
	return H5Fopen(FileName, H5F_ACC_RDWR, H5P_DEFAULT);
}

hid_t GAPS_IO_HDF5_CreateGroup(hid_t fileID, char *GroupName)
{
	hid_t groupID;
	groupID=H5Gcreate2(fileID,GroupName,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	return groupID;
}

hid_t GAPS_IO_HDF5_CreateDataset(hid_t fileORgroupID, char *DatasetName, hid_t DataType , int Rank, int *DimInfo)
{//DimInfo is a 1D-array: DimInfo[rank+1]={rank, num_1, num_2, ... , num_rank}, rank is the dimention of data, num_1~num_rank is the number of each dimention.

	hid_t filespace, datasetID;
	hsize_t dimsf[Rank];

	int i;
	for(i=0;i<Rank;i++)
	{
		(*(dimsf+i)) = (*(DimInfo+i));
	}

	filespace = H5Screate_simple(Rank,dimsf,NULL);

	datasetID = H5Dcreate(fileORgroupID,DatasetName,DataType,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

	H5Sclose(filespace);
	return datasetID;
}

int GAPS_IO_HDF5_WriteData(hid_t datasetID, hid_t DataType, int Rank, int *pCount,int *Offset, void *data)
{
	herr_t status;
	hsize_t count[Rank], offset[Rank];

	int i;
	for (i=0;i<Rank;i++)
	{
		(*(count+i)) = (*(pCount+i));
		(*(offset+i)) = (*(Offset+i));
	}

	hid_t filespace, memspace;
	memspace = H5Screate_simple(Rank,count,NULL);

	filespace = H5Dget_space(datasetID);
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET,offset,NULL,count,NULL);


	status = H5Dwrite(datasetID, DataType,memspace,filespace,H5P_DEFAULT,data);

	H5Sclose(filespace);
	H5Sclose(memspace);
	return 0;
}
#endif

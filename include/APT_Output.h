#define GAPS_APT_CONST_PTC_SAVEDATA_LEN 11
#define GAPS_APT_CONST_DIRSTR_LEN 1024

#ifdef DATA_FORMAT_GAPSIO

typedef struct{
	Gaps_IO_DataFile PTC;//Particle data--Output: die,t,x,y,z,gamma,px,py,pz,ax,ay,az
	Gaps_IO_DataFile E;//Electric field--Output: 
	Gaps_IO_DataFile B;//Magnetic field--Output: 
	Gaps_IO_DataFile Work_E;//Electric field work--Output: 

	Gaps_IO_DataFile *Force;//Output Force data: runtime setting
	Gaps_IO_DataFile *Work_Force;//Output Force work data: runtime setting
	int num_forces;
	char **ForceName;
	int *ForceType;
}Gaps_APT_OutputFiles;

#elif DATA_FORMAT_HDF5

typedef struct{
	hid_t File;
	long offset_time;
	hid_t DataType;
	int Rank;

	hid_t PTC;
	hid_t E;
	hid_t B;
	hid_t Work_E;

	hid_t *Force;//Output Force data: runtime setting
	hid_t *Work_Force;//Output Force work data: runtime setting

	int num_forces;
	char **ForceName;
	int *ForceType;
}Gaps_APT_OutputFiles;

#endif

int GAPS_APT_GenOutputFileName(char *pDirOut,char *pDirIN,char *filename);
int GAPS_APT_ProducePath(char *pDirOut,char *pDirIN);
int GAPS_APT_InitOutputFiles(Gaps_APT_OutputFiles *pOutFiles,Gaps_APT_ParticleGroup *pPtcgrp,char *DIR_Output,Gaps_IO_InputsContainer *pInputs);
int GAPS_APT_EraseOutputFiles(Gaps_APT_OutputFiles *pOutFiles,Gaps_IO_InputsContainer *pInputs);
inline int GAPS_APT_SaveData(Gaps_APT_OutputFiles *pOutFiles,Gaps_APT_ParticleGroup *pPtcgrp,long Step,Gaps_IO_InputsContainer *pInputs);
int GAPS_APT_OutPutParticle(Gaps_APT_OutputFiles *pOutFiles,Gaps_APT_ParticleGroup *pPtcgrp,Gaps_IO_InputsContainer *pInputs);
int GAPS_APT_SendParticleData2SaveRank(Gaps_APT_ParticleGroup *pPtcgrp,int SaveRank,int flag);
int GAPS_APT_WritePtcgrpData2Disk(Gaps_APT_OutputFiles *pOutFiles,double *data,long num_ptc,int data_per_ptc,long offset_ptc,Gaps_IO_InputsContainer *pInputs);
int GAPS_APT_CommPtcNum_Disp(Gaps_APT_ParticleGroup *pPtcgrp);


//Hdf5 functions
#ifdef DATA_FORMAT_HDF5
hid_t GAPS_IO_HDF5_CreateFile(char *FileName);
hid_t GAPS_IO_HDF5_OpenFile(char *FileName);
hid_t GAPS_IO_HDF5_CreateGroup(hid_t fileID, char *GroupName);
hid_t GAPS_IO_HDF5_CreateDataset(hid_t fileORgroupID, char *DatasetName, hid_t DataType , int Rank, int *DimInfo);
int GAPS_IO_HDF5_WriteData(hid_t datasetID, hid_t DataType, int Rank, int *pCount,int *Offset, void *data);
#endif

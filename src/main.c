#include "APT_AllHeaders.h"

#define GAPS_APT_CONST_NUM_MANAGE_LUA_SCRIPTS 2
#define GAPS_APT_CONST_NUM_ARGUMENTS 2
char *ArgOptionList[GAPS_APT_CONST_NUM_ARGUMENTS]={
						"-c",// 0: Directory of Configuration Files
						"-o" // 1: Directory of Output Data
					};

char GAPS_APT_DIR_Input[GAPS_APT_CONST_DIRSTR_LEN];
char GAPS_APT_DIR_Output[GAPS_APT_CONST_DIRSTR_LEN];

Gaps_APT_DiscreteTensor *Global_pDisTensor;
Gaps_APT_Particle **Global_pCurrentPtc;
Gaps_APT_ParticleGroup *Global_pPtcgrp;

Gaps_IO_InputsContainer Inputs;

int main(int argc,char *argv[])
{
	//Set arguments from bash command
	char *Arguments[GAPS_APT_CONST_NUM_ARGUMENTS]={NULL};
	RT_ARG_GetParameterSeq(Arguments,argc,argv,ArgOptionList,GAPS_APT_CONST_NUM_ARGUMENTS);
	GAPS_APT_ProducePath(GAPS_APT_DIR_Input,Arguments[0]);
	GAPS_APT_ProducePath(GAPS_APT_DIR_Output,Arguments[1]);
	//Create output directory
	if(-1==access(GAPS_APT_DIR_Output,0))
	{
		mkdir(GAPS_APT_DIR_Output,0775);
	}

	if(-1==access(GAPS_APT_DIR_Input,0))
	{
		fprintf(stderr,"ERROR: The directory of configuration files is NOT exist!! (Option: -c)\n");
		assert(-1!=access(GAPS_APT_DIR_Input,0));
	}

	//Load Inputs
	int num_lua_files=GAPS_APT_CONST_NUM_MANAGE_LUA_SCRIPTS;	
	char **LuaFiles=(char **)malloc(num_lua_files*sizeof(char *));
	char LuaSetDir[GAPS_APT_CONST_DIRSTR_LEN],LuaRootDir[GAPS_APT_CONST_DIRSTR_LEN];
	strcpy(LuaSetDir,GAPS_APT_DIR_Input);strcat(LuaSetDir,"pkg/lib/Setting.lua");
	strcpy(LuaRootDir,GAPS_APT_DIR_Input);strcat(LuaRootDir,"Config.lua");
	LuaFiles[0] = LuaSetDir;
	LuaFiles[1] = LuaRootDir;
	Gaps_IO_LuaInputEnv LuaEnv;
	GAPS_IO_Load_ConfigEnvironment_MultiFile(&LuaEnv,LuaFiles,num_lua_files);	
	GAPS_IO_LoadLua2C(&LuaEnv,&Inputs);
	//Initial parallelization
	int proc_rank,proc_size;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&proc_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&proc_size);

/*
	
	int i = 0;
	char hostname[256];
	gethostname(hostname, sizeof(hostname));
	printf("PID %d on %s ready for attach\n", getpid(), hostname);
	fflush(stdout);
	while (0 == i)
	sleep(5);

*/

	Gaps_APT_DiscreteTensor DisTensor;
	//Set Discrete field data
	if(Inputs.EMField_Type == -1)
	{
		Global_pDisTensor = &DisTensor;	
		GAPS_APT_Field_Discrete_LoadData(&DisTensor,Inputs.EMField_Discrete_Filename);
	}
	
	//Set Particle pusher
	Gaps_APT_ParticlePusher PusherUsing;
	GAPS_APT_SetParticlePusher(&PusherUsing,&Inputs);
	if(0 == proc_rank){
		printf("PushType is %d\n",Inputs.Pusher_Type);
		printf("EMFieldType is %d\n",Inputs.EMField_Type);
	}
	//Initialize Particle structure: assign RAM, set fields and forces
	Gaps_APT_Particle *CurrentPtc;
	Global_pCurrentPtc =&CurrentPtc;
	Gaps_APT_ParticleGroup Ptcgrp;
	Global_pPtcgrp = &Ptcgrp;// Assign address of Ptcgrp to global variable
	GAPS_APT_InitParticleGroup(&Ptcgrp,&Inputs,proc_rank,proc_size,MPI_DOUBLE,MPI_COMM_WORLD);

	//Initialize Particle Initial Condition
	srand(1+proc_rank);
	GAPS_APT_ParticleInitialCondition(&Ptcgrp,&Inputs);

	//Generate Configuration Information MATLAB file
	if(proc_rank == 0)
	{
		char filename_info_matlab[GAPS_APT_CONST_DIRSTR_LEN],filename_info_python[GAPS_APT_CONST_DIRSTR_LEN];
		strcpy(filename_info_matlab,GAPS_APT_DIR_Output);strcat(filename_info_matlab,"LoadCalInfo.m");
		strcpy(filename_info_python,GAPS_APT_DIR_Output);strcat(filename_info_python,"LoadCalInfo.py");
		GAPS_IO_GenCalInfoMfile(filename_info_matlab,&Inputs);
		GAPS_IO_GenCalInfoPython(filename_info_python,&Inputs);
	}
	//Initialize Output files	#####################
	Gaps_APT_OutputFiles OutputFiles;
	if(proc_rank == 0)
	{
		GAPS_APT_InitOutputFiles(&OutputFiles,&Ptcgrp,GAPS_APT_DIR_Output,&Inputs);
	}

	//Iteration	#####################
	double t=0;
	long Step,idx_ptc;
	double iter_start_mpi,iter_finish_mpi;
	iter_start_mpi=MPI_Wtime();

	for(Step=0;Step<Inputs.num_steps;Step++)
	{
		if(0==proc_rank)
		{
		//	printf("##\tstep = %ld\n",Step);
		}

		GAPS_APT_SaveData(&OutputFiles,&Ptcgrp,Step,&Inputs);

		// Set particle source: add particles or delete particles
		//	GAPS_APT_ParticleSource(&Ptcgrp,&Inputs);

		// Push particles
		GAPS_APT_Push_ParticleGroup(&Ptcgrp,PusherUsing,&Inputs);
	}

	iter_finish_mpi=MPI_Wtime();
	fprintf(stdout,"Rank:%d, Iteration Time: %e s\n",proc_rank,iter_finish_mpi-iter_start_mpi);
	//Free RAM
	GAPS_APT_EraseParticleGroup(&Ptcgrp);
	if(proc_rank==0)
	{
		GAPS_APT_EraseOutputFiles(&OutputFiles,&Inputs);
	}
	if(Inputs.EMField_Type == -1)
	{
		GAPS_APT_EraseDiscreteTensor(&DisTensor);
	}
	//MPI End
	MPI_Finalize();	
	free(LuaFiles);
	return 0;
}

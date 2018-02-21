
	#ifndef   GAPS_IO_LIB_VERS_0001    
		
#define GAPS_IO_LIB_VERS_0001
#include <stdint.h>
#include <unistd.h>
#include <sys/types.h>
typedef struct { 	int64_t  type ;
	int64_t  dim ;
	int64_t *  pdimarray ;
	int64_t  numperstep ;
	FILE *  pfile ;
} Gaps_IO_DataFile;
typedef enum {GAPS_IO_CONST_INT32,GAPS_IO_INT64,GAPS_IO_INT16,GAPS_IO_UINT32,GAPS_IO_UINT64,GAPS_IO_UINT16,GAPS_IO_FLOAT32,GAPS_IO_FLOAT64,GAPS_IO_INT128,GAPS_IO_UINT128,GAPS_IO_FLOAT128,GAPS_IO_INT8,GAPS_IO_UINT8} Gaps_IO_Type;

extern int* GAPS_IO_GlobalTypeLen;int  GAPS_IO_TruncateFile (Gaps_IO_DataFile *  pOutfile ,int64_t  NumSteps ,int64_t  NumData );
int  GAPS_IO_DeleteDataInfo (Gaps_IO_DataFile *  pOutfile );
int  GAPS_IO_DataSeek (Gaps_IO_DataFile *  pOutfile ,int64_t  Step ,int64_t  Offset );
int  GAPS_IO_DataNumStepsAndResidue (Gaps_IO_DataFile *  pOutfile ,int64_t *  pNumsteps ,int64_t *  pResidue );
int  GAPS_IO_FWrite (Gaps_IO_DataFile *  pOutfile ,void *  pData ,int64_t  NumData );
int  GAPS_IO_FRead (Gaps_IO_DataFile *  pOutfile ,void *  pData ,int64_t  NumData );
int  GAPS_IO_InitIFile (Gaps_IO_DataFile *  pOutfile ,char *  pName );
int  GAPS_IO_FileGotoBegin (Gaps_IO_DataFile *  pOutfile );
int  GAPS_IO_InitIFilePointer (Gaps_IO_DataFile *  pOutfile ,FILE *  fp );
int  GAPS_IO_InitOFile (Gaps_IO_DataFile *  pOutfile ,char *  pName );
int  GAPS_IO_InitDataInfo (Gaps_IO_DataFile *  pOutfile ,int64_t  Type ,int64_t  Dim ,int64_t *  pDimarray );

	#else
		
	 #endif


// input_lua
typedef struct {
	lua_State *configenv;
}Gaps_IO_LuaInputEnv ;

int GAPS_IO_Load_ConfigEnvironment_MultiFile(Gaps_IO_LuaInputEnv *pLuaenv, char **pFiles,int num_files);
int GAPS_IO_Load_ConfigEnvironment(Gaps_IO_LuaInputEnv *pLuaenv, char *pFilename);
int GAPS_IO_Load_double(Gaps_IO_LuaInputEnv *pLuaenv,char *pVarname,double *pTarget);
int GAPS_IO_Load_long(Gaps_IO_LuaInputEnv *pLuaenv,char *pVarname,long *pTarget);
int GAPS_IO_Load_string(Gaps_IO_LuaInputEnv *pLuaenv,char *pVarname,char *pTarget);
int GAPS_IO_Load_table(Gaps_IO_LuaInputEnv *pLuaenv,char *pVarname,long Dims,double *pTarget);
int GAPS_IO_Load_function(Gaps_IO_LuaInputEnv *pLuaenv, char *pFuncName, long num_inputs, double *inputs, long num_outputs, double *outputs);

long GAPS_IO_Load_LenTable(Gaps_IO_LuaInputEnv *pLuaenv,char *pName);
double GAPS_IO_Load_InitTable(Gaps_IO_LuaInputEnv *pLuaenv,int index,double key);


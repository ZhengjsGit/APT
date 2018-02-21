
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


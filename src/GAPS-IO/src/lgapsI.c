#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "lgapsI.h"

int GAPS_IO_Load_ConfigEnvironment_MultiFile(Gaps_IO_LuaInputEnv *pLuaenv, char **pFiles,int num_files)
{
	pLuaenv->configenv = luaL_newstate();
	lua_State *L=pLuaenv->configenv;
	luaopen_base(L);
	luaL_openlibs(L);
	luaopen_io(L);
	luaopen_string(L);
	luaopen_math(L);	
	int i;
	for(i=0;i<num_files;i++)
	{
		if (luaL_loadfile(L, pFiles[i]) || lua_pcall(L,0,0,0))
		{
			error(L, "Cannot run configuration file: %s", lua_tostring(L, -1));
		}
	}
	return 0;
}
int GAPS_IO_Load_ConfigEnvironment(Gaps_IO_LuaInputEnv *pLuaenv,char *pFilename)
{
	pLuaenv->configenv = luaL_newstate();
	lua_State *L=pLuaenv->configenv;
	luaopen_base(L);
	luaL_openlibs(L);
	luaopen_io(L);
	luaopen_string(L);
	luaopen_math(L);	
	if (luaL_loadfile(L, pFilename) || lua_pcall(L,0,0,0))
	{
		error(L, "Cannot run configuration file: %s", lua_tostring(L, -1));
	}
	return 0;
}

int GAPS_IO_Load_double(Gaps_IO_LuaInputEnv *pLuaenv,char *pVarname,double *pTarget)
{
	lua_State *L=pLuaenv->configenv;
	lua_getglobal(L, pVarname);
	
	if (lua_isnumber(L,-1))
	{
		*pTarget= (double)lua_tonumber(L,-1);
		lua_settop(L,0);
	}
	else
	{
		fprintf(stderr, "\"%s\" is not a number or it is not defined in Configuration file!\n",pVarname);
		error(L, "This variable is not a number or it is not defined in `Config.lua'!");
	}
	return 0;
}


int GAPS_IO_Load_long(Gaps_IO_LuaInputEnv *pLuaenv,char *pVarname,long *pTarget)
{
	lua_State *L=pLuaenv->configenv;
	lua_getglobal(L, pVarname);
	
	if (lua_isnumber(L,-1))
	{
		*pTarget= (long)lua_tonumber(L,-1);
		lua_settop(L,0);
	}
	else
	{
		fprintf(stderr, "\"%s\" is not a number or it is not defined in Configuration file!\n",pVarname);
		error(L, "This variable is not a number or it is not defined in `Config.lua'!");
	}

}

int GAPS_IO_Load_string(Gaps_IO_LuaInputEnv *pLuaenv,char *pVarname,char *pTarget)
{
	lua_State *L=pLuaenv->configenv;
	lua_getglobal(L, pVarname);
	
	if (lua_isstring(L,-1))
	{
		*pTarget='\0';
		strcat(pTarget,lua_tostring(L,-1));
		lua_settop(L,0);
	}
	else
	{
		fprintf(stderr, "\"%s\" is not a string or it is not defined in Configuration file!\n",pVarname);
		error(L, "This variable is not a string or it is not defined in `Config.lua'!");
	}

}

int GAPS_IO_Load_table(Gaps_IO_LuaInputEnv *pLuaenv,char *pVarname,long Dims,double *pTarget)
{
	lua_State *L=pLuaenv->configenv;
	long i,tabledims;
	tabledims=GAPS_IO_Load_LenTable(pLuaenv,pVarname);

	if (tabledims > Dims)
	{
		printf("Notice! Array-\"%s\" : Length of Inputed table is longer than the declared length of!\n",pVarname);
		tabledims=Dims;
	}
	lua_getglobal(L, pVarname);
	if (lua_istable(L,-1))
	{
		for (i=0;i<tabledims;i++)
		{
			*(pTarget+i) = GAPS_IO_Load_InitTable(pLuaenv, -1, i+1);
			
		}
		lua_settop(L,0);
	}
	else
	{
		fprintf(stderr, "\"%s\" is not a table or it is not defined in Configuration file!\n",pVarname);
		error(L, "This variable is not a table or it is not well defined in `Config.lua'!");
	}

	return tabledims;
}

double GAPS_IO_Load_InitTable(Gaps_IO_LuaInputEnv *pLuaenv, int index, double key)
{
	lua_State *L=pLuaenv->configenv;
	double result;
	lua_pushnumber(L, key);
	lua_gettable(L, index-1);
	if (!lua_isnumber(L, -1))
		error(L, "Invalid component in table");
	result = (double)lua_tonumber(L, -1);
	lua_pop(L,1);
	return result;
}

long GAPS_IO_Load_LenTable(Gaps_IO_LuaInputEnv *pLuaenv, char *pName)
{
	lua_State *L=pLuaenv->configenv;
	lua_getglobal(L,pName);
	long len=0,i;
	lua_pushnil(L);
	while (lua_next(L,-2))
	{
		len++;
		lua_pop(L,1);
	}
	lua_settop(L,0);
	return len;
}

int GAPS_IO_Load_function(Gaps_IO_LuaInputEnv *pLuaenv, char *pFuncName, long num_inputs, double *pInputs, long num_outputs, double *pOutputs)
{
	lua_State *L=pLuaenv->configenv;
	int i;
	lua_getglobal(L,pFuncName);
	for (i=0;i<num_inputs;i++)
	{
		lua_pushnumber(L,*(pInputs+i));
	}
	
	if (lua_pcall(L,num_inputs,num_outputs,0) != 0)
		
	{	fprintf(stderr,"Function \"%s\" doesnot exit or has bugs.\n",pFuncName);
		error("running function");
	}
	for (i=0;i<num_outputs;i++)
	{
		(*(pOutputs+num_outputs-i-1)) = (double)lua_tonumber(L,-(i+1));
	}
	lua_settop(L,0);
	return 0;
}

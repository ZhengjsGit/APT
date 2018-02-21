#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<assert.h>
#include<RT_ARG.h>

int RT_ARG_MapParaKey2Num(char *pKey,char *OptionList[],int numOption)
{
	int len=numOption;
	int i;
	for(i=0;i<len;i++)
	{
		if(!strcmp(pKey,OptionList[i]))
		{
			return i;
		}
	}
}

int RT_ARG_IsOptionKey(char *pOpKey,char *OptionList[],int numOption)
{
	int len=numOption;
	int i;
	int count=0;
	for(i=0;i<len;i++)
	{
		count+=(!strcmp(pOpKey,OptionList[i]));
	}
	return count!=0;
}

int RT_ARG_GetParameterSeq(char *pParameters[],int argc,char *argv[],char *OptionList[],int numOption)
{
	int i=0;
	while(i<argc)
	{
		if(RT_ARG_IsOptionKey(argv[i],OptionList,numOption))
		{
			if(i>=argc-1)
			{
				fprintf(stderr,"ERROR: No option is assigned to '%s'\n",argv[i]);
				assert(i<argc-1);
			}
			pParameters[RT_ARG_MapParaKey2Num(argv[i],OptionList,numOption)]=argv[i+1];
			i++;
		}
		i++;
	}
	return 0;
}


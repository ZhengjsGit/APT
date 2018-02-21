#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/types.h>
typedef struct {
	int64_t type ;
	int64_t dim ;
	int64_t *pdimarray ;
	int64_t numperstep ;
	FILE *pfile ;
} Gaps_IO_DataFile;
typedef enum {GAPS_IO_CONST_INT32,GAPS_IO_INT64,GAPS_IO_INT16,GAPS_IO_UINT32,GAPS_IO_UINT64,GAPS_IO_UINT16,GAPS_IO_FLOAT32,GAPS_IO_FLOAT64,GAPS_IO_INT128,GAPS_IO_UINT128,GAPS_IO_FLOAT128,GAPS_IO_INT8,GAPS_IO_UINT8} Gaps_IO_Type;

int GAPS_IO_GlobalTypeLen[]= {4,8,2,4,8,2,4,8,16,16,16,1,1};
int GAPS_IO_InitDataInfo(Gaps_IO_DataFile *pOutfile ,int64_t Type ,int64_t Dim ,int64_t *pDimarray) {
	(pOutfile)->dim = Dim;
	size_t dimarrlen = (sizeof(int64_t) * Dim) ;
	(pOutfile)->pdimarray = malloc(dimarrlen);
	if(((pOutfile)->pdimarray == NULL)) {
		fprintf(stderr , "Error: unable to allocate the memory\n");
		return -1 ;
	} else {
		0;
	}
	memcpy((pOutfile)->pdimarray , pDimarray , dimarrlen);
	(pOutfile)->type = Type;
	int64_t alll = 1 ;
	long g = 0 ;
	for(g=0 ; (g < Dim) ; g++) {
		alll = (alll * (pDimarray)[g]);
	}(pOutfile)->numperstep = alll;
	return 0 ;
}
int GAPS_IO_InitOFile(Gaps_IO_DataFile *pthis ,char *pName) {
	int64_t type = (pthis)->type ;
	int64_t dim = (pthis)->dim ;
	int64_t *pdimarray = (pthis)->pdimarray ;
	int64_t numperstep = (pthis)->numperstep ;
	FILE *pfile = (pthis)->pfile ;
	pfile = fopen(pName , "w+");
	assert(pfile);
	(pthis)->pfile = pfile;
	int64_t version = 0 ;
	assert((1 == fwrite(& (version) , sizeof(int64_t) , 1 , pfile)));
	assert((1 == fwrite(& (type) , sizeof(int64_t) , 1 , pfile)));
	assert((1 == fwrite(& (dim) , sizeof(int64_t) , 1 , pfile)));
	assert((dim == fwrite(pdimarray , sizeof(int64_t) , dim , pfile)));
	return 0 ;
}
int GAPS_IO_InitIFilePointer(Gaps_IO_DataFile *pthis ,FILE *fp) {
	int64_t type = (pthis)->type ;
	int64_t dim = (pthis)->dim ;
	int64_t *pdimarray = (pthis)->pdimarray ;
	int64_t numperstep = (pthis)->numperstep ;
	FILE *pfile = (pthis)->pfile ;
	pfile = fp;
	(pthis)->pfile = pfile;
	int64_t version ;
	int64_t *pversion = & (version) ;
	assert((1 == fread(pversion , sizeof(int64_t) , 1 , pfile)));
	assert((version == 0));
	assert((1 == fread(& (type) , sizeof(int64_t) , 1 , pfile)));
	assert(((type * sizeof(int)) < sizeof(GAPS_IO_GlobalTypeLen)));
	(pthis)->type = type;
	assert((1 == fread(& (dim) , sizeof(int64_t) , 1 , pfile)));
	(pthis)->dim = dim;
	pdimarray = malloc((dim * sizeof(int64_t)));
	assert(pdimarray);
	assert((dim == fread(pdimarray , sizeof(int64_t) , dim , pfile)));
	(pthis)->pdimarray = pdimarray;
	numperstep = 1;
	long i = 0 ;
	for(i=0 ; i<dim ; i++) {
		numperstep = (numperstep * (pdimarray)[i]);
		assert((numperstep > 0));
	}(pthis)->numperstep = numperstep;
	return 0 ;
}
int GAPS_IO_FileGotoBegin(Gaps_IO_DataFile *pthis) {
	int64_t type = (pthis)->type ;
	int64_t dim = (pthis)->dim ;
	int64_t *pdimarray = (pthis)->pdimarray ;
	int64_t numperstep = (pthis)->numperstep ;
	FILE *pfile = (pthis)->pfile ;
	return fseek(pfile , ((3 + dim) * sizeof(int64_t)) , SEEK_SET) ;
}
int GAPS_IO_InitIFile(Gaps_IO_DataFile *pthis ,char *pName) {
	int64_t type = (pthis)->type ;
	int64_t dim = (pthis)->dim ;
	int64_t *pdimarray = (pthis)->pdimarray ;
	int64_t numperstep = (pthis)->numperstep ;
	FILE *pfile = (pthis)->pfile ;
	pfile = fopen(pName , "r+");
	assert(pfile);
	return GAPS_IO_InitIFilePointer(pthis , pfile) ;
}
int GAPS_IO_FRead(Gaps_IO_DataFile *pthis ,void *pData ,int64_t NumData) {
	int64_t type = (pthis)->type ;
	int64_t dim = (pthis)->dim ;
	int64_t *pdimarray = (pthis)->pdimarray ;
	int64_t numperstep = (pthis)->numperstep ;
	FILE *pfile = (pthis)->pfile ;
	if((NumData < 0)) {
		NumData = numperstep;
	} else {
		0;
	}
	return fread(pData , (GAPS_IO_GlobalTypeLen)[type] , NumData , pfile) ;
}
int GAPS_IO_FWrite(Gaps_IO_DataFile *pthis ,void *pData ,int64_t NumData) {
	int64_t type = (pthis)->type ;
	int64_t dim = (pthis)->dim ;
	int64_t *pdimarray = (pthis)->pdimarray ;
	int64_t numperstep = (pthis)->numperstep ;
	FILE *pfile = (pthis)->pfile ;
	if((NumData < 0)) {
		NumData = numperstep;
	} else {
		0;
	}
	return fwrite(pData , (GAPS_IO_GlobalTypeLen)[type] , NumData , pfile) ;
}
int GAPS_IO_DataNumStepsAndResidue(Gaps_IO_DataFile *pthis ,int64_t *pNumsteps ,int64_t *pResidue) {
	int64_t type = (pthis)->type ;
	int64_t dim = (pthis)->dim ;
	int64_t *pdimarray = (pthis)->pdimarray ;
	int64_t numperstep = (pthis)->numperstep ;
	FILE *pfile = (pthis)->pfile ;
	long cur_pos = ftell(pfile) ;
	assert((cur_pos >= 0));
	GAPS_IO_FileGotoBegin(pthis);
	long head_pos=ftell(pfile);
	assert((0 == fseek(pfile , 0 , SEEK_END)));
	long end_pos = ftell(pfile) ;
	assert((end_pos > 0));
	end_pos = ((end_pos-head_pos) / (GAPS_IO_GlobalTypeLen)[type]);
	(pNumsteps)[0] = (end_pos / numperstep);
	(pResidue)[0] = (end_pos % numperstep);
	fseek(pfile , cur_pos , SEEK_SET);
	return 0 ;
}
int GAPS_IO_DataSeek(Gaps_IO_DataFile *pthis ,int64_t Step ,int64_t Offset) {
	int64_t type = (pthis)->type ;
	int64_t dim = (pthis)->dim ;
	int64_t *pdimarray = (pthis)->pdimarray ;
	int64_t numperstep = (pthis)->numperstep ;
	FILE *pfile = (pthis)->pfile ;
	assert((0 == fseek(pfile , (((GAPS_IO_GlobalTypeLen)[type] * (Offset + (Step * numperstep))) + (sizeof(int64_t) * (3 + dim))) , SEEK_SET)));
	return 0 ;
}
int GAPS_IO_DeleteDataInfo(Gaps_IO_DataFile *pthis) {
	int64_t type = (pthis)->type ;
	int64_t dim = (pthis)->dim ;
	int64_t *pdimarray = (pthis)->pdimarray ;
	int64_t numperstep = (pthis)->numperstep ;
	FILE *pfile = (pthis)->pfile ;
	if(pfile) {
		fclose(pfile);
	} else {
		0;
	}
	free(pdimarray);
	return 0 ;
}
int GAPS_IO_TruncateFile(Gaps_IO_DataFile *pthis ,int64_t NumSteps ,int64_t NumData) {
	int64_t type = (pthis)->type ;
	int64_t dim = (pthis)->dim ;
	int64_t *pdimarray = (pthis)->pdimarray ;
	int64_t numperstep = (pthis)->numperstep ;
	FILE *pfile = (pthis)->pfile ;
	if((NumData < 0)) {
		NumData = numperstep;
	} else {
		0;
	}
	return ftruncate(fileno(pfile) , ((NumSteps * (NumData * (GAPS_IO_GlobalTypeLen)[type])) + (sizeof(int64_t) * (3 + dim)))) ;
}

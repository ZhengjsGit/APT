#define GAPS_APT_CONST_EB_TENSOR_DIM 6

typedef int (*Gaps_APT_Field)(double *pTensor,double *pSpaceTime4,int Order, Gaps_IO_InputsContainer *pInputs);

inline int GAPS_APT_PowerInt(int x,int n);
int GAPS_APT_ResetTensor(double *pTensor,int Order);
int GAPS_APT_Field_Discrete(double *pTensor,double *pSpaceTime4,int Order, Gaps_IO_InputsContainer *pInputs);
inline int GAPS_APT_TensorLen(int Order);
int GAPS_APT_Field_Discrete(double *pTensor,double *pSpaceTime4,int Order, Gaps_IO_InputsContainer *pInputs);

int GAPS_APT_GatherPtcInfo(double *pData);

typedef struct {
	//Global Field Parameters
		long xnum;
		long ynum; 
		long znum;
		long tnum;
		long dim;
		double dx;
		double dy;
		double dz;
		double dt;
		long num_data_per_moment;
	
		long order;
		double boundary[8];
		double **data;
		double *data1D;
} Gaps_APT_DiscreteTensor;

int GAPS_APT_Field_Discrete_LoadData(Gaps_APT_DiscreteTensor *pDisTensor,char *pFilename);
int GAPS_APT_EraseDiscreteTensor(Gaps_APT_DiscreteTensor *pDisTensor);
inline int GAPS_APT_Field_Discrete_GetGridIdx(long *pIdx4,long *pRelIdx4,Gaps_APT_DiscreteTensor *pDisTensor,double *pST4);
double * GAPS_APT_GetFieldVector(double *data,long *pWT_idx,Gaps_APT_DiscreteTensor *pDisTensor);
#include "EMField_Prototype.h"


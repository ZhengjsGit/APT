#include "APT_AllHeaders.h"

inline double GAPS_APT_TensorValue(double *pTensor,int *pIndex,int Dim,int Order)
{
	if(0==Order)
	{
		return *pTensor;
	}
	else if(-1 == Order)
	{
		return pTensor[pIndex[0]];
	}
	else if(1<=Order)
	{
		int i;
		int idx=pIndex[Order-1];
		for(i=Order-1;i>0;i--)
		{
			idx=idx*Dim+pIndex[i-1];	// idx = i0+Dim(i1+Dim(i2+Dim(i3 ...)))
		}
		return pTensor[idx];
	}
	else
	{
		assert(Order>=-1);
		fprintf(stderr,"ERROR: In function GAPS_APT_TensorValue: Order must be larger than -1\n");
		return 0;
	}
}

inline double * GAPS_APT_TensorPointer(double *pTensor,int *pIndex,int Dim,int Order)
{
	if(0==Order)
	{
		return pTensor;
	}
	else if(-1 == Order)
	{
		return pTensor+pIndex[0];
	}
	else if(1<=Order)
	{
		int i;
		int idx=pIndex[Order-1];
		for(i=Order-1;i>0;i--)
		{
			idx=idx*Dim+pIndex[i-1];	// idx = i0+Dim(i1+Dim(i2+Dim(i3 ...)))
		}
		return pTensor+idx;
	}
	else
	{
		assert(Order>=-1);
		fprintf(stderr,"ERROR: In function GAPS_APT_TensorValue: Order must be larger than -1\n");
		return 0;
	}
}

inline double GAPS_APT_CalGamma(double *pP,double *pMass)
{
	return sqrt(1.+(pP[0]*pP[0]+pP[1]*pP[1]+pP[2]*pP[2])/(pMass[0]*pMass[0]));
}

#include "APT_AllHeaders.h"

int GAPS_APT_Field_TokamakRE(double *pTensor,double *pSpaceTime4,int Order,Gaps_IO_InputsContainer *pInputs)
{	

//	printf("EMType=%d\n",pInputs->EMField_Type);
	int MaxOrder = -1;
	double TT=pSpaceTime4[0],XX=pSpaceTime4[1],YY=pSpaceTime4[2],ZZ=pSpaceTime4[3];
	double R0=pInputs->EMField_Tokamak_R0;
	double a =pInputs->EMField_Tokamak_a;
	double q =pInputs->EMField_Tokamak_q;
	double B0=pInputs->EMField_B0;
	double E0=pInputs->EMField_E0;

/*	int num_total_particles=pInputs->num_total_particles;
	double *pGatherData=(double *)calloc(num_total_particles*6,sizeof(double));
	GAPS_APT_GatherPtcInfo(pGatherData);
	int i;
	for(i=0;i<num_total_particles;i++)
	{
		double *test=pGatherData+i*6;
		printf("X = ");
		printf("%e,%e,%e;\t",test[0],test[1],test[2]);
		printf("P = ");
		printf("%e,%e,%e;",test[3],test[4],test[5]);
		printf("\n");
	}
*/
	
	double R2=XX*XX+YY*YY;
	double R=sqrt(R2);
	double B0overR2=B0/(R2*q);
	double r=sqrt((R-R0)*(R-R0)+ZZ*ZZ);
	double phi;
	double deltaB;
	double E0R0overR2=E0*R0/R2;

	if(-1 == Order)
	{
		if(pInputs->EMField_Cal_E)
		{
			//Calculate E and assign values to pTensor[0]~pTensor[2]
			pTensor[0]=-1*E0R0overR2*YY;
			pTensor[1]=E0R0overR2*XX;
			pTensor[2]=0.;
		}	
		if(pInputs->EMField_Cal_B)
		{
			//Calculate B and assign values to pTensor[3]~pTensor[5]
			pTensor[3]=B0overR2*(-q*R0*YY+XX*ZZ);
			pTensor[4]=B0overR2*(q*R0*XX+YY*ZZ);
			pTensor[5]=B0/q*(-1+R0/R);
		}		
	}

//	free(pGatherData);
	return 0; 
}


#include "APT_AllHeaders.h"
//Regular_Boris calc intgral point of x and v
double ** zhengjssetp();
static void zhengjsfreep();
int GAPS_APT_Pusher_Onehalf_Boris(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	/**Step 1: Get pointers of particle data**/
	double *pT = GAPS_APT_GetT1(pPtc);
	double *pX = GAPS_APT_GetX3(pPtc);
	
	double *pP = GAPS_APT_GetP3(pPtc);
	
	double *E = GAPS_APT_GetE3(pPtc);
	double *B = GAPS_APT_GetB3(pPtc);
	
	//check APT_Get function
	double *pMass = GAPS_APT_GetMass1(pPtc);
	double *pCharge= GAPS_APT_GetCharge1(pPtc);

	/**Step 2: Get parameters from pInputs**/
	double dT=pInputs->dT;
	
	/**Step 3: Update pT, pX, pP, E, B through an algorithm**/
	double F_ext[3],E_eff[3];
	//Update and Return E,B
	GAPS_APT_CalEB(E,B,pPtc,pInputs);
	
	//Update and Return sum of all extern forces
	GAPS_APT_MergeExtForce(F_ext,pPtc,pInputs);

	E_eff[0]=pCharge[0]*E[0]+F_ext[0];
	E_eff[1]=pCharge[0]*E[1]+F_ext[1];
	E_eff[2]=pCharge[0]*E[2]+F_ext[2];

	// Core of algorithm
	
	double ** omega =  zhengjssetp(3,3);
	double ** EYE = zhengjssetp(3,3);
	int i,j;
	for(i=0;i<3;i++){
		*(*(omega+i)+i) = 0;
	}

	*(*(omega+0)+1) = -B[2]*dT * pCharge[0] / pMass[0];
	*(*(omega+0)+2) =  B[1]*dT * pCharge[0] / pMass[0];
	*(*(omega+1)+0) =  B[2]*dT * pCharge[0] / pMass[0];
	*(*(omega+1)+2) =  B[0]*dT * pCharge[0] / pMass[0];
	*(*(omega+2)+0) = -B[1]*dT * pCharge[0] / pMass[0];
	*(*(omega+2)+1) = -B[0]*dT * pCharge[0] / pMass[0];
		
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			if(i==j)
				*(*(EYE+i)+j) =1;
			else
				
				*(*(EYE+i)+j) =0;
		}
	}
	
	double **IplusOM3; //declare
	IplusOM3 = zhengjssetp(3,3);//assign space
	M3plusM3(EYE, omega, IplusOM3);//calc tmpM3 = EYE + omega
	
	double **invIplusOM3; //declare
	invIplusOM3 = zhengjssetp(3,3);//assign space
	Invert_LU( invIplusOM3, IplusOM3,3);


	double **IminusOM3; //declare
	IminusOM3 = zhengjssetp(3,3);//assign space
	M3minusM3(EYE, omega, IminusOM3);//calc tmpM3 = EYE + omega
	
	double **R; //declare
	R = zhengjssetp(3,3);//assign space
	M3multiM3(invIplusOM3, IminusOM3, R);//calc tmpV3 = tmpM3 * pP;

	double tmp1V3[3];//declare
	double tmp2V3[3];//declare
	
	M3multiV3(R, pP,tmp1V3);//calc tmpV3 = tmpM3 * pP;
	M3multiV3(invIplusOM3, E ,tmp2V3);//calc tmpV3 = tmpM3 * pP;
	
	//printf("tmpV3%f,%f,%f",tmpV3[0],tmpV3[1],tmpV3[2]);
	//printf("\n");
	pP[0] = tmp1V3[0] + tmp2V3[0] * pCharge[0]/pMass[0] * dT;
	pP[1] = tmp1V3[1] + tmp2V3[1] * pCharge[0]/pMass[0] * dT; 
	pP[2] = tmp1V3[2] + tmp2V3[2] * pCharge[0]/pMass[0] * dT; 
	
	//printf("phP%f,%f,%f",phP[0],phP[1],phP[2]);
	
	//x_k+1
	pX[0]+=dT * pP[0]; 
	pX[1]+=dT * pP[1];
	pX[2]+=dT * pP[2];
	
	*pT +=dT;
	// End: Core of algorithm
	zhengjsfreep(omega,3);
	zhengjsfreep(EYE,3);
	zhengjsfreep(IplusOM3,3);
	zhengjsfreep(invIplusOM3,3);
	zhengjsfreep(IminusOM3,3);
	zhengjsfreep(R,3);
	return 0;
}
//C = A + B 
inline int M3plusM3( double **A, double **B, double **C){
	int i,j;
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			*( *(C+i)+ j) =  *( *(A+i)+ j) +  *( *(B+i)+ j);
		
	return 0;
}

//C = A - B
inline int M3minusM3(double **A, double **B, double **C){
	int i,j;
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			*( *(C+i)+ j) =  *( *(A+i)+ j) -  *( *(B+i)+ j);
		
	return 0;
}

//B = M * A
inline int M3multiV3 (double **M, double *A, double *B){

	*B	= (*(*M)) * (*A) +(*(*M +1)) * (*(A+1)) + (*(*M + 2)) * (*(A+2));
	*(B+1)	= (*(*(M+1))) * (*A) +(*(*(M+1) +1)) * (*(A+1)) + (*(*(M+1) + 2)) * (*(A+2));
	*(B+2)	= (*(*(M+2))) * (*A) +(*(*(M+2) +1)) * (*(A+1)) + (*(*(M+2) + 2)) * (*(A+2));
	return 0;
}

inline double ** zhengjssetp(int M, int N){
	double **A;
	A=(double **)malloc( M*sizeof(double *) ); //分配n行
	int i;
	for( i=0;i<M;i++ )
	{
		    A[i]=(double *)malloc( N*sizeof(double) ); //为每行分配m个数据空间
	}
	return A;
}
static inline void zhengjsfreep(double **A,int d){
	int i;
	for(i=0;i<d;i++){
	free( (void *)A[i] );
	}
	free((void *)A);
}

//M3 = M1 * M2
inline int M3multiM3 (double **M1, double **M2, double **M3){
	
	*(*(M3+0)+0)	= (*(*(M1+0)+0)) * (*(*(M2+0)+0)) +(*(*(M1+0)+1)) * (*(*(M2+1)+0)) + (*(*(M1+0)+2)) * (*(*(M2+2)+0));
	*(*(M3+0)+1)	= (*(*(M1+0)+0)) * (*(*(M2+0)+1)) +(*(*(M1+0)+1)) * (*(*(M2+1)+1)) + (*(*(M1+0)+2)) * (*(*(M2+2)+1));
	*(*(M3+0)+2)	= (*(*(M1+0)+0)) * (*(*(M2+0)+2)) +(*(*(M1+0)+1)) * (*(*(M2+1)+2)) + (*(*(M1+0)+2)) * (*(*(M2+2)+2));
	*(*(M3+1)+0)	= (*(*(M1+1)+0)) * (*(*(M2+0)+0)) +(*(*(M1+1)+1)) * (*(*(M2+1)+0)) + (*(*(M1+1)+2)) * (*(*(M2+2)+0));
	*(*(M3+1)+1)	= (*(*(M1+1)+0)) * (*(*(M2+0)+1)) +(*(*(M1+1)+1)) * (*(*(M2+1)+1)) + (*(*(M1+1)+2)) * (*(*(M2+2)+1));
	*(*(M3+1)+2)	= (*(*(M1+1)+0)) * (*(*(M2+0)+2)) +(*(*(M1+1)+1)) * (*(*(M2+1)+2)) + (*(*(M1+1)+2)) * (*(*(M2+2)+2));
	*(*(M3+2)+0)	= (*(*(M1+2)+0)) * (*(*(M2+0)+0)) +(*(*(M1+2)+1)) * (*(*(M2+1)+0)) + (*(*(M1+2)+2)) * (*(*(M2+2)+0));
	*(*(M3+2)+1)	= (*(*(M1+2)+0)) * (*(*(M2+0)+1)) +(*(*(M1+2)+1)) * (*(*(M2+1)+1)) + (*(*(M1+2)+2)) * (*(*(M2+2)+1));
	*(*(M3+2)+2)	= (*(*(M1+2)+0)) * (*(*(M2+0)+2)) +(*(*(M1+2)+1)) * (*(*(M2+1)+2)) + (*(*(M1+2)+2)) * (*(*(M2+2)+2));
return 0;
}


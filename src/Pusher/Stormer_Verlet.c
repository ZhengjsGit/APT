#include "APT_AllHeaders.h"
//Regular_Boris calc intgral point of x and v
static inline double ** zhengjssetp();
static inline void zhengjsfreep();
static inline int M3plusM3(double **C, double **A, double **B);

static inline int M3minusM3(double **C, double **A, double **B);


static inline int M3multiV3 (double **M, double *A, double *B);


int GAPS_APT_Pusher_Stormer_Verlet(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
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
	//Step 3: Update pT, pX, pP, E, B through an algorithm
	double F_ext[3],E_eff[3];
	//Update and Return E,B
	GAPS_APT_CalEB(E,B,pPtc,pInputs);
	//Update and Return sum of all extern forces
	GAPS_APT_MergeExtForce(F_ext,pPtc,pInputs);

	E_eff[0]=pCharge[0]*E[0]+F_ext[0];
	E_eff[1]=pCharge[0]*E[1]+F_ext[1];
	E_eff[2]=pCharge[0]*E[2]+F_ext[2];
	// Core of algorithm
	
	double dT_half=0.5*dT;	
	double K1[3], K2[3], L1[3], L2[3];
	double ** omega =  zhengjssetp(3,3);
	double ** MatB = zhengjssetp(3,3);
	double ** EYE = zhengjssetp(3,3);
	double pXplusHalf[4];
	double EB[6];
	
	int i,j;
	
	pXplusHalf[0] = 0;
	
	for(i=0;i<3;i++){
		pXplusHalf[i+1] = pX[i] + dT_half * pP[i]; 
	}
//calcuate B at pX + halfT *Vn
	(pPtc->FieldFunc)(EB,pXplusHalf,-1,pInputs);
	
	for(i=0;i<3;i++){
		E[i] = EB[i];
	}
	
	for(i=0;i<3;i++){
		B[i] = EB[i+3];
	}
//done

//init omega
	for(i=0;i<3;i++){
		*(*(omega+i)+i) = 0;
	}

	*(*(omega+0)+1) = -B[2]*dT_half;
	*(*(omega+0)+2) =  B[1]*dT_half;
	*(*(omega+1)+0) =  B[2]*dT_half;
	*(*(omega+1)+2) =  B[0]*dT_half;
	*(*(omega+2)+0) = -B[1]*dT_half;
	*(*(omega+2)+1) = -B[0]*dT_half;
//done
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
		*(*(MatB+i)+j) = *(*(omega+i)+j) / dT_half;
		}

	}
	
	
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			if(i==j)
				*(*(EYE+i)+j) =1;
			else
				
				*(*(EYE+i)+j) =0;
		}
	}
	
	for(i=0;i<3;i++){
		K1[i] = pP[i];
	}
	
	
	double tmpV3[3];//declare
	M3multiV3(MatB, pP,tmpV3);//calc tmpV3 = MatB * pP;
	
	for(i=0;i<3;i++){
		L1[i] = *pCharge / *pMass * ( - tmpV3[i] + E[i]);
	}

	for(i=0;i<3;i++){
		tmpV3[i] = pP[i] + L1[i] * dT_half;
	}
	double tmpV3_1[3];
	M3multiV3(MatB, tmpV3, tmpV3_1);//calc tmpV3_1 = MatB * tmpV3;
	for(i=0;i<3;i++){
		tmpV3[i] =  *pCharge / *pMass *( - tmpV3_1[i] + E[i]);
	}
	
	double **tmpM3; //declare
	tmpM3 = zhengjssetp(3,3);//assign space
	M3plusM3(tmpM3, EYE, omega);//calc tmpM3 = EYE - omega
	
	//printf("tmpV3%f,%f,%f",tmpV3[0],tmpV3[1],tmpV3[2]);
	//printf("\n");
	double **tmpM3_1;
	tmpM3_1 = zhengjssetp(3,3);
	Invert_LU( tmpM3_1, tmpM3,3);
	M3multiV3(tmpM3_1,tmpV3,L2);
	for(i=0;i<3;i++){
		K2[i] =  pP[i] + dT_half * (L1[i] + L2[i]);
	}
	for(i=0;i<3;i++){
		pX[i] = pX[i] + dT_half * (K1[i] + K2[i]);
	}
	for(i=0;i<3;i++){
		pP[i] = pP[i] + dT_half * (L1[i] + L2[i]);
	}
	*pT +=dT;
	// End: Core of algorithm
	zhengjsfreep(omega,3);
	zhengjsfreep(MatB,3);
	zhengjsfreep(EYE,3);
	zhengjsfreep(tmpM3,3);
	zhengjsfreep(tmpM3_1,3);
	return 0;
}
static inline int M3plusM3(double **C, double **A, double **B){
	int i,j;
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			*( *(C+i)+ j) =  *( *(A+i)+ j) +  *( *(B+i)+ j);
		
	return 0;
}

static inline int M3minusM3(double **C, double **A, double **B){
	int i,j;
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			*( *(C+i)+ j) =  *( *(A+i)+ j) -  *( *(B+i)+ j);
		
	return 0;
}
static inline int M3multiV3 (double **M, double *A, double *B){

	*B	= (*(*M)) * (*A) +(*(*M +1)) * (*(A+1)) + (*(*M + 2)) * (*(A+2));
	*(B+1)	= (*(*(M+1))) * (*A) +(*(*(M+1) +1)) * (*(A+1)) + (*(*(M+1) + 2)) * (*(A+2));
	*(B+2)	= (*(*(M+2))) * (*A) +(*(*(M+2) +1)) * (*(A+1)) + (*(*(M+2) + 2)) * (*(A+2));
	return 0;
}

static inline double ** zhengjssetp(int M, int N){
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
	free( A[i] );
	}
	free(A);
}

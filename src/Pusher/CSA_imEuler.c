#include "APT_AllHeaders.h"
//void RFfunc_CanonSymp3D(double *X_iter,double *X_last,double *A,double (*DelA)[3],double *DelPHI, double *F,double dT);
//void RFJacobi_CanonSymp3D(double *X_iter,double *X_last,double *A,double (*DelA)[3],double *DelPHI, double **J,double dT);
int GAPS_APT_CSimEulerRootFindPushP_nonRe(double *P,  double *pDt,double *pTolerance,double *pCharge,double *pMass, Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
int GAPS_APT_CSimEulerRootFindFlow_nonRe(double *F,   double *P_iter,double *P_last, double *DelPHI_mid, double *A_mid , double (*DelA_mid)[3], double *pDt,double *pCharge,double *pMass);
int GAPS_APT_CSimEulerRootFindJacobi_nonRe(double **J,   double *P_iter,double *P_last, double *A_mid, double (*DelA_mid)[3], double (*DDelA_mid)[3][3],double *pDt,double *pCharge,double *pMass);
void Trans4to3(double (*DDelA)[3][3], double (*DelA)[3], double * DelPHI, double * A,   double * DDelA4, double * DelA4,double * pA4);
int GAPS_APT_Pusher_CSA_imEuler(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double Tolerance=pInputs->Pusher_RootFindingTol;
	double dT=pInputs->dT;
	double *pX4 = GAPS_APT_GetX4(pPtc);
	double *pP4 = GAPS_APT_GetP4(pPtc);
	double *pA4 = GAPS_APT_GetA4(pPtc);//find out this step
	//double *pGamma=GAPS_APT_GetGamma1(pPtc);//get gamma

	double *pMass = GAPS_APT_GetMass1(pPtc);
	double *pCharge= GAPS_APT_GetCharge1(pPtc);

	double CanP4[4];

	//	double DelA4[16];
	//Key point, write expression of function A4 and Del A4
	//Field function just what we set in bash script. maybe this one has parameters
	//see? using order

	(pPtc->FieldFunc)(pA4,pX4,1,pInputs);
	//	(pPtc->FieldFunc)(DelA4,pX4,2,pInputs);
	{int i;for(i=0;i<4;i++){CanP4[i]=pP4[i]+pCharge[0]*pA4[i];}}

	//	double DDelA3[3][3][3];
	//	double DelA3[3][3];
	//  double A3[3],DelPHI[3];
	double CanxP6[6]={pX4[1],pX4[2],pX4[3],CanP4[1],CanP4[2],CanP4[3]};


	GAPS_APT_CSimEulerRootFindPushxP_nonRe(CanxP6,&dT,&Tolerance,pCharge,pMass, pPtc, pInputs);

	pX4[1]=CanxP6[0];
	pX4[2]=CanxP6[1];
	pX4[3]=CanxP6[2];
	pX4[0]+=dT;

	pP4[1]=CanxP6[3];
	pP4[2]=CanxP6[4];
	pP4[3]=CanxP6[5];
	pP4[0]=1;
	//new X
	(pPtc->FieldFunc)(pA4,pX4,1,pInputs);

	pP4[1]=CanxP6[3]-pCharge[0]*pA4[1];
	pP4[2]=CanxP6[4]-pCharge[0]*pA4[2];
	pP4[3]=CanxP6[5]-pCharge[0]*pA4[3];

	//	EraseMatrixSpace(J);
	return 0;
}
//Newton iteration main progress
//LocalA DelPHI DelA DDelA should be calculated at half way of X
//missing DDelA
//xp is input arg, initial x and P and we keep iterating on it until it gets out
int GAPS_APT_CSimEulerRootFindPushxP_nonRe(double * xP,double *pDt,double *pTolerance,double *pCharge,double *pMass,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	//xP_last is n_th x and P
	double xP_last[6]={xP[0],xP[1],xP[2],xP[3],xP[4],xP[5]};
	//iteration container for x and P
	double xP_iter[6]={xP[0],xP[1],xP[2],xP[3],xP[4],xP[5]};
	double F[6];
	double Jacobi_1D[36];
	double **Jacobi=(double **)malloc(sizeof(double *)*6);
	Jacobi[0]=&Jacobi_1D[0];
	Jacobi[1]=&Jacobi_1D[6];
	Jacobi[2]=&Jacobi_1D[12];
	Jacobi[3]=&Jacobi_1D[18];
	Jacobi[4]=&Jacobi_1D[24];
	Jacobi[5]=&Jacobi_1D[30];
	//alloc memory DONE

	double tol,tolF;
	long times_iter=0;
	do
	{
		int i,j;
		double delta_xP[6];
		double fold=0.;//0.5*F*F
		double g[6];//Del(0.5*F*F)

		//since we cannt using a static A and DDel A Del A ... (they involve xP n+1)
		//We have to calculate it inside iteration.


		//prepare Field components, declare and alloc memory
		double DDelA4_mid[64];
		double DelA4_mid[16];
		double pA4_mid[4];
		//avg position
		double pXtmp[4];
		pXtmp[0] = 0;
		
		for(i=0;i<3;i++){
			pXtmp[i+1] = ( xP_last[i] + xP_iter[i] )/2.0 ;
		}

		(pPtc->FieldFunc)(DDelA4_mid, pXtmp, 3, pInputs);
		(pPtc->FieldFunc)(DelA4_mid, pXtmp, 2, pInputs);
		(pPtc->FieldFunc)(pA4_mid,pXtmp,1,pInputs);

		double DDelA_mid[3][3][3];
		double DelA_mid[3][3];
		double A_mid[3],DelPHI_mid[3];

		Trans4to3(DDelA_mid, DelA_mid, DelPHI_mid, A_mid, DDelA4_mid, DelA4_mid, pA4_mid);

		//FLOW and JACOBI wants A DelA DelPHI and DDelA at half point of X, so PushxP should satisfied this
		GAPS_APT_CSimEulerRootFindFlow_nonRe(F,xP_iter,xP_last,DelPHI_mid,A_mid,DelA_mid,pDt,pCharge,pMass);

		tolF=sqrt(F[0]*F[0]+F[1]*F[1]+F[2]*F[2]+F[3]*F[3]+F[4]*F[4]+F[5]*F[5]);
		//MISSING DDelA
		GAPS_APT_CSimEulerRootFindJacobi_nonRe(Jacobi,xP_iter,xP_last,A_mid,DelA_mid,DDelA_mid,pDt,pCharge,pMass);
		if(500 < times_iter){
	    for(i=0;i<6;i++){
		   for(j=0;j<6;j++){
		       printf("%.2f, ",Jacobi[i][j]);
		   }
		   printf("\n");
		}
		printf("iteration too many, abort");
		exit(1);
		}
        //printf("----------------------------------------------\n");
		double neg_F[6]={-1*F[0],-1*F[1],-1*F[2],-1*F[3],-1*F[4],-1*F[5]};
		LinearSolver_LU(delta_xP,Jacobi,neg_F,6);
		double xP_new[6];
		double fnew;
		tol=0;
		for(i=0;i<6;i++)
		{
			xP_iter[i]+=delta_xP[i];
			tol+=pow(delta_xP[i],2);
		}
		tol=sqrt(tol);	
		times_iter++;
		//printf("%d\n",times_iter);
	}while(tol>=(*pTolerance) && tolF>=(*pTolerance));
	//Pusha x P
	int i;
	for(i=0;i<6;i++)
	{
		xP[i] =	xP_iter[i];
	}

	free(Jacobi);
	return 0;
}
//What is Flow? Flow is [deltaP - dT * (partial H over partial q)]
//we have to calculate A delA delPhi inside iteration
//actually we using DelA DelPhi in our parameter p_N+1 + p_N divide 2
int GAPS_APT_CSimEulerRootFindFlow_nonRe(double *F,double *xP_iter,double *xP_last,double * DelPHI_mid, double * A_mid,double (*DelA_mid)[3], double *pDt,double *pCharge,double *pMass)
{
	double Charge=(*pCharge);
	double Mass=(*pMass);
	double ChargeMass=Charge/Mass;
	double dT=(*pDt);
	int i;
	//for(i=1;i<4;i++){pPtmp[i] = (xP_last[i-1] + xP_iter[i-1] )/2.0}
	for(i=0;i<3;i++){
		//correct1 
		F[i] = xP_iter[i] - xP_last[i] - dT *( (xP_iter[i+3] + xP_last[i+3])/2.0 - Charge * A_mid[i] );
	}
	for(i=0;i<3;i++)
	{
		int j;
		double sumP=0;
		double avgP=0;
		for(j=0;j<3;j++)
		{	
			avgP =( xP_iter[j+3] + xP_last[j+3] )/ 2.0;
			sumP+= - (avgP - Charge * A_mid[j]) * DelA_mid[j][i];
		}
		F[i+3] = (xP_iter[i+3]-xP_last[i+3] + dT*(ChargeMass*sumP + Charge * DelPHI_mid[i]));			
	}
	return 0;
}

int GAPS_APT_CSimEulerRootFindJacobi_nonRe(double **J,double *P_iter,double *P_last, double *A_mid, double (*DelA_mid)[3] ,double (*DDelA_mid)[3][3], double *pDt,double *pCharge,double *pMass)
{
	double Charge=(*pCharge);
	double Mass=(*pMass);
	double ChargeMass=Charge/Mass;
	double dT=(*pDt);


//correct2 
	int i,j,k;
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			J[i][j] = Delta(i,j) + dT * DelA_mid[i][j] /2.0 ; 
		}
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			J[i][j+3] = - Delta(i,j) * dT /2.0; 
		}
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			double avgP;
			double sumP = 0;
			for(k=0;k<3;k++){
				avgP =( P_iter[k] + P_last[k] )/ 2.0;
				sumP+= -(avgP - A_mid[k])*DDelA_mid[k][j][j] / 2.0 + DelA_mid[k][j] * DelA_mid[k][j] /2.0; 
			}
			J[i+3][j] = sumP * dT; 
		}
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++)
		{
			J[i+3][j+3]=Delta(i,j)+dT*ChargeMass*(-DelA_mid[j][i])/2.0;
		}
	}
	return 0;
}

//IN DDelA4 DelA4 A4
//OUT PHI A3 DelA3[][] DDelA
//all pTensor(dim 4) is one dimensional
void Trans4to3(double (*DDelA3)[3][3], double (*DelA3)[3], double * DelPHI,double * A3, double * DDelA4, double * DelA4,double * pA4){

	//A3 and PHI are come from A4
	int i,j,k;
	for(i=0;i<3;i++)
	{
		A3[i]=pA4[i+1];
		int idx1[2]={0,i+1};
		DelPHI[i]=GAPS_APT_TensorValue(DelA4,idx1,4,2);
		for(j=0;j<3;j++)
		{
			int idx2[2]={i+1,j+1};
			DelA3[i][j]=GAPS_APT_TensorValue(DelA4,idx2,4,2);

			for(k=0;k<3;k++){
				int idx3[3]={i+1,j+1,k+1};
				DDelA3[i][j][k]=GAPS_APT_TensorValue(DDelA4,idx3,4,2);
			}
		}
	}
}

#include "APT_AllHeaders.h"
//void RFfunc_CanonSymp3D(double *X_iter,double *X_last,double *A,double (*DelA)[3],double *DelPHI, double *F,double dT);
//void RFJacobi_CanonSymp3D(double *X_iter,double *X_last,double *A,double (*DelA)[3],double *DelPHI, double **J,double dT);
int GAPS_APT_SVRootFindCoeffience(double *pK1, double *K2, double * L1, double * pL2, double *xP, double *pDt,double *pTolerance,double *pCharge,double *pMass, Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
int GAPS_APT_SVRootFindFlow_K1(double *F, double *K1, double *xP, double *A, double *pDt,double *pCharge,double *pMass);
int GAPS_APT_SVRootFindJacobi_K1(double **J, double *xP, double (*DelA)[3], double *pDt,double *pCharge,double *pMass);
int GAPS_APT_SVRootFindFlow_L2(double *F, double *L1, double * L2, double *xP, double *DelPHI, double *A, double (*DelA)[3], double *pDt,double *pCharge,double *pMass);
int GAPS_APT_SVRootFindJacobi_L2(double **J, double *xP, double (*DelA)[3], double *pDt,double *pCharge,double *pMass);


int GAPS_APT_Pusher_Stormer_Verlet(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
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
	//we have to compmise

	(pPtc->FieldFunc)(pA4,pX4,1,pInputs);
	//	(pPtc->FieldFunc)(DelA4,pX4,2,pInputs);

	//	double DDelA3[3][3][3];
	//	double DelA3[3][3];
	//  double A3[3],DelPHI[3];
	double CanxP6[6]={pX4[1],pX4[2],pX4[3],CanP4[1],CanP4[2],CanP4[3]};

	double *pK1= GAPS_APT_GetCanP3(pPtc);
	double K2[3];
	double L1[3];
	double *pL2= GAPS_APT_GetAclr3(pPtc);
	/*
	{int i;for(i=0;i<3;i++){printf("K1_%d=%.4f----",i,pK1[i]);}}
	printf("\n");
	{int i;for(i=0;i<3;i++){printf("K2_%d=%.4f----",i,K2[i]);}}
	printf("\n");
	{int i;for(i=0;i<3;i++){printf("L1_%d=%.4f----",i,L1[i]);}}
	printf("\n");
	{int i;for(i=0;i<3;i++){printf("L2_%d=%.4f----",i,pL2[i]);}}
	printf("\n####################################################\n");
	*/
	GAPS_APT_SVRootFindCoeffience(pK1, K2, L1, pL2, CanxP6, & dT, &Tolerance, pCharge, pMass, pPtc, pInputs);

	//{int i;for(i=0;i<3;i++){printf("L1_%d=%.4f\n",i,L1[i]);}}
	pP4[0]=1;
	pX4[0]+=dT;
	int i;
	for(i=0;i<3;i++)
	{
		pX4[i+1]= pX4[i+1] + dT / 2.0 * (pK1[i] + K2[i]) ;
		pP4[i+1]= pP4[i+1] + dT / 2.0 * (L1[i] + pL2[i]);
	}
	//new X
	(pPtc->FieldFunc)(pA4,pX4,1,pInputs);
	pP4[1]=pP4[1]-pCharge[0]*pA4[1];
	pP4[2]=pP4[2]-pCharge[0]*pA4[2];
	pP4[3]=pP4[3]-pCharge[0]*pA4[3];
	//	EraseMatrixSpace(J);
	return 0;
}
//Newton iteration main progress
//LocalA DelPHI DelA DDelA should be calculated at half way of X
//missing DDelA
//xp is input arg, initial x and P and we keep iterating on it until it gets out
int GAPS_APT_SVRootFindCoeffience(double *pK1, double *K2, double * L1, double *pL2, double *xP, double *pDt,double *pTolerance,double *pCharge,double *pMass, Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{

	//Global vars
	double DelA4[16];
	double pA4[4];
	double pXtmp[4];

	double DelA[3][3];
	double A[3],DelPHI[3];
	double Charge=(*pCharge);
	double Mass=(*pMass);
	double ChargeMass=Charge/Mass;
	double dT=(*pDt);

	int i,j,k;
	double K1[3] = {pK1[0],pK1[1],pK1[2]};
	double L2[3] = {pL2[0],pL2[1],pL2[2]};
	


	//Part 1 K1
	{
		//iteration container for x and P
		double F[3];
		double Jacobi_1D[9];
		double **Jacobi=(double **)malloc(sizeof(double *)*3);
		Jacobi[0]=&Jacobi_1D[0];
		Jacobi[1]=&Jacobi_1D[3];
		Jacobi[2]=&Jacobi_1D[6];
		//alloc memory DONE

		double tol,tolF;
		long times_iter=0;

		do
		{
			double delta_K1[3];

			//prepare Field components, declare and alloc memory
			pXtmp[0] = 0;
			for(i=0;i<3;i++){
				pXtmp[i+1] = xP[i] + dT /2.0 * K1[i];
			}

			(pPtc->FieldFunc)(DelA4, pXtmp, 2, pInputs);
			(pPtc->FieldFunc)(pA4,pXtmp,1,pInputs);

			//{int i;for(i=0;i<16;i++){printf("DelA_%d=%.2f\n",i,DelA4[i]);}}
			for(i=0;i<3;i++)
			{
				A[i]=pA4[i+1];
				for(j=0;j<3;j++)
				{
					int idx2[2]={i+1,j+1};
					DelA[i][j]=GAPS_APT_TensorValue(DelA4,idx2,4,2);
				}
			}
				//FLOW and JACOBI wants A DelA DelPHI and DDelA at half point of X, so PushxP should satisfied this
			GAPS_APT_SVRootFindFlow_K1(F, K1, xP, A, pDt, pCharge, pMass);

			tolF=sqrt(F[0]*F[0]+F[1]*F[1]+F[2]*F[2]);
			//MISSING DDelA
			GAPS_APT_SVRootFindJacobi_K1(Jacobi, xP, DelA, pDt, pCharge, pMass);
			if(500 < times_iter){
				for(i=0;i<3;i++){
					for(j=0;j<3;j++){
						printf("%.2f, ",Jacobi[i][j]);
					}
					printf("\n");
				}
				printf("iteration K1 too many, abort");
				exit(1);
			}
			//printf("----------------------------------------------\n");
			double neg_F[3]={-1*F[0],-1*F[1],-1*F[2]};
			LinearSolver_LU(delta_K1,Jacobi,neg_F,3);
			tol=0;
			for(i=0;i<3;i++)
			{
				K1[i]+=delta_K1[i];
				tol+=pow(delta_K1[i],2);
			}
			tol=sqrt(tol);	
			times_iter++;
			//printf("%d\n",times_iter);
		}while(tol>=(*pTolerance) && tolF>=(*pTolerance));
		free(Jacobi);
		for(i=0;i<3;i++){pK1[i] = K1[i];}
	}

//Update A value using new K
		pXtmp[0] = 0;
		for(i=0;i<3;i++){
			pXtmp[i+1] = xP[i] + dT /2.0 * K1[i];
		}
		(pPtc->FieldFunc)(DelA4, pXtmp, 2, pInputs);
		(pPtc->FieldFunc)(pA4,pXtmp,1,pInputs);
//reshape
		//{int i;for(i=0;i<16;i++){printf("DelA_%d=%.2f\n",i,DelA4[i]);}}
		//Trans4to3_zjs(DDelA, DelA, DelPHI, A, DDelA4, DelA4, pA4);
		for(i=0;i<3;i++)
		{
			A[i]=pA4[i+1];
			int idx1[2]={0,i+1};
			DelPHI[i]=GAPS_APT_TensorValue(DelA4,idx1,4,2);
			for(j=0;j<3;j++)
			{
				int idx2[2]={i+1,j+1};
				DelA[i][j]=GAPS_APT_TensorValue(DelA4,idx2,4,2);
			}
		}
//DONE
	//Part 2 L1
	{
		for(i=0;i<3;i++)
		{
			double sumP = 0;
			for(j=0;j<3;j++){
				sumP += - 1.0/Mass * (xP[j+3] - Charge * A[j]) * Charge * DelA[j][i];
			}
			L1[i] = - sumP + Charge * DelPHI[i] ;
		}

	}
	//Part 3 L2
	{
		//iteration container for x and P
		double F[3];
		double Jacobi_1D[9];
		double **Jacobi=(double **)malloc(sizeof(double *)*3);
		Jacobi[0]=&Jacobi_1D[0];
		Jacobi[1]=&Jacobi_1D[3];
		Jacobi[2]=&Jacobi_1D[6];
		//alloc memory DONE

		double tol,tolF;
		long times_iter=0;

		do
		{
			int i,j;
			double delta_L2[3];

			GAPS_APT_SVRootFindFlow_L2(F, L1, L2, xP, DelPHI, A, DelA, pDt, pCharge, pMass);

			tolF=sqrt(F[0]*F[0]+F[1]*F[1]+F[2]*F[2]);
			//MISSING DDelA
			GAPS_APT_SVRootFindJacobi_L2(Jacobi, xP, DelA, pDt, pCharge, pMass);
			if(500 < times_iter){
				for(i=0;i<3;i++){
					for(j=0;j<3;j++){
						printf("%.2f, ",Jacobi[i][j]);
					}
					printf("\n");
				}
				printf("iteration L2 too many, abort");
				exit(1);
			}
			//printf("----------------------------------------------\n");
			double neg_F[3]={-1*F[0],-1*F[1],-1*F[2]};
			LinearSolver_LU(delta_L2,Jacobi,neg_F,3);
			tol=0;
			for(i=0;i<3;i++)
			{
				L2[i]+=delta_L2[i];
				tol+=pow(delta_L2[i],2);
			}
			tol=sqrt(tol);	
			times_iter++;
			//printf("%d\n",times_iter);
		}while(tol>=(*pTolerance) && tolF>=(*pTolerance));
		free(Jacobi);
		for(i=0;i<3;i++){pL2[i] = L2[i];}
	}
	//Part 4 K2
	{
		for(i=0;i<3;i++)
		{
			K2[i] = 1.0 / Mass *(xP[i+3] + (L1[i] + L2[i]) *dT /2.0 - Charge * A[i]);			
		}
	}
	return 0;
}

//What is Flow? Flow is [deltaP - dT * (partial H over partial q)]
//we have to calculate A delA delPhi inside iteration
//A relate q + k/2
int GAPS_APT_SVRootFindFlow_K1(double *F, double *K1, double *xP, double *A,  double *pDt,double *pCharge,double *pMass)
{
	double Charge=(*pCharge);
	double Mass=(*pMass);
	double ChargeMass=Charge/Mass;
	double dT=(*pDt);
	int i;
	for(i=0;i<3;i++)
	{
		F[i] =	K1[i] - 1.0 / Mass * (xP[i+3] - Charge * A[i]);
	}

	return 0;
}

int GAPS_APT_SVRootFindJacobi_K1(double **J, double *xP,  double (*DelA)[3], double *pDt,double *pCharge,double *pMass)
{
	double Charge=(*pCharge);
	double Mass=(*pMass);
	double ChargeMass=Charge/Mass;
	double dT=(*pDt);

	int i,j,k;
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			J[i][j] = Delta(i,j) + 1.0/Mass*(dT / 2.0 * Charge * DelA[i][j] );
		}
	}

	return 0;
}
//A related at +h/2
int GAPS_APT_SVRootFindFlow_L2(double *F, double *L1, double *L2, double *xP, double *DelPHI, double *A, double (*DelA)[3], double *pDt,double *pCharge,double *pMass)
{
	double Charge=(*pCharge);
	double Mass=(*pMass);
	double ChargeMass=Charge/Mass;
	double dT=(*pDt);
	int i;
	for(i=0;i<3;i++)
	{
		int j;
		double sumP = 0;
		for(j=0;j<3;j++){
			sumP += - 1.0/Mass * (xP[j+3] + dT / 2.0 *(L1[i] + L2[i]) - Charge * A[j]) * Charge * DelA[j][i];
		}
		F[i] =L2[i] + ( sumP + DelPHI[i]) ;
	}
	return 0;
}
//A at qn + l1 L2 
int GAPS_APT_SVRootFindJacobi_L2(double **J, double *xP, double (*DelA)[3], double *pDt,double *pCharge,double *pMass)
{
	double Charge=(*pCharge);
	double Mass=(*pMass);
	double ChargeMass=Charge/Mass;
	double dT=(*pDt);

	int i,j;
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			J[i][j] = Delta(i,j) - 1.0/Mass*(dT / 2.0 * Charge * DelA[j][i] );
		}
	}

	return 0;
}

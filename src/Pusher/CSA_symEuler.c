#include "APT_AllHeaders.h"
//void RFfunc_CanonSymp3D(double *X_iter,double *X_last,double *A,double (*DelA)[3],double *DelPHI, double *F,double dT);
//void RFJacobi_CanonSymp3D(double *X_iter,double *X_last,double *A,double (*DelA)[3],double *DelPHI, double **J,double dT);
int GAPS_APT_CSEulerRootFindPushP_nonRe(double *P,double *LocalA,double *DelPHI,double (*LocalDelA)[3],double *pDt,double *pTolerance,double *pCharge,double *pMass);
int GAPS_APT_CSEulerRootFindFlow_nonRe(double *F,double *P_iter,double *P_last,double *A,double (*DelA)[3],double *DelPHI,double *pDt,double *pCharge,double *pMass);
int GAPS_APT_CSEulerRootFindJacobi_nonRe(double **J,double *P_iter,double *P_last,double *A,double (*DelA)[3],double *pDt,double *pCharge,double *pMass);

int GAPS_APT_Pusher_CSA_SymEuler(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
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
	double DelA4[16];
	//Key point, write expression of function A4 and Del A4
	//Field function just what we set in bash script. maybe this one has parameters
	//see? using order
	(pPtc->FieldFunc)(pA4,pX4,1,pInputs);
	(pPtc->FieldFunc)(DelA4,pX4,2,pInputs);
	{int i;for(i=0;i<4;i++){CanP4[i]=pP4[i]+pCharge[0]*pA4[i];}}

	double DelA3[3][3];
	double PHI=pA4[0];
	double A3[3],DelPHI[3];
//A3 and PHI are come from A4
	{//Transform 4Vectors to 3D vectors
		int i,j;
		for(i=0;i<3;i++)
		{
			A3[i]=pA4[i+1];
			int idx1[2]={0,i+1};
			DelPHI[i]=GAPS_APT_TensorValue(DelA4,idx1,4,2);
			for(j=0;j<3;j++)
			{
				int idx2[2]={i+1,j+1};
				DelA3[i][j]=GAPS_APT_TensorValue(DelA4,idx2,4,2);
			}
		}
	}//Transform 4Vectors to 3D vectors--end

	double CanP3[3]={CanP4[1],CanP4[2],CanP4[3]};
	GAPS_APT_CSEulerRootFindPushP_nonRe(CanP3,A3,DelPHI,DelA3, &dT,&Tolerance,pCharge,pMass);

	//Push X
	int i;
	double Dx[3];
//	double normDx=0;
	//double gamma_Inv=1./sqrt((pow(CanP3[0]-pCharge[0]*A3[0],2)+pow(CanP3[1]-pCharge[0]*A3[1],2)+pow(CanP3[2]-pCharge[0]*A3[2],2))/(pMass[0]*pMass[0])+1);
	for(i=0;i<3;i++)
	{
		Dx[i]=(dT)*(CanP3[i]-pCharge[0]*A3[i]);
		pX4[i+1]+=Dx[i];
//		normDx+=pow(Dx[i],2);
	}
	pX4[0]+=dT;
	
	(pPtc->FieldFunc)(pA4,pX4,1,pInputs);

	pP4[1]=CanP3[0]-pCharge[0]*pA4[1];
	pP4[2]=CanP3[1]-pCharge[0]*pA4[2];
	pP4[3]=CanP3[2]-pCharge[0]*pA4[3];
	pP4[0]=GAPS_APT_CalGamma(pP4+1,pMass);

//	EraseMatrixSpace(J);
	return 0;
}
//Newton iteration main progress
int GAPS_APT_CSEulerRootFindPushP_nonRe(double *P,double *LocalA,double *DelPHI,double (*LocalDelA)[3],double *pDt,double *pTolerance,double *pCharge,double *pMass)
{
	double P_last[3]={P[0],P[1],P[2]};
	double P_iter[3]={P[0],P[1],P[2]};
	double F[3];
	double Jacobi_1D[9];
	double **Jacobi=(double **)malloc(sizeof(double *)*3);
	Jacobi[0]=&Jacobi_1D[0];
	Jacobi[1]=&Jacobi_1D[3];
	Jacobi[2]=&Jacobi_1D[6];

	double tol,tolF;

	long times_iter=0;
	do
	{
		int i,j;
		double delta_P[3];
		double fold=0.;//0.5*F*F
		double g[3];//Del(0.5*F*F)

		GAPS_APT_CSEulerRootFindFlow_nonRe(F,P_iter,P_last,LocalA,LocalDelA,DelPHI,pDt,pCharge,pMass);

		tolF=sqrt(F[0]*F[0]+F[1]*F[1]+F[2]*F[2]);
		GAPS_APT_CSEulerRootFindJacobi_nonRe(Jacobi,P_iter,P_last,LocalA,LocalDelA,pDt,pCharge,pMass);
		double neg_F[3]={-1*F[0],-1*F[1],-1*F[2]};
		LinearSolver_LU(delta_P,Jacobi,neg_F,3);
		double P_new[3];
		double fnew;
		tol=0;
		for(i=0;i<3;i++)
		{
			P_iter[i]+=delta_P[i];
			tol+=pow(delta_P[i],2);
		}
		tol=sqrt(tol);	
		times_iter++;
	}while(tol>=(*pTolerance) && tolF>=(*pTolerance));
	//Push P
	P[0]=P_iter[0];
	P[1]=P_iter[1];
	P[2]=P_iter[2];
	free(Jacobi);
	return 0;
}
//What is Flow? Flow is [deltaP - dT * (partial H over partial q)]
int GAPS_APT_CSEulerRootFindFlow_nonRe(double *F,double *P_iter,double *P_last,double *A,double (*DelA)[3],double *DelPHI,double *pDt,double *pCharge,double *pMass)
{
	double Charge=(*pCharge);
	double Mass=(*pMass);
	double ChargeMass=Charge/Mass;
	double dT=(*pDt);
	//double gamma=sqrt((pow(P_iter[0]-Charge*A[0],2)+pow(P_iter[1]-Charge*A[1],2)+pow(P_iter[2]-Charge*A[2],2))/(Mass*Mass)+1);
	//double gamma_rev=1./gamma;

	int i;
	for(i=0;i<3;i++)
	{
		int j;
		double sumP=0;
		for(j=0;j<3;j++)
		{
			sumP+=(P_iter[j]-Charge*A[j])*DelA[j][i];
		}
		//F[i] = 	(P_iter[i]-P_last[i]-dT*(ChargeMass*sumP*gamma_rev-Charge*DelPHI[i]));			
		F[i] = 	(P_iter[i]-P_last[i]-dT*(ChargeMass*sumP - Charge*DelPHI[i]));			
	}
	return 0;
}

int GAPS_APT_CSEulerRootFindJacobi_nonRe(double **J,double *P_iter,double *P_last,double *A,double (*DelA)[3],double *pDt,double *pCharge,double *pMass)
{
	double Charge=(*pCharge);
	double Mass=(*pMass);
	double ChargeMass=Charge/Mass;
	double dT=(*pDt);
	//double gamma=sqrt((pow(P_iter[0]-Charge*A[0],2)+pow(P_iter[1]-Charge*A[1],2)+pow(P_iter[2]-Charge*A[2],2))/(Mass*Mass)+1);
	//double gamma_rev=1./gamma;

	int i,j,k;
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			double sumP=0;
			for(k=0;k<3;k++)
			{
				//sumP+=(P_iter[k]-Charge*A[k])*(P_iter[j]-Charge*A[j])*DelA[k][i]*pow(gamma_rev,3);
				//sumP+=(P_iter[k]-Charge*A[k])*(P_iter[j]-Charge*A[j])*DelA[k][i];
			}
			J[i][j]=Delta(i,j)+dT*ChargeMass*(sumP-DelA[j][i]);		
		}
	}
	return 0;
}

/*void RFfunc_CanonSymp3D(double *X_iter,double *X_last,double *A,double (*DelA)[3],double *DelPHI, double *F,double dT)
{
	double xlast[3]={X_last[0],X_last[1],X_last[2]};	
	double plast[3]={X_last[3],X_last[4],X_last[5]};	

	double xn[3]={X_iter[0],X_iter[1],X_iter[2]};	
	double pn[3]={X_iter[3],X_iter[4],X_iter[5]};	
	double gamma=sqrt(pow(pn[0]-A[0],2)+pow(pn[1]-A[1],2)+pow(pn[2]-A[2],2)+1);

	int i;
	for(i=0;i<6;i++)
	{
		if(i<3)
		{
			F[i] = 	-1*(xn[i]-xlast[i]-1*dT*(pn[i]-A[i])/gamma);			
		}
		else
		{
			int n=i-3;
			int j;
			double sumP=0;
			for(j=0;j<3;j++)
			{
				sumP+=(pn[j]-A[j])/gamma*DelA[j][n];
			}

			F[i] = 	-1*(pn[n]-plast[n]-1*dT*(sumP-DelPHI[n]));			
		}
	}
}

void RFJacobi_CanonSymp3D(double *X_iter,double *X_last,double *A,double (*DelA)[3],double *DelPHI, double **J,double dT)
{
	double xlast[3]={X_last[0],X_last[1],X_last[2]};	
	double plast[3]={X_last[3],X_last[4],X_last[5]};	

	double xn[3]={X_iter[0],X_iter[1],X_iter[2]};	
	double pn[3]={X_iter[3],X_iter[4],X_iter[5]};	
	
	double gamma=sqrt(pow(pn[0]-A[0],2)+pow(pn[1]-A[1],2)+pow(pn[2]-A[2],2)+1);

	int i,j,k;
	for(i=0;i<6;i++)
	{
		for(j=0;j<6;j++)
		{
			if(i<3 && j<3)
			{
				J[i][j]=Delta(i,j);		
			}
			if(i<3 && j>=3)
			{
				int m=i,n=j-3;
				J[i][j]=-1*dT*(Delta(m,n)/gamma-(pn[m]-A[m])*(pn[n]-A[n])/pow(gamma,3));		
			}
			if(i>=3 && j<3)
			{
				J[i][j]=0;		
			}
			if(i>=3 && j>=3)
			{
				double sumP=0;
				int m=i-3,n=j-3;
				for(k=0;k<3;k++)
				{
					sumP+=(Delta(k,n)/gamma-(pn[k]-A[k])*(pn[n]-A[n])/pow(gamma,3))*DelA[k][m];
				}
				J[i][j]=Delta(m,n)-dT*sumP;		
			}
				
		}
	}
}
*/

#include "APT_AllHeaders.h"
int GAPS_APT_Pusher_LCCSA_IMP_RootFindPush(double *XP,double dT,double Tolerance,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
int GAPS_APT_Pusher_LCCSA_IMP_RFflow(double *F,double *X_iter,double *X_last,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
int GAPS_APT_Pusher_LCCSA_IMP_RFJacobi(double **J,double *X_iter,double *X_last,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
inline void tmpMatrixMult(double (*M1)[4],double (*M2)[4],double (*M)[4]);

int GAPS_APT_Pusher_LCCSA_IMP(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pCharge=GAPS_APT_GetCharge1(pPtc);
	double Tolerance=pInputs->Pusher_RootFindingTol;
	double dT=pInputs->dT;
	double *pX4 = GAPS_APT_GetX4(pPtc);
	double *pP4 = GAPS_APT_GetP4(pPtc);
	double *pA4 = GAPS_APT_GetA4(pPtc);

	double CanP4[4];
	(pPtc->FieldFunc)(pA4,pX4,1,pInputs);
	{int i;for(i=0;i<4;i++){CanP4[i]=pP4[i]+pCharge[0]*pA4[i];}}
	CanP4[0]*=-1.;

	double XP8[8]={pX4[0],pX4[1],pX4[2],pX4[3],CanP4[0],CanP4[1],CanP4[2],CanP4[3]};
	//Push P X
	GAPS_APT_Pusher_LCCSA_IMP_RootFindPush(XP8,dT,Tolerance,pPtc,pInputs);

	pX4[0]=XP8[0];
	pX4[1]=XP8[1];
	pX4[2]=XP8[2];
	pX4[3]=XP8[3];
	(pPtc->FieldFunc)(pA4,pX4,1,pInputs);

	pP4[0]=-XP8[4]-pCharge[0]*pA4[0];
	pP4[1]= XP8[5]-pCharge[0]*pA4[1];
	pP4[2]= XP8[6]-pCharge[0]*pA4[2];
	pP4[3]= XP8[7]-pCharge[0]*pA4[3];
	return 0;
}

int GAPS_APT_Pusher_LCCSA_IMP_RootFindPush(double *XP,double dT,double Tolerance,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double X_last[8]={XP[0],XP[1],XP[2],XP[3],XP[4],XP[5],XP[6],XP[7]};
	double X_iter[8]={XP[0],XP[1],XP[2],XP[3],XP[4],XP[5],XP[6],XP[7]};

	double **J=GenMatrixSpace(8,8),F[8];
	double tol,tolF;
	do
	{
		double delta_x[8];
		GAPS_APT_Pusher_LCCSA_IMP_RFflow(F,X_iter,X_last,dT,pPtc,pInputs);
		tolF=sqrt(F[0]*F[0]+F[1]*F[1]+F[2]*F[2]+F[3]*F[3]+F[4]*F[4]+F[5]*F[5]+F[6]*F[6]+F[7]*F[7]);
		GAPS_APT_Pusher_LCCSA_IMP_RFJacobi(J,X_iter,X_last,dT,pPtc,pInputs);
		double neg_F[8]={-1*F[0],-1*F[1],-1*F[2],-1*F[3],-1*F[4],-1*F[5],-1*F[6],-1*F[7]};
		LinearSolver_LU(delta_x,J,neg_F,8);

		double XP_new[8];
		double fnew;
		tol=0;
		int i;
		for(i=0;i<8;i++)
		{
			X_iter[i]+=delta_x[i];
			tol+=pow(delta_x[i],2);
		}
		tol=sqrt(tol);	
	}while(tol>=Tolerance && tolF>=Tolerance);
	
	int i;
	for(i=0;i<8;i++)
	{
		XP[i]=X_iter[i];
	}
	EraseMatrixSpace(J);
	return 0;
}

int GAPS_APT_Pusher_LCCSA_IMP_RFflow(double *F,double *X_iter,double *X_last,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pCharge=GAPS_APT_GetCharge1(pPtc);
	double *pMass=GAPS_APT_GetMass1(pPtc);
	double mass_rev=1./(*pMass);
	double qOverMass=(*pCharge)*mass_rev;
	double A[4],DelA[16];
	double q=(*pCharge);

	double xlast[4]={X_last[0],X_last[1],X_last[2],X_last[3]};	
	double plast[4]={X_last[4],X_last[5],X_last[6],X_last[7]};	
	double xn[4]={X_iter[0],X_iter[1],X_iter[2],X_iter[3]};	
	double pn[4]={X_iter[4],X_iter[5],X_iter[6],X_iter[7]};	
	double x_mid[4],p_mid[4];
	{
		int i;
		for(i=0;i<4;i++)
		{
			x_mid[i]=(xlast[i]+xn[i])/2.;	
			p_mid[i]=(plast[i]+pn[i])/2.;	
		}
	}
	(pPtc->FieldFunc)(A,x_mid,1,pInputs);
	(pPtc->FieldFunc)(DelA,x_mid,2,pInputs);

	int i;
	double Vx[4]={-p_mid[0]-q*A[0],p_mid[1]-q*A[1],p_mid[2]-q*A[2],p_mid[3]-q*A[3]};
	double V[4]={p_mid[0]+q*A[0],p_mid[1]-q*A[1],p_mid[2]-q*A[2],p_mid[3]-q*A[3]};
	for(i=0;i<8;i++)
	{
		if(i<4)
		{
			F[i] = 	(xn[i]-xlast[i]-dT*Vx[i]*mass_rev);	
		}
		else
		{
			int n=i-4;
			int j;
			double sumP=0;
			for(j=0;j<4;j++)
			{
				int IJ[2]={j,n};
				sumP+=V[j]*GAPS_APT_TensorValue(DelA,IJ,4,2);
			}
			F[i] = 	(pn[n]-plast[n]-dT*(sumP)*qOverMass);			
		}
	}
	return 0;
}

int GAPS_APT_Pusher_LCCSA_IMP_RFJacobi(double **J,double *X_iter,double *X_last,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pCharge=GAPS_APT_GetCharge1(pPtc);
	double *pMass=GAPS_APT_GetMass1(pPtc);
	double mass_rev=1./(*pMass);
	double qOverMass=(*pCharge)*mass_rev;
	double q=(*pCharge);

	double A[4],DelA1D[16],DelDelA[64];
	double xlast[4]={X_last[0],X_last[1],X_last[2],X_last[3]};	
	double plast[4]={X_last[4],X_last[5],X_last[6],X_last[7]};	
	double xn[4]={X_iter[0],X_iter[1],X_iter[2],X_iter[3]};	
	double pn[4]={X_iter[4],X_iter[5],X_iter[6],X_iter[7]};	
	double x_mid[4],p_mid[4];
	{
		int i;
		for(i=0;i<4;i++)
		{
			x_mid[i]=(xlast[i]+xn[i])/2.;	
			p_mid[i]=(plast[i]+pn[i])/2.;	
		}
	}
	(pPtc->FieldFunc)(A,x_mid,1,pInputs);
	(pPtc->FieldFunc)(DelA1D,x_mid,2,pInputs);
	(pPtc->FieldFunc)(DelDelA,x_mid,3,pInputs);

	double g[4][4],gDelA[4][4],DelA[4][4];
	g[0][0]=1;g[0][1]=0;g[0][2]=0;g[0][3]=0;
	g[1][0]=0;g[1][1]=-1;g[1][2]=0;g[1][3]=0;
	g[2][0]=0;g[2][1]=0;g[2][2]=-1;g[2][3]=0;
	g[3][0]=0;g[3][1]=0;g[3][2]=0;g[3][3]=-1;
	{
		int i,j;
		for(i=0;i<4;i++)
		{
			for(j=0;j<4;j++)
			{
				int IJ[2]={i,j};
				DelA[i][j]=GAPS_APT_TensorValue(DelA1D,IJ,4,2);
			}
		}
	}
	double V[4]={p_mid[0]+q*A[0],p_mid[1]-q*A[1],p_mid[2]-q*A[2],p_mid[3]-q*A[3]};
	tmpMatrixMult(g,DelA,gDelA);

	int i,j,k;
	for(i=0;i<4;i++)
	{
		for(j=0;j<4;j++)
		{
			if (i==j)
				J[i][j]=1+dT/2*DelA[i][j]*qOverMass;
			else
				J[i][j]=dT/2*DelA[i][j]*qOverMass;
		}
	}

	for(i=0;i<4;i++)
	{
		for(j=4;j<8;j++)
		{
			int ii=i,jj=j-4;
			J[i][j]=dT/2*g[ii][jj]*mass_rev;
		}
	}
	for(i=4;i<8;i++)
	{
		for(j=0;j<4;j++)
		{
			double sumVB=0.,sumDgD=0.;
			int ii=i-4,jj=j;

			for(k=0;k<4;k++)
			{
				int IJK[3]={k,ii,jj};
				sumVB+=V[k]*GAPS_APT_TensorValue(DelDelA,IJK,4,3);
			}
			for(k=0;k<4;k++)
			{
				sumDgD+=DelA[k][jj]*gDelA[k][ii]*q;
			}
			J[i][j]=-dT/2*(sumVB+sumDgD)*qOverMass;
		}
	}
	for(i=4;i<8;i++)
	{
		for(j=4;j<8;j++)
		{
			int ii=i-4,jj=j-4;
			if (i==j)
				J[i][j]=1-dT/2*DelA[jj][ii]*qOverMass;
			else
				J[i][j]=-dT/2*DelA[jj][ii]*qOverMass;
		}
	
	}
	return 0;
}

inline void tmpMatrixMult(double (*M1)[4],double (*M2)[4],double (*M)[4])
{
	long i,j,k;
	int dim=4;
	for(i=0;i<dim;i++)
	{
		for(j=0;j<dim;j++)
		{
			double tmp=0;
			for(k=0;k<dim;k++)
			{
				tmp+=M1[i][k]*M2[k][j];	
			}
			M[i][j]=tmp;
		}
	}
}

#include "APT_AllHeaders.h"
inline int PushP_CCS(double *P,double *A,double *DelA_array,double *PP,double dT,double q,double m);

int GAPS_APT_Pusher_LCCSA_SymEuler(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double dT=pInputs->dT;
	double *ptau = GAPS_APT_GetS1(pPtc);
	double *pX4 = GAPS_APT_GetX4(pPtc);
	double *pP4 = GAPS_APT_GetP4(pPtc);
	double *pA4 = GAPS_APT_GetA4(pPtc);

	double *pMass = GAPS_APT_GetMass1(pPtc);
	double *pCharge= GAPS_APT_GetCharge1(pPtc);

	double CanP4[4];
	double DelA[16];
	(pPtc->FieldFunc)(pA4,pX4,1,pInputs);
	(pPtc->FieldFunc)(DelA,pX4,2,pInputs);
	{int i;for(i=0;i<4;i++){CanP4[i]=pP4[i]+pCharge[0]*pA4[i];}}

	//Push CanP
	double CanP4_next[4];
	PushP_CCS(CanP4,pA4,DelA,CanP4_next,dT,pCharge[0],pMass[0]);

	{//Push X
		int i;
		for(i=0;i<4;i++)
		{
			pX4[i]+=dT*(CanP4_next[i]-pCharge[0]*pA4[i])/pMass[0];
		}
	}
	//Update A
	(pPtc->FieldFunc)(pA4,pX4,1,pInputs);

	{int i;for(i=0;i<4;i++){pP4[i]=CanP4_next[i]-pCharge[0]*pA4[i];}}
	ptau[0]+=dT;

	return 0;
}

inline int PushP_CCS(double *P,double *A,double *DelA_array,double *PP,double dT,double q,double m)
{	
	double qOverMdT=dT*q/m;
	double DelA[4][4];
	int idx[2];
	for(idx[0]=0;idx[0]<4;(idx[0])++)
	{
		for(idx[1]=0;idx[1]<4;(idx[1])++)
		{
			DelA[idx[0]][idx[1]]=qOverMdT*GAPS_APT_TensorValue(DelA_array,idx,4,2);	
		}
	}
	double M00=1-DelA[0][0],   M01=DelA[1][0],   M02=DelA[2][0],  M03=DelA[3][0];
	double M10=  DelA[0][1], M11=1-DelA[1][1],  M12=-DelA[2][1], M13=-DelA[3][1];
	double M20=  DelA[0][2],  M21=-DelA[1][2], M22=1-DelA[2][2], M23=-DelA[3][2];
	double M30=  DelA[0][3],  M31=-DelA[1][3],  M32=-DelA[2][3],M33=1-DelA[3][3];
	double V0=P[0]-(A[0]*DelA[0][0] - A[1]*DelA[1][0]- A[2]*DelA[2][0]- A[3]*DelA[3][0])*q ;
	double V1=P[1]+(A[0]*DelA[0][1] - A[1]*DelA[1][1]- A[2]*DelA[2][1]- A[3]*DelA[3][1])*q;
	double V2=P[2]+(A[0]*DelA[0][2] - A[1]*DelA[1][2]- A[2]*DelA[2][2]- A[3]*DelA[3][2])*q;
	double V3=P[3]+(A[0]*DelA[0][3] - A[1]*DelA[1][3]- A[2]*DelA[2][3]- A[3]*DelA[3][3])*q;

	double down;
	down=M01*M13*M22*M30 - M01*M12*M23*M30 - M00*M13*M22*M31 + M00*M12*M23*M31 - M01*M13*M20*M32 + M00*M13*M21*M32 + M01*M10*M23*M32 - M00*M11*M23*M32 + 
	M03*(M12*M21*M30 - M11*M22*M30 - M12*M20*M31 + M10*M22*M31 + M11*M20*M32 - M10*M21*M32) + M01*M12*M20*M33 - M00*M12*M21*M33 - M01*M10*M22*M33 + M00*M11*M22*M33 + 
	M02*(-(M13*M21*M30) + M11*M23*M30 + M13*M20*M31 - M10*M23*M31 - M11*M20*M33 + M10*M21*M33);
	PP[0]=(-(M11*M23*M32*V0) + M11*M22*M33*V0 + M03*M22*M31*V1 - M02*M23*M31*V1 - M03*M21*M32*V1 + M01*M23*M32*V1 + M02*M21*M33*V1 - M01*M22*M33*V1 + M03*M11*M32*V2 - M02*M11*M33*V2 - M03*M11*M22*V3 + 
			M02*M11*M23*V3 + M13*(-(M22*M31*V0) + M21*M32*V0 + M02*M31*V2 - M01*M32*V2 - M02*M21*V3 + M01*M22*V3) + M12*(M23*M31*V0 - M21*M33*V0 - M03*M31*V2 + M01*M33*V2 + M03*M21*V3 - M01*M23*V3))/down;
	PP[1]=(M10*M23*M32*V0 - M10*M22*M33*V0 - M03*M22*M30*V1 + M02*M23*M30*V1 + M03*M20*M32*V1 - M00*M23*M32*V1 - M02*M20*M33*V1 + M00*M22*M33*V1 - M03*M10*M32*V2 + M02*M10*M33*V2 + M03*M10*M22*V3 - 
			M02*M10*M23*V3 + M13*(M22*M30*V0 - M20*M32*V0 - M02*M30*V2 + M00*M32*V2 + M02*M20*V3 - M00*M22*V3) + M12*(-(M23*M30*V0) + M20*M33*V0 + M03*M30*V2 - M00*M33*V2 - M03*M20*V3 + M00*M23*V3))/down;
	PP[2]=(-(M10*M23*M31*V0) + M10*M21*M33*V0 + M03*M21*M30*V1 - M01*M23*M30*V1 - M03*M20*M31*V1 + M00*M23*M31*V1 + M01*M20*M33*V1 - M00*M21*M33*V1 + M03*M10*M31*V2 - M01*M10*M33*V2 - M03*M10*M21*V3 + 
			M01*M10*M23*V3 + M13*(-(M21*M30*V0) + M20*M31*V0 + M01*M30*V2 - M00*M31*V2 - M01*M20*V3 + M00*M21*V3) + M11*(M23*M30*V0 - M20*M33*V0 - M03*M30*V2 + M00*M33*V2 + M03*M20*V3 - M00*M23*V3))/down;
	PP[3]=(M10*M22*M31*V0 - M10*M21*M32*V0 - M02*M21*M30*V1 + M01*M22*M30*V1 + M02*M20*M31*V1 - M00*M22*M31*V1 - M01*M20*M32*V1 + M00*M21*M32*V1 - M02*M10*M31*V2 + M01*M10*M32*V2 + M02*M10*M21*V3 - 
			M01*M10*M22*V3 + M12*(M21*M30*V0 - M20*M31*V0 - M01*M30*V2 + M00*M31*V2 + M01*M20*V3 - M00*M21*V3) + M11*(-(M22*M30*V0) + M20*M32*V0 + M02*M30*V2 - M00*M32*V2 - M02*M20*V3 + M00*M22*V3))/down;
	return 0;	
}

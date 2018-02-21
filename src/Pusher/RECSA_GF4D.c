#include "APT_AllHeaders.h"
int GAPS_APT_Pusher_RECSA_GF4D_Order1(double *pX4,double *CanP4,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
int GAPS_APT_Pusher_RECSA_GF4D_Order2(double *pX4,double *CanP4,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
int GAPS_APT_Pusher_RECSA_GF4D_Order3(double *pX4,double *CanP4,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
int GAPS_APT_Pusher_RECSA_GF4D_Map_1(double *pX4,double *CanP4,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
int GAPS_APT_Pusher_RECSA_GF4D_Map_2(double *pX4,double *CanP4,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
int GAPS_APT_Pusher_RECSA_GF4D_Map_3(double *pX4,double *CanP4,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
int GAPS_APT_Pusher_RECSA_GF4D_Map_4(double *pX4,double *CanP4,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
int GAPS_APT_Pusher_RECSA_GF4D_Map_5(double *pX4,double *CanP4,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
int GAPS_APT_Pusher_RECSA_GF4D_Map_6(double *pX4,double *CanP4,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
int GAPS_APT_Pusher_RECSA_GF4D_Map_7(double *pX4,double *CanP4,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);

int GAPS_APT_Pusher_RECSA_GF4D(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double dT=pInputs->dT;
	double *ptau = GAPS_APT_GetS1(pPtc);
	double *pX4 = GAPS_APT_GetX4(pPtc);
	double *pP4 = GAPS_APT_GetP4(pPtc);
	double *pA4 = GAPS_APT_GetA4(pPtc);
	double *pMass = GAPS_APT_GetMass1(pPtc);
	double *pCharge= GAPS_APT_GetCharge1(pPtc);

	double CanP4[4];
	(pPtc->FieldFunc)(pA4,pX4,1,pInputs);
	{int i;for(i=0;i<4;i++){CanP4[i]=pP4[i]+pCharge[0]*pA4[i];}}

	CanP4[0]*=-1.;

	if(1==pInputs->Pusher_RECSA_GF4D_Order)
	{
		GAPS_APT_Pusher_RECSA_GF4D_Order1(pX4,CanP4,dT,pPtc,pInputs);
	}
	else if(2==pInputs->Pusher_RECSA_GF4D_Order)
	{
		GAPS_APT_Pusher_RECSA_GF4D_Order2(pX4,CanP4,dT,pPtc,pInputs);
	}
	else if(3==pInputs->Pusher_RECSA_GF4D_Order)
	{
		GAPS_APT_Pusher_RECSA_GF4D_Order3(pX4,CanP4,dT,pPtc,pInputs);
	}
	else
	{
		fprintf(stderr,"Wrong ORDER_ECSA4D for RECSA_4D\n");
	}

	(pPtc->FieldFunc)(pA4,pX4,1,pInputs);

	pP4[0]=-CanP4[0]-pCharge[0]*pA4[0];
	pP4[1]=CanP4[1]-pCharge[0]*pA4[1];
	pP4[2]=CanP4[2]-pCharge[0]*pA4[2];
	pP4[3]=CanP4[3]-pCharge[0]*pA4[3];
	return 0;
}
int GAPS_APT_Pusher_RECSA_GF4D_Order1(double *pX4,double *CanP4,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	GAPS_APT_Pusher_RECSA_GF4D_Map_1(pX4,CanP4,dT,pPtc,pInputs);
	GAPS_APT_Pusher_RECSA_GF4D_Map_2(pX4,CanP4,dT,pPtc,pInputs);
	GAPS_APT_Pusher_RECSA_GF4D_Map_3(pX4,CanP4,dT,pPtc,pInputs);
	GAPS_APT_Pusher_RECSA_GF4D_Map_4(pX4,CanP4,dT,pPtc,pInputs);
	GAPS_APT_Pusher_RECSA_GF4D_Map_5(pX4,CanP4,dT,pPtc,pInputs);
	GAPS_APT_Pusher_RECSA_GF4D_Map_6(pX4,CanP4,dT,pPtc,pInputs);
	GAPS_APT_Pusher_RECSA_GF4D_Map_7(pX4,CanP4,dT,pPtc,pInputs);
	return 0;
}

int GAPS_APT_Pusher_RECSA_GF4D_Order2(double *pX4,double *CanP4,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double dT2=dT*0.5;
	GAPS_APT_Pusher_RECSA_GF4D_Map_1(pX4,CanP4,dT2,pPtc,pInputs);
	GAPS_APT_Pusher_RECSA_GF4D_Map_2(pX4,CanP4,dT2,pPtc,pInputs);
	GAPS_APT_Pusher_RECSA_GF4D_Map_3(pX4,CanP4,dT2,pPtc,pInputs);
	GAPS_APT_Pusher_RECSA_GF4D_Map_4(pX4,CanP4,dT2,pPtc,pInputs);
	GAPS_APT_Pusher_RECSA_GF4D_Map_5(pX4,CanP4,dT2,pPtc,pInputs);
	GAPS_APT_Pusher_RECSA_GF4D_Map_6(pX4,CanP4,dT2,pPtc,pInputs);
	GAPS_APT_Pusher_RECSA_GF4D_Map_7(pX4,CanP4,dT, pPtc,pInputs);
	GAPS_APT_Pusher_RECSA_GF4D_Map_6(pX4,CanP4,dT2,pPtc,pInputs);
	GAPS_APT_Pusher_RECSA_GF4D_Map_5(pX4,CanP4,dT2,pPtc,pInputs);
	GAPS_APT_Pusher_RECSA_GF4D_Map_4(pX4,CanP4,dT2,pPtc,pInputs);
	GAPS_APT_Pusher_RECSA_GF4D_Map_3(pX4,CanP4,dT2,pPtc,pInputs);
	GAPS_APT_Pusher_RECSA_GF4D_Map_2(pX4,CanP4,dT2,pPtc,pInputs);
	GAPS_APT_Pusher_RECSA_GF4D_Map_1(pX4,CanP4,dT2,pPtc,pInputs);
	return 0;                                           
}

int GAPS_APT_Pusher_RECSA_GF4D_Order3(double *pX4,double *CanP4,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double a=1./(2.-pow(2.,1./3.));
	double b=1-2*a;
	GAPS_APT_Pusher_RECSA_GF4D_Order2(pX4,CanP4,a*dT,pPtc,pInputs);
	GAPS_APT_Pusher_RECSA_GF4D_Order2(pX4,CanP4,b*dT,pPtc,pInputs);
	GAPS_APT_Pusher_RECSA_GF4D_Order2(pX4,CanP4,a*dT,pPtc,pInputs);
	return 0;
}

int GAPS_APT_Pusher_RECSA_GF4D_Map_1(double *pX4,double *CanP4,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pMass = GAPS_APT_GetMass1(pPtc);
	int i;
	double mass_rev=1./pMass[0];
	for(i=1;i<4;i++)
	{
		pX4[i] += (dT*CanP4[i])*mass_rev;	
	}
	return 0;
}

int GAPS_APT_Pusher_RECSA_GF4D_Map_2(double *pX4,double *CanP4,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pMass = GAPS_APT_GetMass1(pPtc);
	double *pCharge= GAPS_APT_GetCharge1(pPtc);
	double qqOverMass=pCharge[0]*pCharge[0]/pMass[0];
	double mass_rev=1./pMass[0];

	double pA4[4];
	(pPtc->FieldFunc)(pA4,pX4,1,pInputs);
	double pDA4[16];
	(pPtc->FieldFunc)(pDA4,pX4,2,pInputs);
	double tmpGamma=0.;
	int i,j;
	int IJ[2];
	for(i=1;i<=3;i++)
	{
		IJ[0]=i;IJ[1]=0;tmpGamma+=pA4[i]*GAPS_APT_TensorValue(pDA4,IJ,4,2);
		double tmpP1=0.;
		for(j=1;j<=3;j++)
		{
			IJ[0]=j;IJ[1]=i;
			tmpP1+=pA4[j]*GAPS_APT_TensorValue(pDA4,IJ,4,2);
		}
		IJ[0]=0;IJ[1]=i;
		CanP4[i]+=qqOverMass*dT*(pA4[0]*GAPS_APT_TensorValue(pDA4,IJ,4,2)-tmpP1);
	}
	CanP4[0]+=qqOverMass*dT*(-tmpGamma+pA4[0]*pDA4[0]);
	return 0;
}

int GAPS_APT_Pusher_RECSA_GF4D_Map_3(double *pX4,double *CanP4,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pMass = GAPS_APT_GetMass1(pPtc);
	pX4[0]-=dT*CanP4[0]/pMass[0];
	return 0;
}

int GAPS_APT_Pusher_RECSA_GF4D_Map_4(double *pX4,double *CanP4,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{	
	double pA4[4];
	double pDA4[16];
	double pDDA4[64];
	(pPtc->FieldFunc)(pA4,pX4,1,pInputs);
	(pPtc->FieldFunc)(pDA4,pX4,2,pInputs);
	(pPtc->FieldFunc)(pDDA4,pX4,3,pInputs);
	//X={t,x,y,z); P={-gamma-phi,Px,Py,Pz}
	double *pMass = GAPS_APT_GetMass1(pPtc);
	double *pCharge= GAPS_APT_GetCharge1(pPtc);
	double qOverMass=pCharge[0]/pMass[0];
	double qOverMass_2=qOverMass*qOverMass;
	int FlowNo=1;

	double term_Gamma;
	double term_x;
	double term_p[3];

	int IJ[2],IJK[3];

	IJK[0]=FlowNo;IJK[1]=FlowNo;IJK[2]=0;
	double DDAii0=GAPS_APT_TensorValue(pDDA4,IJK,4,3);

	IJ[0]=FlowNo;IJ[1]=0;
	double DA_N0=GAPS_APT_TensorValue(pDA4,IJ,4,2);
	IJ[0]=FlowNo;IJ[1]=FlowNo;
	double DA_NN=GAPS_APT_TensorValue(pDA4,IJ,4,2);

	term_x=dT*(-1.*pA4[FlowNo]*qOverMass+qOverMass_2*dT*0.5*pA4[FlowNo]*DA_NN);
	term_Gamma=dT*(DA_N0*qOverMass-dT/2.*DA_NN*DA_N0*qOverMass_2-dT/2.*pA4[FlowNo]*DDAii0*qOverMass_2);
	{
		int i;
		for(i=1;i<=3;i++)
		{
			IJK[0]=FlowNo;IJK[1]=FlowNo;IJK[2]=i;
			double DDAiiI=GAPS_APT_TensorValue(pDDA4,IJK,4,3);
			IJ[0]=FlowNo;IJ[1]=i;
			double DA_NI=GAPS_APT_TensorValue(pDA4,IJ,4,2);
			term_p[i-1]=dT*(DA_NI*qOverMass-dT/2.*DA_NN*DA_NI*qOverMass_2-dT/2.*pA4[FlowNo]*DDAiiI*qOverMass_2);
		}
	}

	pX4[1]+=term_x;

//push p1
	double C_p1=1.-term_p[0];
	if(C_p1!=0)
	{
		CanP4[1]=CanP4[1]/C_p1;
	}
	else
	{
		fprintf(stderr,"Encounter Singular Point In Flow3\n");
	}
	CanP4[2]+=CanP4[1]*term_p[1];
	CanP4[3]+=CanP4[1]*term_p[2];
	CanP4[0]+=CanP4[1]*term_Gamma;
	return 0;
}

int GAPS_APT_Pusher_RECSA_GF4D_Map_5(double *pX4,double *CanP4,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double pA4[4];
	double pDA4[16];
	double pDDA4[64];
	(pPtc->FieldFunc)(pA4,pX4,1,pInputs);
	(pPtc->FieldFunc)(pDA4,pX4,2,pInputs);
	(pPtc->FieldFunc)(pDDA4,pX4,3,pInputs);
	double *pMass = GAPS_APT_GetMass1(pPtc);
	double *pCharge= GAPS_APT_GetCharge1(pPtc);
	double qOverMass=pCharge[0]/pMass[0];
	double qOverMass_2=qOverMass*qOverMass;
	int FlowNo=2;

	double term_Gamma;
	double term_x;
	double term_p[3];

	int IJ[2],IJK[3];

	IJK[0]=FlowNo;IJK[1]=FlowNo;IJK[2]=0;
	double DDAii0=GAPS_APT_TensorValue(pDDA4,IJK,4,3);

	IJ[0]=FlowNo;IJ[1]=0;
	double DA_N0=GAPS_APT_TensorValue(pDA4,IJ,4,2);
	IJ[0]=FlowNo;IJ[1]=FlowNo;
	double DA_NN=GAPS_APT_TensorValue(pDA4,IJ,4,2);

	term_x=dT*(-1.*pA4[FlowNo]*qOverMass+qOverMass_2*dT*0.5*pA4[FlowNo]*DA_NN);
	term_Gamma=dT*(DA_N0*qOverMass-dT/2.*DA_NN*DA_N0*qOverMass_2-dT/2.*pA4[FlowNo]*DDAii0*qOverMass_2);
	{
		int i;
		for(i=1;i<=3;i++)
		{
			IJK[0]=FlowNo;IJK[1]=FlowNo;IJK[2]=i;
			double DDAiiI=GAPS_APT_TensorValue(pDDA4,IJK,4,3);
			IJ[0]=FlowNo;IJ[1]=i;
			double DA_NI=GAPS_APT_TensorValue(pDA4,IJ,4,2);
			term_p[i-1]=dT*(DA_NI*qOverMass-dT/2.*DA_NN*DA_NI*qOverMass_2-dT/2.*pA4[FlowNo]*DDAiiI*qOverMass_2);
		}
	}

	pX4[2]+=term_x;

//push p1
	double C_p2=1.-term_p[1];
	if(C_p2!=0)
	{
		CanP4[2]=CanP4[2]/C_p2;
	}
	else
	{
		fprintf(stderr,"Encounter Singular Point In Flow4\n");
	}
	CanP4[1]+=CanP4[2]*term_p[0];
	CanP4[3]+=CanP4[2]*term_p[2];
	CanP4[0]+=CanP4[2]*term_Gamma;
	return 0;
}

int GAPS_APT_Pusher_RECSA_GF4D_Map_6(double *pX4,double *CanP4,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double pA4[4];
	double pDA4[16];
	double pDDA4[64];
	(pPtc->FieldFunc)(pA4,pX4,1,pInputs);
	(pPtc->FieldFunc)(pDA4,pX4,2,pInputs);
	(pPtc->FieldFunc)(pDDA4,pX4,3,pInputs);
	double *pMass = GAPS_APT_GetMass1(pPtc);
	double *pCharge= GAPS_APT_GetCharge1(pPtc);
	double qOverMass=pCharge[0]/pMass[0];
	double qOverMass_2=qOverMass*qOverMass;
	int FlowNo=3;

	double term_Gamma;
	double term_x;
	double term_p[3];

	int IJ[2],IJK[3];

	IJK[0]=FlowNo;IJK[1]=FlowNo;IJK[2]=0;
	double DDAii0=GAPS_APT_TensorValue(pDDA4,IJK,4,3);

	IJ[0]=FlowNo;IJ[1]=0;
	double DA_N0=GAPS_APT_TensorValue(pDA4,IJ,4,2);
	IJ[0]=FlowNo;IJ[1]=FlowNo;
	double DA_NN=GAPS_APT_TensorValue(pDA4,IJ,4,2);

	term_x=dT*(-1.*pA4[FlowNo]*qOverMass+qOverMass_2*dT*0.5*pA4[FlowNo]*DA_NN);
	term_Gamma=dT*(DA_N0*qOverMass-dT/2.*DA_NN*DA_N0*qOverMass_2-dT/2.*pA4[FlowNo]*DDAii0*qOverMass_2);
	{
		int i;
		for(i=1;i<=3;i++)
		{
			IJK[0]=FlowNo;IJK[1]=FlowNo;IJK[2]=i;
			double DDAiiI=GAPS_APT_TensorValue(pDDA4,IJK,4,3);
			IJ[0]=FlowNo;IJ[1]=i;
			double DA_NI=GAPS_APT_TensorValue(pDA4,IJ,4,2);
			term_p[i-1]=dT*(DA_NI*qOverMass-dT/2.*DA_NN*DA_NI*qOverMass_2-dT/2.*pA4[FlowNo]*DDAiiI*qOverMass_2);
		}
	}

	pX4[3]+=term_x;

	double C_p3=1.-term_p[2];
	if(C_p3!=0)
	{
		CanP4[3]=CanP4[3]/C_p3;
	}
	else
	{
		fprintf(stderr,"Encounter Singular Point In Flow4\n");
	}
	CanP4[1]+=CanP4[3]*term_p[0];
	CanP4[2]+=CanP4[3]*term_p[1];
	CanP4[0]+=CanP4[3]*term_Gamma;
	return 0;
}

int GAPS_APT_Pusher_RECSA_GF4D_Map_7(double *pX4,double *CanP4,double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double pA4[4];
	double pDA4[16];
	double pDDA4[64];
	(pPtc->FieldFunc)(pA4,pX4,1,pInputs);
	(pPtc->FieldFunc)(pDA4,pX4,2,pInputs);
	(pPtc->FieldFunc)(pDDA4,pX4,3,pInputs);
	double *pMass = GAPS_APT_GetMass1(pPtc);
	double *pCharge= GAPS_APT_GetCharge1(pPtc);
	double qOverMass=pCharge[0]/pMass[0];
	double qOverMass_2=qOverMass*qOverMass;
	int FlowNo=0;

	double term_Gamma;
	double term_x;
	double term_p[3];

	int IJ[2],IJK[3];

	IJK[0]=FlowNo;IJK[1]=FlowNo;IJK[2]=0;
	double DDAii0=GAPS_APT_TensorValue(pDDA4,IJK,4,3);

	IJ[0]=FlowNo;IJ[1]=0;
	double DA_N0=GAPS_APT_TensorValue(pDA4,IJ,4,2);
	IJ[0]=FlowNo;IJ[1]=FlowNo;
	double DA_NN=GAPS_APT_TensorValue(pDA4,IJ,4,2);

	term_x=dT*(-1.*pA4[FlowNo]*qOverMass+qOverMass_2*dT*0.5*pA4[FlowNo]*DA_NN);
	term_Gamma=dT*(DA_N0*qOverMass-dT/2.*DA_NN*DA_N0*qOverMass_2-dT/2.*pA4[FlowNo]*DDAii0*qOverMass_2);
	{
		int i;
		for(i=1;i<=3;i++)
		{
			IJK[0]=FlowNo;IJK[1]=FlowNo;IJK[2]=i;
			double DDAiiI=GAPS_APT_TensorValue(pDDA4,IJK,4,3);
			IJ[0]=FlowNo;IJ[1]=i;
			double DA_NI=GAPS_APT_TensorValue(pDA4,IJ,4,2);
			term_p[i-1]=dT*(DA_NI*qOverMass-dT/2.*DA_NN*DA_NI*qOverMass_2-dT/2.*pA4[FlowNo]*DDAiiI*qOverMass_2);
		}
	}

	pX4[0]+=term_x;

	double C_p0=1.-term_Gamma;
	if(C_p0!=0)
	{
		CanP4[0]=CanP4[0]/C_p0;
	}
	else
	{
		fprintf(stderr,"Encounter Singular Point In Flow4\n");
	}
	CanP4[1]+=CanP4[0]*term_p[0];
	CanP4[2]+=CanP4[0]*term_p[1];
	CanP4[3]+=CanP4[0]*term_p[2];
	return 0;
}

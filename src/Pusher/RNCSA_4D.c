#include "APT_AllHeaders.h"
int GAPS_APT_Pusher_RNCSA_4D_MapT(double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
int GAPS_APT_Pusher_RNCSA_4D_MapX(double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
int GAPS_APT_Pusher_RNCSA_4D_MapY(double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
int GAPS_APT_Pusher_RNCSA_4D_MapZ(double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);

int GAPS_APT_Pusher_RNCSA_4D(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pXP8 = GAPS_APT_GetX4(pPtc);
	double *pMass = GAPS_APT_GetMass1(pPtc);

	double dT=pInputs->dT;
	double P3[3]={pXP8[5],pXP8[6],pXP8[7]};
	double F_ext1[3],gamma1=pXP8[4];
	GAPS_APT_MergeExtForce(F_ext1,pPtc,pInputs);
	
	GAPS_APT_Pusher_RNCSA_4D_MapX(0.5*dT,pPtc,pInputs);
	GAPS_APT_Pusher_RNCSA_4D_MapY(0.5*dT,pPtc,pInputs);
	GAPS_APT_Pusher_RNCSA_4D_MapZ(0.5*dT,pPtc,pInputs);
	GAPS_APT_Pusher_RNCSA_4D_MapT(dT,pPtc,pInputs);
	GAPS_APT_Pusher_RNCSA_4D_MapZ(0.5*dT,pPtc,pInputs);
	GAPS_APT_Pusher_RNCSA_4D_MapY(0.5*dT,pPtc,pInputs);
	GAPS_APT_Pusher_RNCSA_4D_MapX(0.5*dT,pPtc,pInputs);
	
	//Update and Return sum of all extern forces
	double F_ext2[3],gamma2=pXP8[4];
	GAPS_APT_MergeExtForce(F_ext2,pPtc,pInputs);

	//Add external force effects using midpoint method
//	double F_ext[3]={(F_ext1[0]+F_ext2[0])*.5,(F_ext1[1]+F_ext2[1])*.5,(F_ext1[2]+F_ext2[2])*.5};
//	double gamma_F_ext[3]={(gamma1*F_ext1[0]+gamma2*F_ext2[0])*.5,(gamma1*F_ext1[1]+gamma2*F_ext2[1])*.5,(gamma1*F_ext1[2]+gamma2*F_ext2[2])*.5};
//	double P3_next[3]={(P3[0]+pXP8[5])*.5,(P3[1]+pXP8[6])*.5,(P3[2]+pXP8[7])*.5};
//	pXP8[4]+=dT/pMass[0]*(P3_next[0]*F_ext[0]+P3_next[1]*F_ext[1]+P3_next[2]*F_ext[2]);
//	pXP8[5]+=dT*gamma_F_ext[0];
//	pXP8[6]+=dT*gamma_F_ext[1];
//	pXP8[7]+=dT*gamma_F_ext[2];
	return 0;
}

int GAPS_APT_Pusher_RNCSA_4D_MapT(double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pX4= GAPS_APT_GetX4(pPtc);
	double *pP4= GAPS_APT_GetP4(pPtc);
	double *pMass = GAPS_APT_GetMass1(pPtc);
	double *pCharge= GAPS_APT_GetCharge1(pPtc);

	double z0=pX4[0];
	pX4[0]+=dT*pP4[0];
	double z1=pX4[0];
	
	//Numerical integral
	int n =pInputs->Pusher_NIntegral_N;
	int Order = -1; // Order of EB is -1
	int idxInt =0;
	double pST0[4]={ z0,pX4[1],pX4[2],pX4[3]};
	double pST1[4]={ z1,pX4[1],pX4[2],pX4[3]};
	double IntE[3];

	int pIndex[1];
	pIndex[0]=0; //Ex
	IntE[0]= GAPS_APT_IntegralField(n,pPtc,Order,pIndex,idxInt,pST0,pST1,pInputs);
	pIndex[0]=1; //Ey
	IntE[1]= GAPS_APT_IntegralField(n,pPtc,Order,pIndex,idxInt,pST0,pST1,pInputs);
	pIndex[0]=2; //Ez
	IntE[2]= GAPS_APT_IntegralField(n,pPtc,Order,pIndex,idxInt,pST0,pST1,pInputs);

	pP4[1]+=IntE[0]*pCharge[0];
	pP4[2]+=IntE[1]*pCharge[0];
	pP4[3]+=IntE[2]*pCharge[0];
	return 0;
}

int GAPS_APT_Pusher_RNCSA_4D_MapX(double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pX4= GAPS_APT_GetX4(pPtc);
	double *pP4= GAPS_APT_GetP4(pPtc);
	double *pMass = GAPS_APT_GetMass1(pPtc);
	double *pCharge= GAPS_APT_GetCharge1(pPtc);

	double z0=pX4[1];
	pX4[1]+=dT*pP4[1]/pMass[0];
	double z1=pX4[1];
//	printf("dT_in=%e,z0=%e,z1=%e\n",dT,z0,z1);
	
	//Numerical integral
	int n =pInputs->Pusher_NIntegral_N;
	int Order = -1; // Order of EB is -1
	int idxInt =1;
	double pST0[4]={ pX4[0],z0,pX4[2],pX4[3]};
	double pST1[4]={ pX4[0],z1,pX4[2],pX4[3]};
	int pIndex[1]; 

	//gamma
	pIndex[0]=0;//Ex
	pP4[0]+= pCharge[0]*GAPS_APT_IntegralField(n,pPtc,Order,pIndex,idxInt,pST0,pST1,pInputs);

	//py
	pIndex[0]=5;//Bz
	double tmp;
	pP4[2]+= pCharge[0]*GAPS_APT_IntegralField(n,pPtc,Order,pIndex,idxInt,pST0,pST1,pInputs);
//	tmp= pCharge[0]*GAPS_APT_IntegralField(n,pPtc,Order,pIndex,idxInt,pST0,pST1,pInputs);
//	printf("intBz_x=%e,solve=%e\n",tmp,dT*pP4[1]/pMass[0]);


	//pz
	pIndex[0]=4;//By
	pP4[3]-= pCharge[0]*GAPS_APT_IntegralField(n,pPtc,Order,pIndex,idxInt,pST0,pST1,pInputs);

	return 0;
}

int GAPS_APT_Pusher_RNCSA_4D_MapY(double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pX4= GAPS_APT_GetX4(pPtc);
	double *pP4= GAPS_APT_GetP4(pPtc);
	double *pMass = GAPS_APT_GetMass1(pPtc);
	double *pCharge= GAPS_APT_GetCharge1(pPtc);

	double z0=pX4[2];
	pX4[2]+=dT*pP4[2]/pMass[0];
	double z1=pX4[2];
	
	//Numerical integral
	int n =pInputs->Pusher_NIntegral_N;
	int Order = -1; // Order of EB is -1
	int idxInt =2;
	double pST0[4]={ pX4[0],pX4[1],z0,pX4[3]};
	double pST1[4]={ pX4[0],pX4[1],z1,pX4[3]};
	int pIndex[1]; 

	//gamma
	pIndex[0]=1;//Ey
	pP4[0]+= pCharge[0]*GAPS_APT_IntegralField(n,pPtc,Order,pIndex,idxInt,pST0,pST1,pInputs);

	//px
	pIndex[0]=5;//Bz
	pP4[1]-= pCharge[0]*GAPS_APT_IntegralField(n,pPtc,Order,pIndex,idxInt,pST0,pST1,pInputs);

	//pz
	pIndex[0]=3;//Bx
	pP4[3]+= pCharge[0]*GAPS_APT_IntegralField(n,pPtc,Order,pIndex,idxInt,pST0,pST1,pInputs);
	return 0;
}

int GAPS_APT_Pusher_RNCSA_4D_MapZ(double dT,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pX4= GAPS_APT_GetX4(pPtc);
	double *pP4= GAPS_APT_GetP4(pPtc);
	double *pMass = GAPS_APT_GetMass1(pPtc);
	double *pCharge= GAPS_APT_GetCharge1(pPtc);

	double z0=pX4[3];
	pX4[3]+=dT*pP4[3]/pMass[0];
	double z1=pX4[3];
	
	//Numerical integral
	int n =pInputs->Pusher_NIntegral_N;
	int Order = -1; // Order of EB is -1
	int idxInt =3;
	double pST0[4]={ pX4[0],pX4[1],pX4[2],z0};
	double pST1[4]={ pX4[0],pX4[1],pX4[2],z1};
	int pIndex[1]; 

	//gamma
	pIndex[0]=2;//Ez
	pP4[0]+= pCharge[0]*GAPS_APT_IntegralField(n,pPtc,Order,pIndex,idxInt,pST0,pST1,pInputs);

	//px
	pIndex[0]=4;//By
	pP4[1]+= pCharge[0]*GAPS_APT_IntegralField(n,pPtc,Order,pIndex,idxInt,pST0,pST1,pInputs);

	//py
	pIndex[0]=3;//Bx
	pP4[2]-= pCharge[0]*GAPS_APT_IntegralField(n,pPtc,Order,pIndex,idxInt,pST0,pST1,pInputs);
	return 0;
}

// Add Pushers Here!
// This file contains definitions related to physical models
// Constants
typedef int (*Gaps_APT_ParticlePusher) (Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);

int GAPS_APT_SetParticlePusher(Gaps_APT_ParticlePusher *pPusher,Gaps_IO_InputsContainer *pInputs);
inline int GAPS_APT_UpdatePtcData_EB(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
inline int GAPS_APT_UpdatePtcData_A4(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
inline int GAPS_APT_CalEB(double *E,double *B,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
inline int GAPS_APT_CalA(double *A,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
inline int GAPS_APT_Push_ParticleGroup(Gaps_APT_ParticleGroup *pPtcgrp,Gaps_APT_ParticlePusher pPusher,Gaps_IO_InputsContainer *pInputs);
inline int GAPS_APT_SavePhaseSpace2Cache(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
inline int GAPS_APT_CalAcceleration(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
// Pushers
//void RVPA_cay (double *Px, double *Pp, double *Pv,ForcePointer *Pforce);
//void GAPS_APT_Pusher_CanonSymp4D(SParticle *single_ptc);//4D explicit canonical sympletic algorithm
//void GAPS_APT_Pusher_CanonSymp3D(SParticle *single_ptc);//3D implicit canonical sympletic algorithm
//void GAPS_APT_Pusher_RVPA_Cay4D(SParticle *single_ptc);//Covariant Relativistic Volume-preserving Algorithm -- Cayley Transformation
//void GAPS_APT_Pusher_RVPA_exp(SParticle *single_ptc);
//void GAPS_APT_Pusher_RungeKutta(SParticle *single_ptc);
//void GAPS_APT_Pusher_RNonCS_4D(SParticle *single_ptc);
//void GAPS_APT_Pusher_RVPA_Cay3D_xboost (SParticle *single_ptc);
//void GAPS_APT_Pusher_RECSA_4D(SParticle *single_ptc);
//void GAPS_APT_Pusher_CanonSymp4D_IMP(SParticle *single_ptc);
//void GAPS_APT_Pusher_Canon_NonCovTest_4D(SParticle *single_ptc);

// Aux-functions for Pushers

//void p_minus2p_plus(double dT,double *Ppm,double *PBn, double *Ppp);		// Cayley transformation part for RVPA_cay
//void PushP_CCS(double *P,double *A,double **DelA,double *PP);
//void RFfunc_CanonSymp3D(double *X_iter,double *X_last,double *A,double **DelA,double *DelPHI, double *F);
//void RFJacobi_CanonSymp3D(double *X_iter,double *X_last,double *A,double **DelA,double *DelPHI, double **J);
//void p_minus2p_plus_exp(double dT,double *Pm,double *B, double *Pp);
//void LorentzFlow(SParticle *ptc,double t, double *y, long Dim, double *F);
//
//void Map0_RNCS4D(double *X,double DeltaTau,Int_EB_Field IntEB_using);
//void MapX_RNCS4D(double *X,double DeltaTau,Int_EB_Field IntEB_using);
//void MapY_RNCS4D(double *X,double DeltaTau,Int_EB_Field IntEB_using);
//void MapZ_RNCS4D(double *X,double DeltaTau,Int_EB_Field IntEB_using);
//
//void AuxP_RVPA_Cay3D_xboost(SParticle *single_ptc,double *tprime,double *xprime,double *gamma_prime,double *pprime,double *mathcalP);
//
//
//void Map1_RECSA_4D(double *X,double *P,double DeltaTau,SParticle *single_ptc);
//void Map2_RECSA_4D(double *X,double *P,double DeltaTau,SParticle *single_ptc);
//void Map3_RECSA_4D(double *X,double *P,double DeltaTau,SParticle *single_ptc);
//void Map4_RECSA_4D(double *X,double *P,double DeltaTau,SParticle *single_ptc);
//void Map5_RECSA_4D(double *X,double *P,double DeltaTau,SParticle *single_ptc);
//void Map6_RECSA_4D(double *X,double *P,double DeltaTau,SParticle *single_ptc);
//void Map7_RECSA_4D(double *X,double *P,double DeltaTau,SParticle *single_ptc);
//void Map_RECSA_4D_1(double *X,double *P,double DeltaTau,SParticle *single_ptc);
//void Map_RECSA_4D_2(double *X,double *P,double DeltaTau,SParticle *single_ptc);
//void Map_RECSA_4D_3(double *X,double *P,double DeltaTau,SParticle *single_ptc);
//
//
//void RFfunc_CanonSymp4D_IMP(SParticle *single_ptc,double *X_iter,double *X_last,double *F);
//void RFJacobi_CanonSymp4D_IMP(SParticle *single_ptc,double *X_iter,double *X_last,double **J);
//
//
//void PushP_NonCov_Euler(double *P,double *A,double **DelA,double *PP);
//void PushX_NonCov_Break_P(double *X,double *P,double *A,double *PP,double *XX);


//Generated headers
#include "Pusher_Prototype.h"

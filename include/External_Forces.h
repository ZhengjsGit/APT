typedef int (*Gaps_APT_ExtForce)(double *pField,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);

inline int GAPS_APT_CalExtForce(double *pForce,int Type,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
inline int GAPS_APT_UpdatePtcData_ExtForce(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);
inline int GAPS_APT_MergeExtForce(double *F_ext,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);

inline int GAPS_APT_CalForceWorks(Gaps_APT_Particle *pPtc, Gaps_IO_InputsContainer *pInputs);

//Generated header
#include "ExtForce_Prototype.h"

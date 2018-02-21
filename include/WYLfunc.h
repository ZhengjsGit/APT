// Functions defined in this file can be used everywhere

inline double LorentzMetric(int i,int j);
inline double Delta(int i,int j);
void LinspaceInt (long a,long b,long n,long *p);// Like "linspace()" in Matlab, this function returns int data
void ProducePath (char *from, char *to);

// This file contains definitions of functions that implement fundamental calculation of vectors(3-D double arrays).

// 3D-Vector operators
void V3plus   (double *A,double *B,double *C);		//C=A+B
void V3minus  (double *A,double *B,double *C);		//C=A-B
void V3cross  (double *A,double *B,double *C);		//C=A X B
void V3Smult  (double s,double *A,double *B);		//B=s*A
double V3INmult  (double *A,double *B);			//A.B
double V3norm   (double *A);				//|A|


//Coordinent transformation
void CART2CYLD (double *pf,double *pt);
void CART2TORD (double *pf,double *pt,double R0);
void CYLD2TORD (double *pc,double *pt,double R0);
void CYLD2CART (double *pf,double *pt);
void TORD2CART (double *pf,double *pt,double R0);
void TORD2CYLD (double *pf,double *pt,double R0);



//Matrix and its operations
void M3multV3 (double (*M)[3], double *A, double *B);	//MA=B
void M3multM3 (double (*M1)[3], double (*M2)[3], double (*M3)[3]);	//M1 M2 = M3


// Get local frame by local magnetic field
void LocalFrame_B(double *B, double *e1,double *e2, double *b);

// Tranform (p_\parallel,p_\perp,phase) to (px,py,pz)
void p_GC_TO_p_CART(double *p0,double *B,double *p);

int SetNum4String(char *str, long num);
double GenRandNum_Uniform(double down,double up);
double func_normal(double x, double miu,double sigma);
double GenRandNum_Norm(double miu,double sigma);
double func_circ(double R, double Rmax,double Rmin);
double GenRandNum_Circ(double Rmin,double Rmax);
double func_para(double x, double m);
double GenRandNum_Para(double m);


//Linear Algebra

double ** GenMatrixSpace(long M,long N);
int EraseMatrixSpace(double **A);
double *** GenTensorSpace3(long Dim);
inline void EraseTensorSpace3(double ***A);

//void LUdecomp(double **L,double **U,double *P,double **A,long N);
//int LinearSolver_LUD(double *x, double **A, double *y, long N);
void LorentzBoost(double *b, double **A,int Inv,double SmallEnough);
void LinearTransform(double **M,double *V,double *Vp,long dim);
void MatrixMult(double **M1,double **M2,double **M,long dim);
void MatrixTranspose(double **M,double **N,long dim);

int LUPDecompose_wiki(double **A,long N,double Tol,long *P);
int LUPSolve_wiki(double **A,long *P,double *b,long N,double *x);
int LUPInvert_wiki(double **A,long *P,long N,double **IA);
double LUPDeterminant_wiki(double **A,long *P,long N);
int LinearSolver_LU(double *x, double **A, double *b, long N);
double Determinant_LU(double **A, long N);
int Invert_LU(double ** InvA,double **A, long N);


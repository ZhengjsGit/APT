// Functions defined in this file can be used everywhere
#include <stdio.h>
#include<string.h>
#include<stdlib.h>
#include <math.h>
#include "WYLfunc.h" 
#define EXTRAME_SMALL 1e-50

inline double LorentzMetric(int i,int j)
{
	if (i==j)
	{
		if(i==0)
		{
			return 1.;
		}
		else
		{
			return -1.;
		}
	}
	else
	{
		return 0.;
	}
}
inline double Delta(int i,int j)
{
	if (i==j)
	{
		return 1.;
	}
	else
	{
		return 0.;
	}
}

void LinspaceInt(long a,long b,long n,long *p){

	long  *temp=p,i;
	if (n == 0)
		n =1;
	if (n <= b-a+1){
		for(i=0;i<n-1;i++){
			
			*(temp+i) =floor( a+(double) (b-a)/(n-1)*i);
		}	
		*(temp+n-1)=b;
	}
	else
		for(i=0;i<n;i++)
			{
			*(temp+i) =a+i;
		}	
}

void ProducePath (char *from, char *to)
{
	char * position;
	int length,pp;

	position=strrchr(from, '/');
	
	if (position)
	{
		pp = (int)(position-from);
		length=strlen(from);

		if(pp+1<length)
		{
			strcpy(to, from);
			strcat(to,"/");
		}
		else
			strcpy(to,from);
	}
	else
	{
		strcpy(to,from);
		strcat(to,"/");
	}


}


// 3D-Vector operators

void V3plus(double *A,double *B,double *C){

	*C = (*A)+(*B);
	(*(C+1)) = (*(A+1))+(*(B+1));
	(*(C+2)) = (*(A+2))+(*(B+2));
}
void V3minus(double *A,double *B,double *C){

	*C = (*A) - (*B);
	(*(C+1)) = (*(A+1)) - (*(B+1));
	(*(C+2)) = (*(A+2)) - (*(B+2));
}
void V3cross(double *A,double *B,double *C){

	double A1,A2,A3,B1,B2,B3;
	A1 = (*A);
	A2 = (*(A+1));
	A3 = (*(A+2));
	B1 = (*B);
	B2 = (*(B+1));
	B3 = (*(B+2));

	(*C) = A2*B3-A3*B2;
	(*(C+1)) = A3*B1-A1*B3;
	(*(C+2)) = A1*B2-A2*B1;
}
void V3Smult(double a,double *A,double *B){

	*B = (*A) * a; 
	*(B+1) = (*(A+1)) * a; 
	*(B+2) = (*(A+2)) * a; 
}
double V3INmult(double *A, double *B){
	double t; 
	t = (*A) * (*B)+ (*(A+1)) * (*(B+1)) + (*(A+2)) * (*(B+2));
	return t; 
}
double V3norm(double *A){

	return sqrt(pow(*A,2)+pow(*(A+1),2)+pow(*(A+2),2));
}



//Coordinent transformation
void CART2CYLD (double *pf, double *pt)
{
	double x=(*pf),y=(*(pf+1)),R,Theta;

	R = sqrt(pow(x,2)+pow(y,2));
	if (R == 0)
	{
		Theta = 0;
		printf("Singular point occors when CART to CYLD\n");
	}
	else if(y >= 0)
		Theta = acos(x / R);
	else
		Theta = 2*M_PI-acos(x/R);
		
	*pt = R;
	*(pt+1)=Theta;
	*(pt+2)= (*(pf+2));
}
void CART2TORD (double *pf,double *pt,double R0){
	double temp[3],*ptemp=temp;
	CART2CYLD(pf,ptemp);
	CYLD2TORD(ptemp,pt,R0);
}
void CART2SPHER (double *pf, double *pt)
{
	double x=(*(pf)),y=(*(pf+1)),z=(*(pf+2));
	double r=sqrt(x*x+y*y+z*z),l=sqrt(x*x+y*y);

//	printf("Cart2spher: x=%e,y=%e,z=%e,r=%e,l=%e\n",x,y,z,r,l);
	if (r==0)
		printf("Error! Reach to the center sigular point of spher-cord\n");
	else
	{
		*pt = r;

		*(pt+1)=acos(z/r);
			
		if (l == 0)
		{
			*(pt+2) = 0;
			printf("Phi-Singular point occors when CART to CYLD\n");
		}
		else 
		{	if(y >= 0)
			*(pt+2)= acos(x / l);
			else
			*(pt+2) = 2*M_PI-acos(x/l);

		}

	}
}

void SPHER2CART (double *pf, double *pt)
{
	double r=(*(pf)),theta=(*(pf+1)),phi=(*(pf+2));

	*pt = r*sin(theta)*cos(phi);
	*(pt+1)=r*sin(theta)*sin(phi);
	*(pt+2)=r*cos(theta);
}


void CYLD2TORD (double *pc,double *pt,double R0){

	*pt = sqrt(pow((*pc) - R0,2)+pow((*(pc+2)),2));
	if ((*pt) == 0)
		*(pt+1) = 0;
	else if((*(pc+2)) >= 0)
		*(pt+1) = acos(((*pc) - R0)/(*pt));
	else
		*(pt+1) = 2*M_PI- acos(((*pc) - R0)/(*pt));
		
	*(pt+2) = -1*(*(pc+1));	
}
void CYLD2CART (double *pf,double *pt){
	
	*pt     = (*pf) * cos(*(pf+1));
	*(pt+1) = (*pf) * sin(*(pf+1));
	*(pt+2) = (*(pf+2));
}
void TORD2CART (double *pf,double *pt,double R0){
	double temp[3],*ptemp=temp;
	TORD2CYLD(pf,ptemp,R0);
	CYLD2CART(ptemp,pt);	
}
void TORD2CYLD (double *pf,double *pt,double R0){
	
	*pt = R0 + (*pf) * cos(*(pf+1));
	*(pt+1) = -1 * (*(pf+2));	
	*(pt+2) = (*pf) * sin(*(pf+1));
}



//Matrix and its operations
void M3multV3 (double (*M)[3], double *A, double *B){

	*B	= (*(*M)) * (*A) +(*(*M +1)) * (*(A+1)) + (*(*M + 2)) * (*(A+2));
	*(B+1)	= (*(*(M+1))) * (*A) +(*(*(M+1) +1)) * (*(A+1)) + (*(*(M+1) + 2)) * (*(A+2));
	*(B+2)	= (*(*(M+2))) * (*A) +(*(*(M+2) +1)) * (*(A+1)) + (*(*(M+2) + 2)) * (*(A+2));
}

void M3multM3 (double (*M1)[3], double (*M2)[3], double (*M3)[3]){
	
	*(*(M3+0)+0)	= (*(*(M1+0)+0)) * (*(*(M2+0)+0)) +(*(*(M1+0)+1)) * (*(*(M2+1)+0)) + (*(*(M1+0)+2)) * (*(*(M2+2)+0));
	*(*(M3+0)+1)	= (*(*(M1+0)+0)) * (*(*(M2+0)+1)) +(*(*(M1+0)+1)) * (*(*(M2+1)+1)) + (*(*(M1+0)+2)) * (*(*(M2+2)+1));
	*(*(M3+0)+2)	= (*(*(M1+0)+0)) * (*(*(M2+0)+2)) +(*(*(M1+0)+1)) * (*(*(M2+1)+2)) + (*(*(M1+0)+2)) * (*(*(M2+2)+2));
	*(*(M3+1)+0)	= (*(*(M1+1)+0)) * (*(*(M2+0)+0)) +(*(*(M1+1)+1)) * (*(*(M2+1)+0)) + (*(*(M1+1)+2)) * (*(*(M2+2)+0));
	*(*(M3+1)+1)	= (*(*(M1+1)+0)) * (*(*(M2+0)+1)) +(*(*(M1+1)+1)) * (*(*(M2+1)+1)) + (*(*(M1+1)+2)) * (*(*(M2+2)+1));
	*(*(M3+1)+2)	= (*(*(M1+1)+0)) * (*(*(M2+0)+2)) +(*(*(M1+1)+1)) * (*(*(M2+1)+2)) + (*(*(M1+1)+2)) * (*(*(M2+2)+2));
	*(*(M3+2)+0)	= (*(*(M1+2)+0)) * (*(*(M2+0)+0)) +(*(*(M1+2)+1)) * (*(*(M2+1)+0)) + (*(*(M1+2)+2)) * (*(*(M2+2)+0));
	*(*(M3+2)+1)	= (*(*(M1+2)+0)) * (*(*(M2+0)+1)) +(*(*(M1+2)+1)) * (*(*(M2+1)+1)) + (*(*(M1+2)+2)) * (*(*(M2+2)+1));
	*(*(M3+2)+2)	= (*(*(M1+2)+0)) * (*(*(M2+0)+2)) +(*(*(M1+2)+1)) * (*(*(M2+1)+2)) + (*(*(M1+2)+2)) * (*(*(M2+2)+2));
}

void LinearTransform(double **M,double *V,double *Vp,long dim)
{
	long i,j;
	for(i=0;i<dim;i++)
	{
		double tmp=0;
		for(j=0;j<dim;j++)
		{
			tmp+=M[i][j]*V[j];
		}
		Vp[i]=tmp;
	}
}
void MatrixMult(double **M1,double **M2,double **M,long dim)
{
	long i,j,k;
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
void MatrixTranspose(double **M,double **N,long dim)
{
	long i,j;
	for(i=0;i<dim;i++)
	{
		for(j=0;j<dim;j++)
		{
			N[i][j]=M[j][i];
		}
	}
}




// Get local frame by local magnetic field
void LocalFrame_B(double *B, double *e1,double *e2, double *b)
{
	double normB,normB_xy,Bx=B[0],By=B[1],Bz=B[2];
	
	normB=V3norm(B);
	normB_xy=sqrt(B[0]*B[0]+B[1]*B[1]);

	if (normB==0)
	{
		printf("In file WYLfunc.c: There is no magnetic field. Cannot set P_init_type=\"GC\"\n");
	}       
	else if (normB_xy==0)
	{
		*e1	=1;
		*(e1+1)	=0;
		*(e1+2)	=0;
		*e2	=0;
		*(e2+1)	=1;
		*(e2+2)	=0;
		*b	=0;
		*(b+1)	=0;
		*(b+2)	=1;
	}
	else
	{
		*e1 	= Bz*Bx/(normB*normB_xy);
		*(e1+1)	= Bz*By/(normB*normB_xy);
    		*(e1+2)	= -1*normB_xy/normB;
    		*e2	= -1*By/normB_xy;
    		*(e2+1)	= Bx/normB_xy;
    		*(e2+2)	= 0;
    		*b	= Bx/normB;
    		*(b+1)	= By/normB;
    		*(b+2)	= Bz/normB;
	}
		
}

// Tranform (p_\parallel,p_\perp,phase) to (px,py,pz)
void p_GC_TO_p_CART(double *p0,double *B,double *p)
{
	double p_para=p0[0],p_perp=p0[1],p_phase=p0[2];
	double normB,temp[3];
	double e1[3],e2[3],b[3],Vp_perp[3];

	LocalFrame_B(B,e1,e2,b);
	V3Smult(cos(p_phase),e1,Vp_perp);
	V3Smult(sin(p_phase),e2,temp);
	V3plus(temp,Vp_perp,Vp_perp);
	V3Smult(p_perp,Vp_perp,Vp_perp);

	V3Smult(p_para,b,temp);
	V3plus(temp,Vp_perp,p);	

}

// string 
int SetNum4String(char *str, long num)
{
	char tmp[30];
	sprintf(tmp,"%ld",num);
	strcat(str,tmp);
	return 0;
}	

//Generate random number
double GenRandNum_Uniform(double down,double up){
	return down + (up-down) * ((rand()+0.5) / ((double) RAND_MAX));
}

double func_normal(double x, double miu,double sigma)
{
    return 1.0/sqrt(2*M_PI)/sigma*exp(-1*(x-miu)*(x-miu)/(2*sigma*sigma));
}

double GenRandNum_Norm(double miu,double sigma){
	double x,y,dScope;
	double min=miu-5*sigma;
	double max=miu+5*sigma;
	do{
        	x=GenRandNum_Uniform(min,max);
        	y=func_normal(x,miu,sigma);
        	dScope=GenRandNum_Uniform(0.0,func_normal(miu,miu,sigma));
    	}while(dScope>y);
    	return x;
}

double func_circ(double R, double Rmax,double Rmin)
{
    return 2*R/(Rmax*Rmax-Rmin*Rmin);
}

double GenRandNum_Circ(double Rmin,double Rmax){
	double x,y,dScope;
	double min=Rmin;
	double max=Rmax;
	do{
        	x=GenRandNum_Uniform(min,max);
        	y=func_circ(x,Rmax,Rmin);
        	dScope=GenRandNum_Uniform(0.0,func_circ(Rmax,Rmax,Rmin));
    	}while(dScope>y);
    return x;
}

double func_para(double x, double m)
{
    return 4*x*(m*m-x*x)/(m*m*m*m);
}
double GenRandNum_Para(double m){
	double x,y,dScope;
	do{
        	x=GenRandNum_Uniform(0.,m);
        	y=func_para(x,m);
        	dScope=GenRandNum_Uniform(0.0,func_para(sqrt(1./3)*m,m));
    	}while(dScope>y);
    return x;
}

//
//
// Linear Algebra
double ** GenMatrixSpace(long M,long N)
{
	double *Data=(double *)malloc(sizeof(double)*(M*N));	
	double **A=(double **)malloc(sizeof(double *)*M);	
	long i;

	for(i=0;i<M;i++)
	{
		long offset=i*N;
		A[i]=Data+offset;
	}
	return A;
}

int EraseMatrixSpace(double **A)
{
//		printf("fuck_free_f_0\n");
	//	printf("A[0][0]=%e\n",A[0][0]);
	free(*A);
//		printf("fuck_free_f_1\n");
	free(A);
	return 0;
}
double *** GenTensorSpace3(long Dim)
{
	double *Data=(double *)malloc(sizeof(double )*((long)pow(Dim,3)));	
	double **Ajk=(double **)malloc(sizeof(double *)*Dim*Dim);	
	double ***A=(double ***)malloc(sizeof(double **)*Dim);	
	long i;

	for(i=0;i<Dim;i++)
	{
		long offset=Dim*i;
		A[i] = Ajk+offset;
	}
	for(i=0;i<Dim*Dim;i++)
	{
		long offset=Dim*i;
		Ajk[i]=Data+offset;
	}
	return A;
}

inline void EraseTensorSpace3(double ***A)
{
	free(A[0][0]);
	free(A[0]);
	free(A);
}

//Lorentz bost matrix
void LorentzBoost(double *b, double **A,int Inv,double SmallEnough)
{
	double b1=b[0];
	double b2=b[1];
	double b3=b[2];
	double beta=V3norm(b);
	double beta2=V3INmult(b,b);
	double r=1./sqrt(1-beta2);

	
	if(beta<=SmallEnough)
	{
		int i,j;
		for(i=0;i<4;i++)
		{
			for(j=0;j<4;j++)
			{
				if(i==j)
					A[i][j]=1;	
				else
					A[i][j]=0;	
			}
		}
	}
	else if(0==Inv)
	{
		A[0][0]=r;	
		A[0][1]=-r*b1;	
		A[0][2]=-r*b2;	
		A[0][3]=-r*b3;	
	
		A[1][0]=-r*b1;	
		A[1][1]=1.+(r-1.)*b1*b1/beta2;	
		A[1][2]=(r-1.)*b1*b2/beta2;	
		A[1][3]=(r-1.)*b1*b3/beta2;	
	
		A[2][0]=-r*b2;	
		A[2][1]=(r-1.)*b1*b2/beta2;	
		A[2][2]=1.+(r-1.)*b2*b2/beta2;	
		A[2][3]=(r-1.)*b2*b3/beta2;	
	
		A[3][0]=-r*b3;	
		A[3][1]=(r-1.)*b1*b3/beta2;	
		A[3][2]=(r-1.)*b2*b3/beta2;	
		A[3][3]=1.+(r-1.)*b3*b3/beta2;	
	}
	else
	{
		A[0][0]=r;	
		A[0][1]=r*b1;	
		A[0][2]=r*b2;	
		A[0][3]=r*b3;	
	
		A[1][0]=r*b1;	
		A[1][1]=1.+(r-1.)*b1*b1/beta2;	
		A[1][2]=(r-1.)*b1*b2/beta2;	
		A[1][3]=(r-1.)*b1*b3/beta2;	
	
		A[2][0]=r*b2;	
		A[2][1]=(r-1.)*b1*b2/beta2;	
		A[2][2]=1.+(r-1.)*b2*b2/beta2;	
		A[2][3]=(r-1.)*b2*b3/beta2;	
	
		A[3][0]=r*b3;	
		A[3][1]=(r-1.)*b1*b3/beta2;	
		A[3][2]=(r-1.)*b2*b3/beta2;	
		A[3][3]=1.+(r-1.)*b3*b3/beta2;	
	}
}

//INPUT: A - array of pointers to rows of a square matrix having dimension N
//       Tol - small tolerance number to detect failure when the matrix is near degenerate
//OUTPUT: Matrix A is changed, it contains both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
//        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1 
//        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N, 
//        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S
int LUPDecompose_wiki(double **A,long N,double Tol,long *P){
	long i,j,k,imax; double maxA,*ptr,absA;
	for(i=0;i<=N;i++)
		P[i]=i; //Unit permutation matrix, P[N] initialized with N 
	for(i=0;i<N;i++){
		maxA=0.0;
		imax=i;
		for(k=i;k<N;k++)
			if((absA=fabs(A[k][i]))>maxA){ maxA=absA; imax=k; }
		if(maxA<Tol)return(0); //failure, matrix is degenerate
		if(imax!=i){
			j=P[i]; P[i]=P[imax]; P[imax]=j; //pivoting P
			ptr=A[i]; A[i]=A[imax]; A[imax]=ptr; //pivoting rows of A
			P[N]++; //counting pivots starting from N (for determinant)
		}
		for(j=i+1;j<N;j++){
			A[j][i]/=A[i][i];
			for(k=i+1;k<N;k++)
				A[j][k]-=A[j][i]*A[i][k];
		}
	}
	return 0;  //decomposition done 
}
//INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
//OUTPUT: x - solution vector of A*x=b
int LUPSolve_wiki(double **A,long *P,double *b,long N,double *x){
	long i;
	for(i=0;i<N;i++){
		x[i]=b[P[i]];
		long k;
		for(k=0;k<i;k++)
			x[i]-=A[i][k]*x[k];
	}
	for(i=N-1;i>=0;i--){
		long k;
		for(k=i+1;k<N;k++)
			x[i]-=A[i][k]*x[k];
		x[i]=x[i]/A[i][i];
	}
	return 0;
}
//INPUT: A,P filled in LUPDecompose; N - dimension
//OUTPUT: IA is the inverse of the initial matrix
int LUPInvert_wiki(double **A,long *P,long N,double **IA){
	long i,k,j;
	for(j=0;j<N;j++){
		for(i=0;i<N;i++){
			if(P[i]==j)IA[i][j]=1.0;
			else IA[i][j]=0.0;
			for(k=0;k<i;k++)
				IA[i][j]-=A[i][k]*IA[k][j];
		}
		for(i=N-1;i>=0;i--){
			for(k=i+1;k<N;k++)
				IA[i][j]-=A[i][k]*IA[k][j];
			IA[i][j]=IA[i][j]/A[i][i];
		}
	}
	return 0;
}
//INPUT: A,P filled in LUPDecompose; N - dimension. 
//OUTPUT: Function returns the determinant of the initial matrix
double LUPDeterminant_wiki(double **A,long *P,long N){
	double det=A[0][0];
	long i;
	for(i=1;i<N;i++)
		det*=A[i][i];                
	if((P[N]-N)%2==0)return(det); 
	else return(-det);
}

int LinearSolver_LU(double *x, double **A, double *b, long N)
{
	double Tol=EXTRAME_SMALL;
	long *P=(long *)malloc((N+1)*sizeof(long));
	LUPDecompose_wiki(A,N,Tol,P);
	LUPSolve_wiki(A,P,b,N,x);
	free(P);
	return 0;
}
double Determinant_LU(double **A, long N)
{
	double Tol=EXTRAME_SMALL;
	long *P=(long *)malloc((N+1)*sizeof(long));
	LUPDecompose_wiki(A,N,Tol,P);
	double det=LUPDeterminant_wiki(A,P,N);
	free(P);
	return det;
}
int Invert_LU(double ** InvA,double **A, long N)
{
	double Tol=EXTRAME_SMALL;
	long *P=(long *)malloc((N+1)*sizeof(long));
	LUPDecompose_wiki(A,N,Tol,P);
	LUPInvert_wiki(A,P,N,InvA);
	free(P);
	return 0;
}


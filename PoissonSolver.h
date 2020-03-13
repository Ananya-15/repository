#include <iostream>
//#include "PrintMat.h"
#include <math.h>
using namespace std;
//#include "LidDrivenCavity.h"
#include <cblas.h>

 #define F77NAME(x) x##_ //Definition for using lapack solver 
 extern "C" {
 // LAPACK routine for solving systems of linear equations
 void F77NAME(dsysv)(const char& uplo, const int& n, const int& nrhs,const double * A, const int& lda,int *ipiv, double* B, 
 const int& ldb, double * work, const int& lwork, int& info);
 }


class PoissonSolver //main class functions: Should create matrix A, should solve to find s1 - create two functions which can do this 

{ //Inherited class which is dependent on class LidDrivenCavity 

//private: 
//int Nx,Ny; 
//double dx,dy; 

  public: 
     PoissonSolver();  //Constructor and destructors defined 
    ~PoissonSolver();
    
//    void PoissonInitialise(int ny,int nx,double deltax,double deltay){
//        Nx=nx-2; 
//        Ny=ny-2;
//        dx=deltax; 
//        dy=deltay;
//    }
  public:
    double* MatrixPoisson(int Nx,int Ny, double dx, double dy){
             int nx=Nx-2; //Define own values here 
                int ny=Ny-2; 
                int var=nx*ny; //Temporary variable to obtain matrix A 
                double* A=new double[var*var]; 
                for (int i=0;i<var;i++){
                   A[i*var+i]=(2/(dx*dx))+(2/(dy*dy)); //Diagonal 
                   if (i<var-1 && ((i+1)%ny)!=0){
                       A[(i+1)*var+i]=-1/(dy*dy); 
                   }
                   if (i<var-(ny)){
                       A[(i+ny)*var+i]=-1/(dx*dx); 
                   }
               }
        //printmat(var,var,A); 
        
        return A; 
    }
    
    double* matrixSolve(double* A, int var, double* v1, double* s1){ //Inputs A, order of matrix A var, v1 and s1. Outputs altered s1 
        double wkopt; 
               int lwork=-1; 
               int info=0; 
               int nrhs=var; 
               int* ipiv=new int[var];
               
               //Will overwrite w1 //Copy and change it for now 
               double* varmat=new double[var]; 
               cblas_dcopy(var,v1,1,varmat,1); //Will change so that varmat=w1 
               F77NAME(dsysv)('U',var,nrhs,A,var,ipiv,varmat,var,&wkopt,lwork,info); 
               
               lwork=(int)wkopt; 
               double* work = new double[lwork]; 
               F77NAME(dsysv)('U',var,nrhs,A,var,ipiv,varmat,var,work,lwork,info); //Will produce psi values in varmat 
               
               //Copy varmat back to psi1 
               cblas_dcopy(var,varmat,1,s1,1); 
               
               return s1; 
    }
   
   
   //Apply lapack, solve to obtain inner streamfunction values. Return streamfunctions values
   
   
   
    
};
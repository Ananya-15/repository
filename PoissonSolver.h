#include <iostream>
using namespace std;
//#include "LidDrivenCavity.h"


class PoissonSolver 

{ //Inherited class which is dependent on class LidDrivenCavity 
  public: 
     PoissonSolver();  //Constructor and destructors defined 
    ~PoissonSolver();
    
   //Obtain Matrix A, should be banded stored //Ensure storage occurs in column major form 
       double* BMatrixPoisson(int Nx, int Ny, double dx, double dy){
       Nx=Nx-2; 
       Ny=Ny-2; 
       int var=Nx*Ny; 
       double* A=new double[var*var]; 
       for (int i=0;i<var;i++){
           A[i*var+i]=(2/(dx*dx))+(2/(dy*dy));
//           if (i<var && ((i+1)%Ny)!=0){
//               A((i+1)*var+i)=
           }
           
           return A; 
 
       }
   }; 
   
   
   //Apply lapack, solve to obtain inner streamfunction values. Return streamfunctions values
   
   
   
    

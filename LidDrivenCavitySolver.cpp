#include <iostream>
#include <iomanip> 
using namespace std;

#include "LidDrivenCavity.h"
//#include "PrintMat.h"
//#include "PoissonSolver.h"


int main(int argc, char **argv)
{
    // Create a new instance of the LidDrivenCavity class
   
    
    string inputs[(argc-1)/2]={"Lx","Ly","Nx","Ny","Px","Py","dt","T","Re"};

    
    int n=argc/2; //Record no. of inputs 
    string *val=new string [n]; //Create array values, storing all the input values, to be passed as inputs to class function 
    for (int i=0;i<n;i++){
        val[i]=argv[2*(i+1)]; 
    }
     LidDrivenCavity* solver = new LidDrivenCavity();
     solver->Initialise(val); 
     

     //Solver tests 
     //Check A matrix is printe corectly: To be used in PoissonSolver class 
     //Only adding values to lower triangle of A, important to specify when using Lapack or scalapack 
//       int Nx=5-2; //Define own values here 
//       int Ny=5-2; 
//       double dx=0.1667; 
//       double dy=0.1667; 
//       int var=Nx*Ny; 
//       double* A=new double[var*var]; 
//       for (int i=0;i<var;i++){
//           A[i*var+i]=(2/(dx*dx))+(2/(dy*dy));
//           if (i<var-1 && ((i+1)%Ny)!=0){
//               A[(i+1)*var+i]=-1/(dy*dy); 
//           }
//           if (i<var-(Ny)){
//               A[(i+Ny)*var+i]=-1/(dx*dx); 
//           }
//        printmat(var,var,A); 
           

     
     


    // Configure the solver here...
    // ...
//    solver->Initialise();
//
//    // Run the solver
    solver->Integrate();

	return 0;
}


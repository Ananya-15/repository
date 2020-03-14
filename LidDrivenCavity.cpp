
#include "LidDrivenCavity.h"
#include "PoissonSolver.h"
#include <iostream>
#include <math.h>
#include "PrintMat.h"
#include <cblas.h>

 #define F77NAME(x) x##_ //Definition for using lapack solver 
 extern "C" {
 // LAPACK routine for solving systems of linear equations
 void F77NAME(dsysv)(const char& uplo, const int& n, const int& nrhs,const double * A, const int& lda,int *ipiv, double* B, 
 const int& ldb, double * work, const int& lwork, int& info);
 }
 
 //Create default destructor and constructor 
LidDrivenCavity::LidDrivenCavity()=default; //Constructor 


LidDrivenCavity::~LidDrivenCavity()=default;  //Destructor 


void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
   Lx=xlen; 
   Ly=ylen; 
}

void LidDrivenCavity::SetGridSize(int nx, int ny)
{
    Nx=nx; 
    Ny=ny; 
}

void LidDrivenCavity::SetTimeStep(double deltat)
{
    dt=deltat; 
}

void LidDrivenCavity::SetFinalTime(double finalt)
{
    T=finalt; 
}

void LidDrivenCavity::SetReynoldsNumber(double re)
{
    Re=re;
}

void LidDrivenCavity::Initialise(string* val) //Initialise all variables //Initialise all arrays required for the operation
{    
    
       Lx=stod(val[0]); 
       Ly=stod(val[1]); 
       Nx=stoi(val[2]); 
       Ny=stoi(val[3]); 
       Px=stoi(val[4]); 
       Py=stoi(val[5]);
       dt=stod(val[6]);
       T=stoi(val[7]);  
       Re=stod(val[8]); 

      v=new double[Nx*Ny]; //Initialize vorticity values 
      vnew=new double[Nx*Ny];//Voriticity at new timestep value 
      s=new double[Nx*Ny]; 
      v1=new double[(Nx-2)*(Ny-2)]; //Initialize inner vorticity and streamfunction values 
      s1=new double[(Nx-2)*(Ny-2)]; 
      nx=Nx-2; 
      ny=Ny-2; 
      var=nx*ny; //
      A=new double[var*var]; 
      
      for (int i=0;i<var;i++){
          for (int j=0;j<var;j++){
              A[i*var+j]=0; 
          }
      }
        dx=Lx/(Nx-1.0); 
        dy=Lx/(Ny-1.0);
        U=1.0; 
       
       
       if (dt>=Re*dx*dy/4){
           cout << "Chosen value of dt too large" << endl; 
//           break; 
       }

}

 //Additional defined class functions 
 
// 
// double* LidDrivenCavity::PoissonMatrix(int nx,int ny, int var, double dx, double dy){ //Gets the function to return a pointer of size A
//           //Create matrix A //get ad array length error for some reason :( //Temporary variable to obtain matrix A 
//
//        for (int i=0;i<var;i++){
//            A[i*var+i]=(2/(dx*dx))+(2/(dy*dy));
//           if (i<var-1 && ((i+1)%ny)!=0){
//               A[(i+1)*var+i]=-1/(dy*dy); 
//           }
//           else {
//               A[(i+1)*var+i]=0; 
//           }
//           if (i<var-(ny)){
//              A[(i+ny)*var+i]=-1/(dx*dx); 
//           }
//           else {
//               A[(i+ny)*var+i]=0;
//               
//           }
//         }
//     return A; 
// }
 
 double* LidDrivenCavity::SolveMatrix(double* A, double* w1, double* s1,int var){ //Returns s1 
 
 
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
           
           delete[] varmat; 
           delete [] v1; 
 }

 void LidDrivenCavity::BoundaryConditions(int Nx, int Ny, double* v, double dx, double dy, double dt, double U){//){
      for (int i=0;i<Nx;i++){
                v[(i+1)*(Ny-1)+i]=(2/pow(dy,2))*(s[(i+1)*(Ny-1)]-s[(i+1)*(Ny-1)-1])-(2*U/dy); //Top
                 v[i*Ny]=(2/pow(dy,2))*(s[i*Ny+1]-s[i*Ny+2]); //Bottom 
            }
       for (int i=0;i<Ny;i++){
        v[i]=(2/(pow(dx,2)))*(s[i]-s[i+Ny]); //left 
        v[Ny*(Nx-1)+i]=(2/pow(dx,2))*(s[(Nx-1)*Ny+i]-s[(Nx-2)*Ny+i]); //Right 
        }
 }

void LidDrivenCavity::InnerVorticity(double* v, double* s, int Nx, int Ny, double dx, double dy){
      for (int i=1;i<Nx-1;i++){
                for (int j=1;j<Ny-1;j++){
                    v[i*Ny+j]=-(s[(i*Ny)+j-1]-2*s[i*Ny+j]+s[(i*Ny)+j+1])/(pow(dy,2))-(s[(i-1)*Ny+j]-2*s[i*Ny+j]+s[(i+1)*Ny+j])/(pow(dx,2)); 
                }
      }
 }
 
 void LidDrivenCavity::NextInnerVorticity(double* v, double* s, int Nx, int Ny, double dx, double dy, double dt, double Re){
        for (int i=1;i<Nx-1;i++){
            for (int j=1;j<Ny-1;j++){ //Create temporary variables to ensure ease of reading 
                double temp1=(s[(i)*Ny+j+1]-s[i*Ny+j-1])*(v[(i+1)*Ny+j]-v[(i-1)*Ny+j])/(4*dy*dx); 
                double temp2=(s[(i+1)*Ny+j]-s[(i-1)*Ny+j])*(v[i*Ny+j+1]-v[i*Ny+j-1])/(4*dy*dx); 
                double temp3=(1/Re)*((v[(i*Ny)+j-1]-2*v[i*Ny+j]+v[(i*Ny)+j+1])/(pow(dy,2))+(v[(i-1)*Ny+j]-2*v[i*Ny+j]+v[(i+1)*Ny+j])/(pow(dx,2))); 
                vnew[i*Ny+j]=v[i*Ny+j]+dt*(temp2-temp1+temp3); 
                
            }
        }

 }
 


void LidDrivenCavity::RecoverInnerVorticity(double*vnew, double*v1, int Nx, int Ny){
     for(int i=1; i<Nx-1;i++){
                for (int j=1;j<Ny-1;j++){
                    v1[(i-1)*(Ny-2)+j-1]=vnew[(i)*Ny+j]; 
                }
            }
            
}


 void LidDrivenCavity::CopyVorticity(double*arr1,double* arr2, int var){ //For first one vnew=v, for second one, v=vnew 
      cblas_dcopy(var,arr1,1,arr2,1); 
 }


void LidDrivenCavity::UpdateInnerStream(double*s, double*s1, int Nx, int Ny){
            for(int i=1; i<Nx-1;i++){
                for (int j=1;j<Ny-1;j++){
                    s[(i)*Ny+j]=s1[(i-1)*(Ny-2)+j-1]; 
                }
               }

}



void LidDrivenCavity::Integrate()
{    
  PoissonSolver* psolver=new PoissonSolver(); //Create new instance psolver of class Poisson Solver to implement solver calculations  

  psolver->Initialise(Nx,Ny,dx,dy); //Initialises values within PoissonSolver class 
  double t=0; //First time step value  
  
    while (t<T){
            cout << "Time step is: " << t << endl << endl; 
            
             LidDrivenCavity::BoundaryConditions(Nx,Ny,v,dx,dy,dt,U); //Update with BCs //
             
             
            
            cout << "Updated vorticity boundar conditions: " << endl; 
            cout << endl << endl; 
            printmat(Ny,Nx,v); 
            
            
             LidDrivenCavity::InnerVorticity(v,s,Nx,Ny,dx,dy); //Inner vorticity values at current timestep 
             
             LidDrivenCavity::CopyVorticity(v,vnew,Nx*Ny); 
            cout << "Inner voritcity values at current timestep" << endl;    
            cout << endl << endl; 
            printmat(Ny,Nx,v);            
  
            LidDrivenCavity::NextInnerVorticity(v,s,Nx,Ny,dx,dy,dt,Re);  //Inner vorticity values at next timestep, updates vnew
        
            
            cout << "Updated vorticity at new timestep: " << endl; 
            printmat(Ny,Nx,vnew); 
            cout << endl << endl; 
            
            LidDrivenCavity::RecoverInnerVorticity(vnew,v1,Nx,Ny); //Extract to obain inner voriticity values and save in matrix v1 
            
            LidDrivenCavity::CopyVorticity(vnew,v,Nx*Ny); //Update vorticity so that v=vnew 
            
            A=psolver->MatrixPoisson(); //Call function to obtain matrix A
               
            cout << "Works here 1" << endl;
            
            cout << "Updated vorticity at new timestep: " << endl; 
            printmat(Ny,Nx,v); 
            cout << endl << endl; 
            
             s1=psolver->matrixSolve(A,var,v,s1); //Time to implement lapack 
            
            //s1=LidDrivenCavity::SolveMatrix(A,v1,s1,var); 
              cout << "Works here 2" << endl; 
               
            LidDrivenCavity::UpdateInnerStream(s,s1,Nx,Ny);  //Need to input s1 back in streamfunction 

              cout << "Works here 3" << endl; 
          t+=dt; 
      }
    cout << "Final streamfunction values: " << endl; 
    printmat(Nx,Ny,v); //Check if final values make sense 

  //delete [] A; 
  delete [] v1; 
  delete [] s1; 
  delete [] v; 
  delete [] s; 
     
}
      
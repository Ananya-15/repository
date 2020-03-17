
#include "LidDrivenCavity.h" //LidDrivenCavity ut with banded 
//#include "PoissonSolver.h"
#include <iostream>
#include <math.h>
#include "PrintMat.h"
#include <cblas.h>

 #define F77NAME(x) x##_ //Definition for using lapack solver 
 extern "C" {
 // LAPACK routine for solving systems of linear equations
 void F77NAME(dgbtrf)(const int& n, const int& m, const int& kl, const int& ku, const double * Ab, const int& ldab,int *ipiv, int& info);
 }
  #define F77NAME(x) x##_ //Definition for using lapack solver 
 extern "C" {
 // LAPACK routine for solving systems of linear equations
 void F77NAME(dgbtrs)(const char& trans, const int& n, const int& kl, const int& ku, const int& nrhs,const double * Ab, const int& lda,int *ipiv, double* B, 
 const int& ldb, int& info);
 }
 
 LidDrivenCavity::LidDrivenCavity()=default; //Constructor 
//{
//    
//}


LidDrivenCavity::~LidDrivenCavity() //Destructor 
{
}
 
 double* LidDrivenCavity::PoissonMatrix(double dx, double dy){
           //Create matrix A //get ad array length error for some reason :( 
        int Ku=ny; //Comes from problem definition 
        int Kl=ny; 
        int rows=3*Kl+1;
        var=nx*ny; 
        int mid=(Ku+Kl+1)/2+Kl; //Position of mid diagonal 
        
        for (int i=0;i<var;i++){
            for (int j=0;j<rows;j++){
                A[i*rows+j]=0; //Initialise and set everything as zero 
            }
        }
        
        for (int i=0;i<var;i++){ //Going across each row 
           
                 A[i*rows+mid]=(2/(dx*dx))+(2/(dy*dy)); //Diagoanl value 
                
                   if (i>ny-1){
                        A[i*rows+Kl]=-1/(dx*dx);
                   }
                   if (i>0 && ((i)%ny)!=0){
                       A[i*rows+mid-1]=-1/(dy*dy); 
                   }
//
                    if (i<var-1 && ((i+1)%ny)!=0){
                       A[i*rows+mid+1]=-1/(dy*dy); 
                       }
                  if (i<var-ny){
                        A[i*rows+rows-1]=-1/(dx*dx);
                   }
        }
         printmat(rows,var,A);
         cout << endl << endl; 
         int info=0; 
         int* ipiv=new int[var];
         //cout << "After factorisation: " << endl; 
         F77NAME(dgbtrf)(var,var,Ku,Kl,A,rows,ipiv,info); 
        // printmat(rows,var,A); 
     return A; 
 }
 
 double* LidDrivenCavity::SolveMatrix(double* A, double* w1, double* s1,int var,int Nx,int Ny){ //Returns s1 
 
         int info=0; 
         int* ipiv=new int[var];
           
           int nrhs=1; 
           int ldab=3*ny+1; 
           
           double* varmat=new double[var]; 
           cblas_dcopy(var,v1,1,varmat,1);
           
          // printmat(Nx,Ny,varmat); 
           cout << endl << endl; 
           
           F77NAME(dgbtrs)('N',var,ny,ny,nrhs,A,ldab,ipiv,varmat,var,info); 
           
           ///printmat(Nx,Ny,varmat); 
           
           cblas_dcopy(var,varmat,1,s1,1); 

//           
           return s1; 
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
       T=stod(val[7]);  
       Re=stod(val[8]); 

      v=new double[Nx*Ny]; //Initialize vorticity values 
      vnew=new double[Nx*Ny];
      s=new double[Nx*Ny]; 
      nx=Nx-2; 
      ny=Ny-2; 
      var=nx*ny; //
      v1=new double[var]; //Initialize inner vorticity and streamfunction values 
      s1=new double[var]; 
      A=new double[(3*ny+1)*(var)]; 
      
      
        dx=Lx/(Nx-1.0); 
        dy=Lx/(Ny-1.0);
        U=1.0; 
       
       
       if (dt>=Re*dx*dy/4){
          throw logic_error("Chosen timestep value too small");
       }

}

 //Additional defined class functions 

 void LidDrivenCavity::BoundaryConditions(int Nx, int Ny, double* v, double dx, double dy, double dt, double U){//){
      for (int i=0;i<Nx;i++){
                v[(i+1)*(Ny-1)+i]=(2/pow(dy,2))*(s[(i+1)*(Ny-1)]-s[(i+1)*(Ny-1)-1])-(2*U/dy); //Top
                 v[i*Ny]=(2/pow(dy,2))*(s[i*Ny]-s[i*Ny+1]); //Bottom 
            }
       for (int i=0;i<Ny;i++){
        v[i]=(2/(pow(dx,2)))*(s[i]-s[i+Ny]); //left 
        v[Ny*(Nx-1)+i]=(2/pow(dx,2))*(s[(Nx-1)*Ny+i]-s[(Nx-2)*Ny+i]); //Right 
        }
 }

void LidDrivenCavity::InnerVorticity(double* v, double* s, int Nx, int Ny, double dx, double dy){
      for (int i=1;i<Nx-1;i++){
                for (int j=1;j<Ny-1;j++){
                    v[i*Ny+j]=-((s[(i*Ny)+j-1]-2*s[i*Ny+j]+s[(i*Ny)+j+1])/pow(dy,2))-((s[(i-1)*Ny+j]-2*s[i*Ny+j]+s[(i+1)*Ny+j])/pow(dx,2)); 
                }
      }
 }
 
 void LidDrivenCavity::NextInnerVorticity(double* v, double* s, int Nx, int Ny, double dx, double dy, double dt, double Re){
        for (int i=1;i<Nx-1;i++){
            for (int j=1;j<Ny-1;j++){ //Create temporary variables to ensure ease of reading 
                double temp1=(s[i*Ny+j+1]-s[i*Ny+j-1])*(v[(i+1)*Ny+j]-v[(i-1)*Ny+j])/(4*dy*dx); 
                double temp2=(s[(i+1)*Ny+j]-s[(i-1)*Ny+j])*(v[i*Ny+j+1]-v[i*Ny+j-1])/(4*dy*dx); 
                double temp3=(1/Re)*(((v[(i*Ny)+j-1]-2*v[i*Ny+j]+v[(i*Ny)+j+1])/pow(dy,2))+((v[(i-1)*Ny+j]-2*v[i*Ny+j]+v[(i+1)*Ny+j])/pow(dx,2))); 
                vnew[i*Ny+j]=v[i*Ny+j]+(dt*(temp2-temp1+temp3)); 
                
            }
        }

 }

void LidDrivenCavity::RecoverInnerVorticity(double*v, double*v1, int Nx, int Ny){
     for(int i=1; i<Nx-1;i++){
                for (int j=1;j<Ny-1;j++){
                    v1[(i-1)*ny+j-1]=v[(i)*Ny+j]; 
                }
            }
            
}

 void LidDrivenCavity::CopyVorticity(double*arr1,double* arr2, int var){ //For first one vnew=v, for second one, v=vnew 
      cblas_dcopy(var,arr1,1,arr2,1); 
 }

void LidDrivenCavity::UpdateInnerStream(double*s, double*s1, int Nx, int Ny){
            for(int i=1; i<Nx-1;i++){
                for (int j=1;j<Ny-1;j++){
                    s[(i)*Ny+j]=s1[(i-1)*ny+j-1]; 
                }
               }

}



void LidDrivenCavity::Integrate()
{    
  //PoissonSolver* psolver=new PoissonSolver(); //Create new instance psolver of class Poisson Solver to implement solver calculations  

  //psolver->Initialise(Nx,Ny,dx,dy); //Initialises values within PoissonSolver class 
  A=LidDrivenCavity::PoissonMatrix(dx,dy); //Call function to obtain matrix A
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
  
            LidDrivenCavity::NextInnerVorticity(v,s,Nx,Ny,dx,dy,dt,Re);  //Inner vorticity values at next timestep 
        
            
            cout << "Updated vorticity at new timestep: " << endl; 
            printmat(Ny,Nx,v); 
            cout << endl << endl; 
            
            LidDrivenCavity::CopyVorticity(vnew,v,Nx*Ny); 
            
            LidDrivenCavity::RecoverInnerVorticity(v,v1,Nx,Ny); //Extract to obain inner voriticity values and save in matrix v1 
            
               
           // cout << "Works here 1" << endl;
            
             s1=LidDrivenCavity::SolveMatrix(A,v1,s1,var,Nx,Ny);  //Time to implement lapack 
            
            //s1=LidDrivenCavity::SolveMatrix(A,v1,s1,var); 
             // cout << "Works here 2" << endl; 
               
            LidDrivenCavity::UpdateInnerStream(s,s1,Nx,Ny);  //Need to input s1 back in streamfunction 

             // cout << "Works here 3" << endl; 
             
             cout << "Stream function values at next timestep: " << endl ; 
             printmat(Ny,Nx,s); 
             cout << endl << endl; 
          t+=dt; 
      }
    //cout << "Final streamfunction values: " << endl; 
    printmat(Ny,Nx,s); //Check if final values make sense 

  delete [] A; 
  delete [] v1; 
  delete [] s1; 
  delete [] v; 
  delete [] s; 
     
}
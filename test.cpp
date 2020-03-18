
#include "LidDrivenCavity.h" //LidDrivenCavity ut with banded 
//#include "PoissonSolver.h"
#include <iostream>
#include <math.h>
#include "PrintMat.h"
#include <cblas.h>
#include <iomanip> 
#include <fstream>


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
 
 void LidDrivenCavity::PoissonMatrix(double dx, double dy){
           //Create matrix A //get ad array length error for some reason :( 
        int Ku=ny; //Comes from problem definition 
        int Kl=ny; 
        int rows=3*Kl+1;
        var=nx*ny; 
        int mid=(Ku+Kl+1)/2+Kl; //Position of mid diagonal 
        
        for (int i=0;i<var*rows;i++){
                A[i]=0; //Initialise and set everything as zero 
            }
        
        for (int i=0;i<var;i++){ //Going across each column
           
                 A[i*rows+mid]=(2/(dx*dx))+(2/(dy*dy)); //Diagoanl value 
                
                 for (int j=0;j<rows-1;j++){
                       if (j<Kl){
                           A[i*rows+j]='*'; 
                       }
                   }
                   if (i>ny-1){
                        A[i*rows+Kl]=-1/(dx*dx);
                   }
                   else {
                       A[i*rows+Kl]='*'; 
                   }
                   if (i>0 && ((i)%ny)!=0){
                       A[i*rows+mid-1]=-1/(dy*dy); 
                   }
//
                    if (i<var-1 && ((i+1)%ny)!=0){
                       A[i*rows+mid+1]=-1/(dy*dy); 
                       }
//              
                  if (i<var-ny){
                        A[i*rows+rows-1]=-1/(dx*dx);
                   }
                   else {
                       A[i*rows+rows-1]='*'; 
                   }
        }
                   
        
         printmat(rows,var,A);
         cout << endl << endl; 
         int info; 
      
         //cout << "After factorisation: " << endl; 
         F77NAME(dgbtrf)(var,var,Ku,Kl,A,rows,ipiv,info); 
         if (info!=0){
             throw logic_error("Info is not equal to zero!");
         }
//         F77NAME(dgbtrf)(rows,var,Ku,Kl,A,rows,ipiv,info); 
        // printmat(rows,var,A); 
     //return A; 
 }
 
 void LidDrivenCavity::SolveMatrix(double* A, double* w1, double* s1,int var,int Nx,int Ny,int* ipiv){ //Returns s1 
 
         int info; 
       ///  int* ipiv=new int[var];
           
           int nrhs=1; 
           int ldab=3*ny+1; 
           
           double* varmat=new double[var]; 
           cblas_dcopy(var,v1,1,varmat,1);
           
          // printmat(Nx,Ny,varmat); 
           cout << endl << endl; 
           
           F77NAME(dgbtrs)('N',var,ny,ny,nrhs,A,ldab,ipiv,varmat,var,info); 
           
          if (info!=0){
             throw logic_error("Info is not equal to zero!");
         }
           
           ///printmat(Nx,Ny,varmat); 
           
           cblas_dcopy(var,varmat,1,s1,1); 

//           
          // return s1; 
           delete[] varmat;
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
      ipiv=new int[var];
      v1=new double[var]; //Initialize inner vorticity and streamfunction values 
      s1=new double[var]; 
      A=new double[(3*ny+1)*(var)]; 
      
      for (int i=0;i<Nx*Ny;i++){
          v[i]=0; 
          s[i]=0; 
      }
      for (int i=0;i<var;i++){
          v1[i]=0; 
          s1[i]=0; 
      }
      
        dx=Lx/(Nx-1.0); 
        dy=Ly/(Ny-1.0);
        U=1.0; 
       
       
       if (dt>=Re*dx*dy/4){
          throw logic_error("Chosen timestep value too large");
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
                //v[i*Ny+j]=v[i*Ny+j]+(dt*(temp2-temp1+temp3)); 
                
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

 void LidDrivenCavity::CopyVorticity(double*arr1,double* arr2, int val){ //For first one vnew=v, for second one, v=vnew 
      cblas_dcopy(val,arr1,1,arr2,1); 
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
  //A=
  LidDrivenCavity::PoissonMatrix(dx,dy); //Call function to obtain matrix A
  double t=0; //First time step value  
  
    while (t<=T){
            cout << "Time step is: " << t << endl << endl; 
            
             LidDrivenCavity::BoundaryConditions(Nx,Ny,v,dx,dy,dt,U); //Update with BCs //
             
          //   LidDrivenCavity::CopyVorticity(v,vnew,Nx*Ny); //Shouldn't be doing this, vnew shouldn't have these boundary conditions 
            
//            cout << "Updated vorticity boundary conditions: " << endl; 
//            cout << endl << endl; 
//            printmat(Ny,Nx,v); 
             LidDrivenCavity::InnerVorticity(v,s,Nx,Ny,dx,dy); //Inner vorticity values at current timestep 
             
            
//            cout << "Inner voritcity values at current timestep" << endl;    
//            cout << endl << endl; 
//            printmat(Ny,Nx,v);            
  
            LidDrivenCavity::NextInnerVorticity(v,s,Nx,Ny,dx,dy,dt,Re);  //Inner vorticity values at next timestep 
        
        
             for(int i=1; i<Nx-1;i++){
                for (int j=1;j<Ny-1;j++){
                    v[(i)*Ny+j]=vnew[(i)*Ny+j]; 
                }
               }

//            cout << "Updated voriticty with Vnew" << endl; 
//            cout << endl << endl; 
//            cout << "Updated vorticity at new timestep: " << endl; 
//            printmat(Ny,Nx,v); //Should have bcs from previous timestep 
//            cout << endl << endl; 
//            
            
           // LidDrivenCavity::CopyVorticity(vnew,v,Nx*Ny); //Ensures v=vnew for the next timestep 
             
            
            
            
            LidDrivenCavity::RecoverInnerVorticity(v,v1,Nx,Ny); //Extract to obain inner voriticity values and save in matrix v1 
            
               
           // cout << "Works here 1" << endl;
            
             //s1=
             LidDrivenCavity::SolveMatrix(A,v1,s1,var,Nx,Ny,ipiv);  //Time to implement lapack 
            
            //s1=LidDrivenCavity::SolveMatrix(A,v1,s1,var); 
             // cout << "Works here 2" << endl; 
               
            LidDrivenCavity::UpdateInnerStream(s,s1,Nx,Ny);  //Need to input s1 back in streamfunction 

             // cout << "Works here 3" << endl; 
             
//             cout << "Stream function values at next timestep: " << endl ; 
//             printmat(Ny,Nx,s); 
//             cout << endl << endl; 
         t+=dt; 
      }
//    cout << "Final streamfunction values: " << endl; 
//    printmat(Ny,Nx,s); //Check if final values make sense 
   
  ofstream vOut("StreamFunctionVal.txt",ios::out); 
   for (int i=0;i<Nx*Ny;i++){
    vOut.precision(5); //n is no. of significant figures 
    vOut << setw(7) << s[i] << endl; //Sets space between values 
   // vOut < setfill('.') << vairbale //Generally, leave empty!
   }

    vOut.close(); 

  delete [] A; 
  delete [] v1; 
  delete [] s1; 
  delete [] v; 
  delete [] s; 
     
}
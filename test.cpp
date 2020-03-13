
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
 
 
 double* LidDrivenCavity::PoissonMatrix(int Nx, int Ny, double dx, double dy){
           //Create matrix A //get ad array length error for some reason :( 
        int nx=Nx-2; //lower 
        int ny=Ny-2; 
        int Ku=ny; //Comes from problem definition 
        int Kl=ny; 
        int rows=Ku+Kl+Kl+1;
        int var=nx*ny; 
        int mid=(Ku+Kl+1)/2+Kl; //Position of mid diagonal 

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
//         }
        }
        
              //printmat(rows,var,A); 
              cout << endl << endl; 
         
        
         int info=0; 
         int* ipiv=new int[var];
         F77NAME(dgbtrf)(var,var,Ku,Kl,A,rows,ipiv,info); 
      // printmat(rows,var,A); 
     return A; 
     delete[] A; 
 }
 
 double* LidDrivenCavity::SolveMatrix(double* A, double* w1, double* s1,int var,int Nx,int Ny){ //Returns s1 
 
           int nx=Nx-2; 
           int ny=Ny-2; 
//           double wkopt; 
//           int lwork=-1; 
         int info=0; 
         int* ipiv=new int[var];
           
           int nrhs=var; 
           int ldab=2*ny+ny+1; 
           
           double* varmat=new double[var]; 
           cblas_dcopy(var,v1,1,varmat,1);
           
          // printmat(Nx,Ny,varmat); 
           cout << endl << endl; 
           
           F77NAME(dgbtrs)('N',var,ny,ny,nrhs,A,ldab,ipiv,varmat,var,info); 
           
           ///printmat(Nx,Ny,varmat); 
           
           cblas_dcopy(var,varmat,1,s1,1); 
           
           //Will overwrite w1 //Copy and change it for now 
//           double* varmat=new double[var]; 
//           cblas_dcopy(var,v1,1,varmat,1); //Will change so that varmat=w1 
//           F77NAME(dsysv)('U',var,nrhs,A,var,ipiv,varmat,var,&wkopt,lwork,info); 
//           
//           lwork=(int)wkopt; 
//           double* work = new double[lwork]; 
//           F77NAME(dsysv)('U',var,nrhs,A,var,ipiv,varmat,var,work,lwork,info); //Will produce psi values in varmat 
//           
//           //Copy varmat back to psi1 
//           cblas_dcopy(var,varmat,1,s1,1); 
//           
           return s1; 
//           
//           delete[] s1; 
 }

 double* LidDrivenCavity::BoundaryConditions(int Nx, int Ny, double* v, double dx, double dy, double dt, double U){
      for (int i=0;i<Nx;i++){
                v[(i+1)*(Ny-1)+i]=(2/pow(dy,2))*(s[(i+1)*(Ny-1)]-s[(i+1)*(Ny-1)-1])-(2*U/dy); //Top
                 v[i*Ny]=(2/pow(dy,2))*(s[i*Ny+1]-s[i*Ny+2]); //Bottom 
            }
       for (int i=0;i<Ny;i++){
        v[i]=(2/(pow(dx,2)))*(s[i]-s[i+Ny]); //left 
        v[Ny*(Nx-1)+i]=(2/pow(dx,2))*(s[(Nx-1)*Ny+i]-s[(Nx-2)*Ny+i]); //Right 
        }
    return v; 
 }

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
     //  cout << Lx << endl; 
       Ly=stod(val[1]); 
     //  cout << Ly << endl; 
       Nx=stoi(val[2]); 
       cout << Nx << endl; 
        Ny=stoi(val[3]); 
       cout << Ny << endl; 
        Px=stoi(val[4]); 
      // cout << Px << endl; 
        Py=stoi(val[5]);
     // cout << Py << endl;  
        dt=stod(val[6]);
      // cout << dt << endl; 
        T=stoi(val[7]); 
     //  cout << T << endl; 
        Re=stod(val[8]); 
      // cout << Re << endl; 

      v=new double[Nx*Ny]; //Initialize vorticity values 
      s=new double[Nx*Ny]; 
      v1=new double[(Nx-2)*(Ny-2)]; //Initialize inner vorticity and streamfunction values 
      s1=new double[(Nx-2)*(Ny-2)]; 
      A=new double[(2*(Nx-2)+(Ny-2)+1)*(Nx-2)*(Ny-2)]; //Set default A values to be zero //5 because 2*Kl+Ku+1=5 
      
      for (int i=0;i<(Nx-2)*(Ny-2);i++){
          for (int j=0;j<(2*(Nx-2)+(Ny-2)+1);j++){
              A[i*(Nx-2)*(Nx-2)+j]=0; 
          }
      }
      
      //A=PoissonSolver::MatrixPoisson(Nx,Ny,dx,dy)
      
            cout << "Works here 1" << endl; 

}
      


void LidDrivenCavity::Integrate()
{    
       double dx=Lx/(Nx-1.0); 
       double dy=Lx/(Ny-1.0);
       double U=1.0; 
       double t=0; //First time step value  
       
       if (dt>=Re*dx*dy/4){
           cout << "Chosen value of dt too large" << endl; 
       }
       
       A=LidDrivenCavity::PoissonMatrix(Nx,Ny,dx,dy); //Return Matrix A for poisson solver 
       

        while (t<0.3){
            cout << "Time step is: " << t << endl << endl; 
            
             v=LidDrivenCavity::BoundaryConditions(Nx,Ny,v,dx,dy,dt,U); //Update with BCs 
            
//            cout << "Updated vorticity boundar conditions: " << endl; 
//            cout << endl << endl; 
//            printmat(Ny,Nx,v); 

            
            //Inner vorticity values at current timestep 
            for (int i=1;i<Nx-1;i++){
                for (int j=1;j<Ny-1;j++){
                    v[i*Ny+j]=-(s[(i*Ny)+j-1]-2*s[i*Ny+j]+s[(i*Ny)+j+1])/(pow(dy,2))-(s[(i-1)*Ny+j]-2*s[i*Ny+j]+s[(i+1)*Ny+j])/(pow(dx,2)); 
                }
            }
            
            cout << "Works fine here 1" << endl; 
//            cout << "Inner voritcity values at current timestep" << endl; 
//            cout << endl << endl; 
//            printmat(Ny,Nx,v); 
            //Inner vorticity values at next timestep 
            for (int i=1;i<Nx-1;i++){
                for (int j=1;j<Ny-1;j++){ //Create temporary variables to ensure ease of reading 
                    double temp1=(s[(i)*Ny+j+1]-s[i*Ny+j-1])*(v[(i+1)*Ny+j]-v[(i-1)*Ny+j])/(4*dy*dx); 
                    double temp2=(s[(i+1)*Ny+j]-s[(i-1)*Ny+j])*(v[i*Ny+j+1]-v[i*Ny+j-1])/(4*dy*dx); 
                    double temp3=(1/Re)*((v[(i*Ny)+j-1]-2*v[i*Ny+j]+v[(i*Ny)+j+1])/(pow(dy,2))+(v[(i-1)*Ny+j]-2*v[i*Ny+j]+v[(i+1)*Ny+j])/(pow(dx,2))); 
                    v[i*Ny+j]=v[i*Ny+j]+dt*(temp2-temp1+temp3); 
                    
                }
            }
            
//            cout << "Updated vorticity at new timestep: " << endl; 
//            printmat(Ny,Nx,v); 
//            cout << endl << endl; 
            
            
            //Extract to obain inner voriticity values and svae in matrix v1 
            for(int i=1; i<Nx-1;i++){
                for (int j=1;j<Ny-1;j++){
                    v1[(i-1)*(Ny-2)+j-1]=v[(i)*Ny+j]; 
                }
            }
            cout << "Works fine here 2" << endl; 
            
            
            //Call function from Poisson Solver class to obtain s1
            //PoissonSolver::PoissonInitialise(Ny,Nx,dx,dy);  //Can't do this for some reason, trouble understanding inherited classes 

                //printmat(var,var,A); 
               
               //Time to implement lapack 
           
               s1=LidDrivenCavity::SolveMatrix(A,v1,s1,(Nx-2)*(Ny-2),Nx,Ny);
               
//               printmat(Nx,Ny,v1); 
//               
//               cout << "Streamfunction: " << endl << endl; 
//               
//               printmat(Nx,Ny,s1);
//               
//               cout << "Works fine here 3" << endl; 
  
               
               //Need to input s1 back in streamfunction 
               
               
              for(int i=1; i<Nx-1;i++){
                for (int j=1;j<Ny-1;j++){
                    s[(i)*Ny+j]=s1[(i-1)*(Ny-2)+j-1]; 
                }
               }
         
//         cout << endl << endl; ;      
          t+=dt; 
     }
    //printmat(Nx,Ny,v); 
     
}
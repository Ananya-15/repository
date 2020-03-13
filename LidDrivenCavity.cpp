
#include "LidDrivenCavity.h"
//#include "PoissonSolver.h"
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
 
 LidDrivenCavity::LidDrivenCavity()=default; //Constructor 
//{
//    
//}


LidDrivenCavity::~LidDrivenCavity() //Destructor 
{
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
      A=new double[(Nx-2)*(Ny-2)*(Nx-2)*(Ny-2)]; 
      nx=Nx-2; 
      ny=Ny-2; 
      var=nx*ny;
      
      for (int i=0;i<var;i++){
          for (int j=0;j<var;j++){
              A[i*var+j]=0; 
          }
      }
      
      //A=PoissonSolver::MatrixPoisson(Nx,Ny,dx,dy)
      
            //cout << "Works here 1" << endl; 

}

 //Additional defined class functions 
 
 
 double* LidDrivenCavity::PoissonMatrix(int nx,int ny, int var, double dx, double dy){ //Gets the function to return a pointer of size A
           //Create matrix A //get ad array length error for some reason :( //Temporary variable to obtain matrix A 

        for (int i=0;i<var;i++){
            A[i*var+i]=(2/(dx*dx))+(2/(dy*dy));
           if (i<var-1 && ((i+1)%ny)!=0){
               A[(i+1)*var+i]=-1/(dy*dy); 
           }
           else {
               A[(i+1)*var+i]=0; 
           }
           if (i<var-(ny)){
              A[(i+ny)*var+i]=-1/(dx*dx); 
           }
           else {
               A[(i+ny)*var+i]=0;
               
           }
         }
     return A; 
     delete[] A; 
 }
 
 double* LidDrivenCavity::SolveMatrix(double* A, double* w1, double* s1,int var){ //Returns s1 
 
 
           double wkopt; 
           int lwork=-1; 
           int info=0; 
           int nrhs=var; 
           int* ipiv=new int[var];
           
           //Will overwrite w1 //Copy and change it for now 
//           double* varmat=new double[var]; 
//           cblas_dcopy(var,v1,1,varmat,1); //Will change so that varmat=w1 
           F77NAME(dsysv)('U',var,nrhs,A,var,ipiv,v1,var,&wkopt,lwork,info); 
           
           lwork=(int)wkopt; 
           double* work = new double[lwork]; 
           F77NAME(dsysv)('U',var,nrhs,A,var,ipiv,v1,var,work,lwork,info); //Will produce psi values in varmat 
           
           //Copy varmat back to psi1 
           cblas_dcopy(var,v1,1,s1,1); 
           
           return s1; 
           
           delete[] s1; 
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
    
    delete[] v; 
 }

 double* LidDrivenCavity::InnerVorticity(double* v, double* s, int Nx, int Ny, double dx, double dy){
      for (int i=1;i<Nx-1;i++){
                for (int j=1;j<Ny-1;j++){
                    v[i*Ny+j]=-(s[(i*Ny)+j-1]-2*s[i*Ny+j]+s[(i*Ny)+j+1])/(pow(dy,2))-(s[(i-1)*Ny+j]-2*s[i*Ny+j]+s[(i+1)*Ny+j])/(pow(dx,2)); 
                }
            }
            return v; 
            delete[] v; 
 }
 
 double* LidDrivenCavity::NextInnerVorticity(double* v, double* s, int Nx, int Ny, double dx, double dy, double dt, double Re){
        for (int i=1;i<Nx-1;i++){
            for (int j=1;j<Ny-1;j++){ //Create temporary variables to ensure ease of reading 
                double temp1=(s[(i)*Ny+j+1]-s[i*Ny+j-1])*(v[(i+1)*Ny+j]-v[(i-1)*Ny+j])/(4*dy*dx); 
                double temp2=(s[(i+1)*Ny+j]-s[(i-1)*Ny+j])*(v[i*Ny+j+1]-v[i*Ny+j-1])/(4*dy*dx); 
                double temp3=(1/Re)*((v[(i*Ny)+j-1]-2*v[i*Ny+j]+v[(i*Ny)+j+1])/(pow(dy,2))+(v[(i-1)*Ny+j]-2*v[i*Ny+j]+v[(i+1)*Ny+j])/(pow(dx,2))); 
                v[i*Ny+j]=v[i*Ny+j]+dt*(temp2-temp1+temp3); 
                
            }
        }
        return v;
        delete[] v; 
 }

double* LidDrivenCavity::RecoverInnerVorticity(double*v, double*v1, int Nx, int Ny){
     for(int i=1; i<Nx-1;i++){
                for (int j=1;j<Ny-1;j++){
                    v1[(i-1)*(Ny-2)+j-1]=v[(i)*Ny+j]; 
                }
            }
        return v1; 
        delete[] v1; 
            
}

double* LidDrivenCavity::UpdateInnerStream(double*s, double*s1, int Nx, int Ny){
            for(int i=1; i<Nx-1;i++){
                for (int j=1;j<Ny-1;j++){
                    s[(i)*Ny+j]=s1[(i-1)*(Ny-2)+j-1]; 
                }
               }
               return s; 
               delete[] s; 
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
       
        //Return Matrix A for poisson solver 
       //printmat((Nx-2)*(Ny-2),(Nx-2)*(Ny-2),A);
       while (t<T){
            cout << "Time step is: " << t << endl << endl; 
            
             v=LidDrivenCavity::BoundaryConditions(Nx,Ny,v,dx,dy,dt,U); //Update with BCs 
            
            cout << "Updated vorticity boundar conditions: " << endl; 
            cout << endl << endl; 
            printmat(Ny,Nx,v); 
             v=LidDrivenCavity::InnerVorticity(v,s,Nx,Ny,dx,dy); 
             
            cout << "Inner voritcity values at current timestep" << endl; 
            cout << endl << endl; 
            printmat(Ny,Nx,v); 
            //Inner vorticity values at next timestep 
//            
            //Inner vorticity values at current timestep 
            v=LidDrivenCavity::NextInnerVorticity(v,s,Nx,Ny,dx,dy,dt,Re); 
        
            
            cout << "Updated vorticity at new timestep: " << endl; 
            printmat(Ny,Nx,v); 
            cout << endl << endl; 
            
            
            //Extract to obain inner voriticity values and save in matrix v1 
            v1=LidDrivenCavity::RecoverInnerVorticity(v,v1,Nx,Ny); 
            
        
            //Call function from Poisson Solver class to obtain s1
            //PoissonSolver::PoissonInitialise(Ny,Nx,dx,dy);  //Can't do this for some reason, trouble understanding inherited classes 

                //printmat(var,var,A); 
               
               //Time to implement lapack 
            A=LidDrivenCavity::PoissonMatrix(nx,ny,var,dx,dy);
            
//            cout << "Matrix A is: " << endl; 
//            printmat(var,var,A); 
//            cout << endl << endl; 
           
            s1=LidDrivenCavity::SolveMatrix(A,v1,s1,var); 
               
               //Need to input s1 back in streamfunction 
               
            s=LidDrivenCavity::UpdateInnerStream(s,s1,Nx,Ny); 

         
//         cout << endl << endl; ;      
          t+=dt; 
          
         // delete[] A; 
      }
    printmat(Nx,Ny,s); //Check if final values make sense 

     
}
      
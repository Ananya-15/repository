
#include "LidDrivenCavity.h"
//#include "PoissonSolver.h"
#include <iostream>
#include <math.h>
#include "PrintMat.h"


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
       cout << Lx << endl; 
       Ly=stod(val[1]); 
       cout << Ly << endl; 
       Nx=stoi(val[2]); 
       cout << Nx << endl; 
        Ny=stoi(val[3]); 
       cout << Ny << endl; 
        Px=stoi(val[4]); 
       cout << Px << endl; 
        Py=stoi(val[5]);
      cout << Py << endl;  
        dt=stod(val[6]);
       cout << dt << endl; 
        T=stoi(val[7]); 
       cout << T << endl; 
        Re=stod(val[8]); 
       cout << Re << endl; 

      v=new double[Nx*Ny]; //Initialize vorticity values 
      s=new double[Nx*Ny]; 
      v1=new double[(Nx-2)*(Ny-2)]; //Initialize inner vorticity and streamfunction values 
      s1=new double[(Nx-2)*(Ny-2)]; 
      
      double* A=new double[(Nx-2)*(Ny-2)*(Nx-2)*(Ny-2)]; 
      
//      A= PoissonSolver::BMatrixPoisson(4,4,0.05,0.05); 
//      printmat(4-2,4-2,A); 

}
      


void LidDrivenCavity::Integrate()
{    
      double dx=Lx/(Nx-1.0); 
       double dy=Lx/(Ny-1.0);
       double U=1.0; 
       
       if (dt>=Re*dx*dy/4){
           cout << "Chosen value of dt too large" << endl; 
       }
    //Boundary conditions 
    for (int i=0;i<Nx;i++){
        v[(i+1)*(Ny-1)+i]=(2/pow(dy,2))*(s[(i+1)*(Ny-1)]-s[(i+1)*(Ny-1)-1])-(2*U/dy); //Top
        v[i*Ny]=(2/pow(dy,2))*(s[i*Ny+1]-s[i*Ny+2]); //Bottom 
    }
    for (int i=0;i<Ny;i++){
        v[i]=(2/(pow(dx,2)))*(s[i]-s[i+Ny]); //left 
        v[Ny*(Nx-1)+i]=(2/pow(dx,2))*(s[(Nx-1)*Ny+i]-s[(Nx-2)*Ny+i]); //Right 
    }

     cout << endl << endl; 
    
    //Inner vorticity values at current timestep 
    for (int i=1;i<Nx-1;i++){
        for (int j=1;j<Ny-1;j++){
            v[i*Ny+j]=-(s[(i*Ny)+j-1]-2*s[i*Ny+j]-s[(i*Ny)+j+1])/(pow(dy,2))-(s[(i-1)*Ny+j]-2*s[i*Ny+j]-s[(i+1)*Ny+j])/(pow(dx,2)); 
        }
    }
    
    //Inner vorticity values at next timestep 
    for (int i=1;i<Nx-1;i++){
        for (int j=1;j<Ny-1;j++){ //Create temporary variables to ensure ease of reading 
            double temp1=(s[(i)*Ny+j+1]-s[i*Ny+j-1])*(v[(i+1)*Ny+j]-v[(i-1)*Ny+j])/(4*dy*dx); 
            double temp2=(s[(i+1)*Ny+j]-s[(i-1)*Ny+j])*(v[i*Ny+j+1]-v[i*Ny+j-1])/(4*dy*dx); 
            double temp3=(1/Re)*((v[(i*Ny)+j-1]-2*v[i*Ny+j]-v[(i*Ny)+j+1])/(pow(dy,2))+(v[(i-1)*Ny+j]-2*v[i*Ny+j]-v[(i+1)*Ny+j])/(pow(dx,2))); 
            v[i*Ny+j]=v[i*Ny+j]+dt*(temp2-temp1+temp3); 
            
        }
    }
    
    //printmat(Nx*Ny,1,v); 
    
    
    
}

#include "LidDrivenCavity.h"
#include <iostream>

LidDrivenCavity::LidDrivenCavity() //Constructor 
{
    
}

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

void LidDrivenCavity::Initialise(string *val) //Initialise all parameters
{    
    
       double Lx=stod(val[0]); 
       cout << Lx << endl; 
       double Ly=stod(val[1]); 
       cout << Ly << endl; 
       int Nx=stoi(val[2]); 
       cout << Nx << endl; 
       int Ny=stoi(val[3]); 
       cout << Ny << endl; 
       int Px=stoi(val[4]); 
       cout << Px << endl; 
       int Py=stoi(val[5]);
      cout << Py << endl;  
       double dt=stod(val[6]);
       cout << dt << endl; 
       double T=stoi(val[7]); 
       cout << T << endl; 
       double Re=stod(val[8]); 
       cout << Re << endl; 

}

void LidDrivenCavity::Integrate()
{
}
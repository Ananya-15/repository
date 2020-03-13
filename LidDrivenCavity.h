//#pragma once
//#include "PoissonSolver.h"
#include <string>
using namespace std;



class LidDrivenCavity
{
public:
    LidDrivenCavity();
    ~LidDrivenCavity(); 
    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);
    
    //initialise and integrate perform solver calculations to determine streamfunction and vorticity values through time 
    void Initialise(string *val);
    void Integrate();
    
    //Functions which should belong to Poisson Solver class 
     double* PoissonMatrix(int nx,int ny, int var, double dx, double dy); 
     double* SolveMatrix(double* A, double* w1, double* s1, int var); 
     //double* SolveMatrix(double* A, double* w1, double* s1,int var,int Nx,int Ny); 
     
    //Calculates boundary conditions for the lid cavity problem 
    void BoundaryConditions(int Nx, int Ny, double* v, double dx, double dy, double dt, double U); 

    //Calculates inner vorticty at current timestep 
    void InnerVorticity(double* v, double* s, int Nx, int Ny, double dx, double dy); 
    
    //Calculates inner voritcity at next timestep 
    void NextInnerVorticity(double* v, double* s, int Nx, int Ny, double dx, double dy, double dt, double Re); 
    
    //Obtains inner voriticyt values for calculation of streamfunction at th enext time step 
    void RecoverInnerVorticity(double*v, double*v1, int Nx, int Ny); 
    
    //Updates streamfunction values for the next time step 
    void UpdateInnerStream(double*s, double*s1, int Nx, int Ny); 
    

private:
    double* v = nullptr; //vorticity
    double* s = nullptr; //streamfunction 
    double *v1= nullptr; //Define inner vorticity and streamfunctions
    double *s1= nullptr; 
    double* A=nullptr; 
  //  PoissonSolver* Poison=nullptr; 
    
    
    
    double dt;
    double dx; 
    double dy; 
    double U; 
    double T;
    double Lx;
    double Ly;
    double Re;
    
    int    Nx;
    int    Ny;
    int    nx; 
    int    ny; 
    int    var; 
    double Px; 
    double Py; 
};

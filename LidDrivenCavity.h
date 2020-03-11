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

    void Initialise(string *val);
    void Integrate();

    // Add any other public functions

private:
    double* v = nullptr; //vorticity
    double* s = nullptr; //streamfunction 
    double *v1= nullptr; //Define inner vorticity and streamfunctions
    double *s1= nullptr; 

    double dt;
    double T;
    double Lx;
    double Ly;
    double Re;

    
protected: //Allow protected access so that it can be accessed by Poisson solver class 
    
    int    Nx;
    int    Ny;
    double Px; 
    double Py; 
};

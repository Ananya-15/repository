//#pragma once
//#include "PoissonSolver.h"
#include <string>
#include <iostream>
//#include <mpi.h>

using namespace std;



class LidDrivenCavityExp
{
public:
    LidDrivenCavityExp();
    ~LidDrivenCavityExp(); 
    void SetDomainSize();
    void SetGridSize();
//    void SetPartitions(int px, int py); 
//    void SetGlobalDomainSize();
    
    //initialise and integrate perform solver calculations to determine streamfunction and vorticity values through time 
    void SubDomainInfo(); 
    void AssignGlobal(string* val);
    void Initialise();
    void Integrate();
    void BoundaryVectors();
    
    int CheckParallel(); 
    void CommunicateBound(); //Function to communicate local edge values of stramfunction and voritcity for each subdomain 
    
    
    //Functions which should belong to Poisson Solver class 
    // double* PoissonMatrix(double dx, double dy); 
     double* SolveMatrix(double* A, double* w1, double* s1, int var); 
//     double* SolveMatrix(double A, double* w1, double* s1,int var,int Nx,int Ny); 
     
    //Calculates boundary conditions for the lid cavity problem 
    void BoundaryConditions(); 
    //int rank, int Nx, int Ny, double* v, double dx, double dy, double dt, double U

    //Calculates inner vorticty at current timestep 
    
    int* GetIndex();
    
    void InnerVorticity(int* ind, int size1, int size2); 
    
    //Calculates inner voritcity at next timestep 
    void NextInnerVorticity(int* ind, int size1,int size2); 
    
    void CopyVorticity(double*arr1,double* arr2, int var); 
    
    //Obtains inner voriticyt values for calculation of streamfunction at th enext time step 
    void RecoverInnerVorticity(double*v, double*v1, int Nx, int Ny); 
    
    //Updates streamfunction values for the next time step 
    void UpdateInnerStream(double*s, double*s1, int Nx, int Ny); 
    
    
    

private:
    //Define pointers used throughout the process 
    double* v = nullptr; //vorticity
    double* vnew=nullptr; //Vorticity at new timestep 
    double* s = nullptr; //streamfunction 
    double *v1= nullptr; //Define inner vorticity and streamfunctions
    double *s1= nullptr; 
    double* A=nullptr; 
    
    double* vtemp=nullptr; 
    double* vtemp1=nullptr; 
    double* vtemp2=nullptr; 
    double* vtemp3=nullptr; 
    double* stemp=nullptr; 
    double* stemp1=nullptr; 
    double* stemp2=nullptr;
    double* stemp3=nullptr;  

    // Global variable definition //Some values may be redundant 
    
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
//    int    nx; 
//    int    ny; 
//    int    var; 
    int Px; //Partition number chosen by user 
    int Py; 
    
    int Nxloc;  //No. of columns and rows determined using fortran numroc function  //definitions for parallelisation 
    int Nyloc; 
    double Lxloc; 
    double Lyloc; 
    
    double* v_left=nullptr;
    double* s_left=nullptr;
    double* y_left=nullptr; //Dummy for storing values during communication //Using y to store streamfunction values (Innervorticity) and x to store vorticity values
    double* x_left=nullptr;
    
    double* v_right=nullptr;
    double* s_right=nullptr;
    double* y_right=nullptr; //Dummy for storing values during communication 
    double* x_right=nullptr;
    
    
    double* v_top=nullptr;
    double* s_top=nullptr;
    double* y_top=nullptr; //Dummy for storing values during communication 
    double* x_top=nullptr;
    
    
    double* v_bot=nullptr;
    double* s_bot=nullptr;
    double* y_bot=nullptr; //Dummy for storing values during communication 
    double* x_bot=nullptr;
//    int nxloc; 
//    int nyloc; 
    //Subdomain array size and length definition 
    
     int rank, nprocs, retval_rank, retval_size;
     int ctxt;
     int myrow, mycol, row_loc, col_loc; 
     
    
};

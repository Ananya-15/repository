//#pragma once

#include <string>
#include <iostream>

//Created by Ananya Dubey. LidDrivenCavtyExp.h contains class definition of LidDrivenCavityExp. Member functions and global variables are declared here, but 
//full definition has been provided in LidDrivenCavityExp.cpp 

using namespace std;



class LidDrivenCavityExp
{
public:
    LidDrivenCavityExp(); //Defualt constructors and destructors 
    ~LidDrivenCavityExp(); 
    
     void SubDomainInfo(); //Use cblacs functions to determine local grid size, grid length and position of each subdomain on the global domain 
    void AssignGlobal(string* val); //Assign values required and used by all processes 
    void SetDomainSize(); //Function to assign Lxloc and Lyloc 
    void SetGridSize(); //Function to assign no. of rows and columns in each subdomain 
    int CheckParallel(); //Function to check any errors thrown by MPI. If error occurs, program terminates with an error message 
    
    
   
    void Initialise(); //Defines size of pointers and initializes other variables to be used by other member functions 
    void Integrate(); //Used to determine how streamfunction and vorticity change over time, and obtain final values 
    //void BoundaryVectors();
    void BoundaryVectorsGen(double* arr, double* arr_left, double* arr_right, double* arr_bot, double* arr_top); //Function assigns 
    
   
    void CommunicateStreamBound(); //Function to communicate local edge values of stramfunction and voritcity for each subdomain 
    void CommunicateVortBound();
    void Communicate(double* arr_top,double* buf_top, double*arr_bot, double* buf_bot, double* arr_left, double* buf_left, double* arr_right, double* buf_right); 
    void OutputValues(); 
    
    
    //Functions which should belong to Poisson Solver class 
    // double* PoissonMatrix(double dx, double dy); 
     double* SolveMatrix(double* A, double* w1, double* s1, int var); 
//     double* SolveMatrix(double A, double* w1, double* s1,int var,int Nx,int Ny); 
     
    //Calculates boundary conditions for the lid cavity problem 
    void BoundaryConditions(); 
    //int rank, int Nx, int Ny, double* v, double dx, double dy, double dt, double U

    //Calculates inner vorticty at current timestep 
    
    int* GetIndex();
    
    //void UpdateTemps(int* ind);
    
    void InnerVorticity(int* ind, int size1, int size2); 
    
    //Calculates inner voritcity at next timestep 
    void NextInnerVorticity(int* ind, int size1,int size2); 
    
    //void CopyVorticity(double*arr1,double* arr2, int var); 
    
    //Obtains inner voriticyt values for calculation of streamfunction at th enext time step 
   // void RecoverInnerVorticity(double*v, double*v1, int Nx, int Ny); 
    
    //Updates streamfunction values for the next time step 
    void UpdateInnerStream(double*s, double*s1, int Nx, int Ny); 
    
    void UpdateVorticity(int* ind,int size1,int size2); 
    
    double* Jacobi(double* s,double*v,int* ind,int size1, int size2, int Nyloc, int Nxloc); 
    
    
    

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

    //Subdomain array size and length definition 
    
     int rank, nprocs, retval_rank, retval_size;
     int ctxt;
     int myrow, mycol, row_loc, col_loc; 
     
    
};

#include <iostream>
#include <cblas.h>
#include <mpi.h> 
//#include "PoissonSolver.h"
#include "LidDrivenCavityExp.h"
#include <cassert>

using namespace std; 

extern "C" {
  void Cblacs_pinfo(int*, int*); //Definition of functions required to split up domains based on user inputs 
  void Cblacs_get(int, int, int*);
  void Cblacs_gridinit(int*, char const*, int, int);
  void Cblacs_gridinfo(int, int*, int*, int*, int*);
  void Cblacs_barrier(int , char*);
  void Cblacs_gridexit(int);
  void Cblacs_exit(int);
  
  int numroc_(int const& n, int const& nb, int const& iproc, int const& isproc, int const& nprocs);
}



int main(int argc, char **argv)
{   
        MPI_Init(&argc, &argv);
        
        int n=argc/2; 
        string *val=new string [n];
        
             for (int i=0;i<n;i++){   //Take user input, create array val with only the numeric values 
                    val[i]=argv[2*(i+1)]; 
             }        
            
         LidDrivenCavityExp* solver = new LidDrivenCavityExp(); //Instantiate 
         
         
         
         solver->AssignGlobal(val); //Assign variables such as Re, Nx, Ny, etc. which will be common to all subdomains 
         
         
         solver->SubDomainInfo(); //returns process grid coordinates 
    //     
         solver->SetDomainSize(); //Determine Lx and Ly of subdomain (Lxloc and Lyloc)
    //     
         solver->SetGridSize(); //Set local grid size for each subdomain (Defined as Nxloc and Nyloc)
    //     
         solver->Initialise(); //Define array size of vorticity and streamfunction variables 
    //     
         solver->Integrate(); //Determines streamfunction and vorticity values at the final timestep 
         
         solver->OutputValues(); //Gather values from all processes, print out to data file to be used for matlab plots 
              
         
         MPI_Finalize(); //Close parallelisation process 
}
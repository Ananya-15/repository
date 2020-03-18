#include <iostream>
//#include <math.h>
//#include "PrintMat.h"
#include <cblas.h>
#include <mpi.h> 
//#include "PoissonSolver.h"
#include "LidDrivenCavityExp.h"
#include <cassert>

using namespace std; 

extern "C" {
  void Cblacs_pinfo(int*, int*);
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
     for (int i=0;i<n;i++){
                    val[i]=argv[2*(i+1)]; 
     }        
    
     int rank, nprocs, retval_rank, retval_size; //Check if parallelisation successful 
     
     retval_rank = MPI_Comm_rank(MPI_COMM_WORLD, &rank); // zero-based //Initial rank is defined as zero 
     retval_size = MPI_Comm_size(MPI_COMM_WORLD, &nprocs); //Obtain rank and process number 
     if (retval_rank == MPI_ERR_COMM || retval_size == MPI_ERR_COMM) {
     std::cout << "Invalid communicator" << std::endl;
     return 1;
     }    
    
     LidDrivenCavityExp* solver = new LidDrivenCavityExp(); //Instantiate 
     
     solver->AssignGlobal(val); //Assign variables such as Re, Nx, Ny, etc. 
     
     solver->SubDomainInfo(rank,nprocs); //returns process grid coordinates 
     
     solver->SetDomainSize(); //Determine Lx and Ly of subdomain 
     
     solver->SetGridSize(); //Set local grid size Nx and Ny 
     
     solver->Initialise(rank); //Define array size of voriticity and streamfunction variables 
     
     cout << "After boundary condition imposition" << endl;
     
     solver->Integrate(rank); 
     
     

     
     
     
    
//     if (myrow>=0){ //myrow and mycol are local coordinates for processes 
//         cout << "Rank " << rank << " has coordinates " << myrow << ", " << mycol << endl << endl; 
//     }
     
      //Instantiate new class within each process 
      
     
     
     
     
     
    
     
     
     
//     
//     solver-> 
          
     
     MPI_Finalize(); 
}
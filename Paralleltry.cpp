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
    
    int rank, nprocs, ctxt, retval_rank, retval_size; //nprocs=no. of processors 
    
    double Lx, Ly, dt, Re; //Command line arguments 
    
    int Nx, Ny, Px, Py, T; 
    
    int myrow, mycol, row_loc, col_loc; //Row_loc and col_loc = Px and Py //myrow and mycol show process coordinates, where no. of processes are 
    //considered to be in a 2D array format 
    
    retval_rank = MPI_Comm_rank(MPI_COMM_WORLD, &rank); // zero-based //Initial rank is defined as zero 
     retval_size = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
     if (retval_rank == MPI_ERR_COMM || retval_size == MPI_ERR_COMM) {
     std::cout << "Invalid communicator" << std::endl;
     return 1;
     }

    int n=argc/2; 
    string *val=new string [n];
     for (int i=0;i<n;i++){
                    val[i]=argv[2*(i+1)]; 
            }
              Lx=stod(val[0]); 
              Ly=stod(val[1]); 
              Nx=stoi(val[2]); 
              Ny=stoi(val[3]); 
              Px=stoi(val[4]); 
              Py=stoi(val[5]);
              dt=stod(val[6]);
              T=stod(val[7]); 
              Re=stod(val[8]); 
    LidDrivenCavityExp* solver = new LidDrivenCavityExp(); //Instantiate 
    
   //  solver->AssignGlobal(argc,argv); //Assigns global value necessary for calculation 
     
     row_loc=Px; //No. local rows and columns expected due to partitioning 
     col_loc=Py; 
     
     double xlen=Lx/Px; 
     double ylen=Ly/Py; 
     
//     cout << row_loc << endl; 
//     cout << col_loc << endl; 
     
     Cblacs_pinfo(&rank,&nprocs); 
     assert(nprocs>=row_loc*col_loc); //Assertion statement to validate user input np for chosen array size nx and ny 
     
     Cblacs_get(0,0,&ctxt); //Attains system context 
     Cblacs_gridinit(&ctxt,"Col-major",row_loc,col_loc); //Creates process grid for a given context 
     
     Cblacs_gridinfo(ctxt,&row_loc,&col_loc,&myrow,&mycol); //Gridinfo 
     
     int nr=numroc_(Ny,1,myrow,0,row_loc);
     int nc=numroc_(Nx,1,mycol,0,col_loc); 
     
     cout << "Local row length: " << nr << endl; 
     cout << "Local Column length: " << nc << endl;  
     
//     if (myrow>=0){ //myrow and mycol are local coordinates for processes 
//         cout << "Rank " << rank << " has coordinates " << myrow << ", " << mycol << endl << endl; 
//     }
     
      //Instantiate new class within each process 
      
     solver->AssignGlobal(val);
     
     solver->SetDomainSize(xlen,ylen); 
     
     solver->SetGridSize(nr,nc);
    
     solver->Initialise(rank); 
     cout << "After boundary condition imposition" << endl;
     solver->Integrate(rank); 
     
     
//     
//     solver-> 
          
     
     MPI_Finalize(); 
}
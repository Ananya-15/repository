#include "LidDrivenCavityExp.h"
//#include "PoissonSolver.h"
#include <iostream>
#include <math.h>
#include "PrintMat.h"
#include <cblas.h>
#include <mpi.h> 
#include <cassert>


using namespace std; 
 #define F77NAME(x) x##_ //Definition for using lapack solver 
 extern "C" {
 // LAPACK routine for solving systems of linear equations
 void F77NAME(dsysv)(const char& uplo, const int& n, const int& nrhs,const double * A, const int& lda,int *ipiv, double* B, 
 const int& ldb, double * work, const int& lwork, int& info);
 }
 
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
 
 //Create default destructor and constructor 
LidDrivenCavityExp::LidDrivenCavityExp()=default; //Constructor 


LidDrivenCavityExp::~LidDrivenCavityExp()=default;  //Destructor 

void LidDrivenCavityExp::AssignGlobal(string *val){

              Lx=stod(val[0]); 
              Ly=stod(val[1]); 
              Nx=stoi(val[2]); 
              Ny=stoi(val[3]); 
              Px=stoi(val[4]); 
              Py=stoi(val[5]);
              dt=stod(val[6]);
              T=stod(val[7]); 
              Re=stod(val[8]);
}

void LidDrivenCavityExp::SubDomainInfo(int rank,int nprocs){
     
     row_loc=Px; //No. local rows and columns expected due to partitioning in process 
     col_loc=Py; 
     
     Cblacs_pinfo(&rank,&nprocs); 
     assert(nprocs>=row_loc*col_loc); //Assertion statement to validate user input np for chosen array size nx and ny 
     
     Cblacs_get(0,0,&ctxt); //Attains system context 
     Cblacs_gridinit(&ctxt,"Col-major",row_loc,col_loc); //Creates process grid for a given context 
     
     Cblacs_gridinfo(ctxt,&row_loc,&col_loc,&myrow,&mycol);
     
    
}


void LidDrivenCavityExp::SetDomainSize()
{  
   Lxloc=Lx/Px; 
   Lyloc=Ly/Py; 
}



void LidDrivenCavityExp::SetGridSize() //Solution created assuming Ny rows and Nx columns 
{
      Nyloc=numroc_(Ny,1,myrow,0,row_loc);
      Nxloc=numroc_(Nx,1,mycol,0,col_loc); 
      
     cout << "Local row length: " << Nyloc << endl; //Verify parallelisation split by checking local subdomain size 
     cout << "Local Column length: " << Nxloc << endl;  
     
}


void LidDrivenCavityExp::Initialise(int rank) //Use initialise to initialise all pointer values 
{    
      nxloc=Nxloc-2; 
      nyloc=Nyloc-2; 
      var=nxloc*nyloc; //
      v=new double[Nxloc*Nyloc]; //Initialize vorticity values 
      vnew=new double[Nxloc*Nyloc];//Voriticity at new timestep value 
      s=new double[Nxloc*Nyloc]; //Initialize streamfunction value 
      v1=new double[var]; //Initialize inner vorticity and streamfunction values 
      s1=new double[var]; 
      
      for (int i=0;i<Nxloc;i++){
          for (int j=0;j<Nyloc;j++){
              v[i*Nyloc+j]=1; 
          }
      }
      
      for (int i=0;i<Nxloc*Nyloc;i++){
          s[i]=0; 
      }
      
//      cout << "Initialised voriticity for rank " << rank << " is: " << endl << endl; 
//      printmat(Nyloc,Nxloc,v); 
//      cout << endl << endl; 
      
        dx=Lx/(Nx-1.0); //Determine other global variables used across several class member functions 
        dy=Lx/(Ny-1.0); 
        U=1.0;
        
        if (dt>=Re*dx*dy/4){
           throw logic_error("Chosen value of dt too large"); 
//     
       }
      
}

void LidDrivenCavityExp::BoundaryVectors(double* arr, int Nxloc; int Nyloc){
    col_left=
    col_right=
    col_top=
    col_bottom= 
    
}

 void LidDrivenCavityExp::BoundaryConditions(int rank,int Nxloc, int Nyloc, double* v, double dx, double dy, double dt, double U){//){
      
      if (rank % Py==0) {
         if (Nyloc==1){ //Communication necessary
            

          
          
            
        }
         else {
              for (int i=0;i<Nxloc;i++){
                         v[i*Nyloc]=(2.0/pow(dy,2))*(s[i*Nyloc]-s[i*Nyloc+1]); //Bottom 
                         cout << "Value is: " << Ly << endl; 
                         
                    }
                    
              }
         }

       if (rank % (Py-1)==0) {
         if (Nyloc==1){ //Communication necessary
          
        }
         else {
              for (int i=0;i<Nxloc;i++){
                        v[(i+1)*(Nyloc-1)+i]=(2/pow(dy,2))*(s[(i+1)*(Nyloc-1)]-s[(i+1)*(Nyloc-1)-1])-(2*U/dy); //Top
                    }
                    
              }
         }
//         
//         
         if  (rank<Py){
             if (Nxloc==1){ //Communication necessary, call communicator class function
          
                }
              else {
               for (int i=0;i<Nyloc;i++){
                v[i]=(2/(pow(dx,2)))*(s[i]-s[i+Nyloc]); //left 
                }
              }
         }
//         
         if  (rank<(Px*Py) && rank>=(Px-1)*Py){
             if (Nxloc==1){ //Communication necessary, call communicator class functionccccccc
          
                }
              else {
               for (int i=0;i<Nyloc;i++){
                v[Nyloc*(Nxloc-1)+i]=(2/pow(dx,2))*(s[(Nxloc-1)*Nyloc+i]-s[(Nxloc-2)*Nyloc+i]); //Right 
                }
              }
         }
    if (rank==3){
    cout << "Rank is: " << rank << endl; 
    printmat(Nyloc,Nxloc,v);    
//    cout << "Streamfunctio: " <<endl; 
//    printmat(Nyloc,Nxloc,s);  
    }
 }
    
           


//void LidDrivenCavity::InnerVorticity(double* v, double* s, int Nx, int Ny, double dx, double dy){
//      for (int i=1;i<Nx-1;i++){
//                for (int j=1;j<Ny-1;j++){
//                    v[i*Ny+j]=-(s[(i*Ny)+j-1]-2*s[i*Ny+j]+s[(i*Ny)+j+1])/(pow(dy,2))-(s[(i-1)*Ny+j]-2*s[i*Ny+j]+s[(i+1)*Ny+j])/(pow(dx,2)); 
//                }
//      }
// }
 
// void LidDrivenCavity::NextInnerVorticity(double* v, double* s, int Nx, int Ny, double dx, double dy, double dt, double Re){
//        for (int i=1;i<Nx-1;i++){
//            for (int j=1;j<Ny-1;j++){ //Create temporary variables to ensure ease of reading 
//                double temp1=(s[(i)*Ny+j+1]-s[i*Ny+j-1])*(v[(i+1)*Ny+j]-v[(i-1)*Ny+j])/(4*dy*dx); 
//                double temp2=(s[(i+1)*Ny+j]-s[(i-1)*Ny+j])*(v[i*Ny+j+1]-v[i*Ny+j-1])/(4*dy*dx); 
//                double temp3=(1/Re)*((v[(i*Ny)+j-1]-2*v[i*Ny+j]+v[(i*Ny)+j+1])/(pow(dy,2))+(v[(i-1)*Ny+j]-2*v[i*Ny+j]+v[(i+1)*Ny+j])/(pow(dx,2))); 
//                vnew[i*Ny+j]=v[i*Ny+j]+dt*(temp2-temp1+temp3); 
//                
//            }
//        }
//
// }
// 


//void LidDrivenCavity::RecoverInnerVorticity(double*vnew, double*v1, int Nx, int Ny){
//     for(int i=1; i<Nx-1;i++){
//                for (int j=1;j<Ny-1;j++){
//                    v1[(i-1)*(Ny-2)+j-1]=vnew[(i)*Ny+j]; 
//                }
//            }
//            
//}


// void LidDrivenCavity::CopyVorticity(double*arr1,double* arr2, int var){ //For first one vnew=v, for second one, v=vnew 
//      cblas_dcopy(var,arr1,1,arr2,1); 
// }


//void LidDrivenCavity::UpdateInnerStream(double*s, double*s1, int Nx, int Ny){
//            for(int i=1; i<Nx-1;i++){
//                for (int j=1;j<Ny-1;j++){
//                    s[(i)*Ny+j]=s1[(i-1)*(Ny-2)+j-1]; 
//                }
//               }
//
//}



void LidDrivenCavityExp::Integrate(int rank)
{    
  //PoissonSolver* psolver=new PoissonSolver(); //Create new instance psolver of class Poisson Solver to implement solver calculations  

  //psolver->Initialise(Nx,Ny,dx,dy); //Initialises values within PoissonSolver class 
  double t=0; //First time step value  
  
   // while (t<T){
            cout << "Time step is: " << t << endl << endl; 
            
             LidDrivenCavityExp::BoundaryConditions(rank,Nxloc,Nyloc,v,dx,dy,dt,U); //Update with BCs //
//             
//             
//            
//            cout << "Updated vorticity boundar conditions: " << endl; 
//            cout << endl << endl; 
//            printmat(Ny,Nx,v); 
//            
//            
//             LidDrivenCavity::InnerVorticity(v,s,Nx,Ny,dx,dy); //Inner vorticity values at current timestep 
//             
//             LidDrivenCavity::CopyVorticity(v,vnew,Nx*Ny); 
//            cout << "Inner voritcity values at current timestep" << endl;    
//            cout << endl << endl; 
//            printmat(Ny,Nx,v);            
//  
//            LidDrivenCavity::NextInnerVorticity(v,s,Nx,Ny,dx,dy,dt,Re);  //Inner vorticity values at next timestep, updates vnew
//        
//            
//            cout << "Updated vorticity at new timestep: " << endl; 
//            printmat(Ny,Nx,vnew); 
//            cout << endl << endl; 
//            
//            LidDrivenCavity::RecoverInnerVorticity(vnew,v1,Nx,Ny); //Extract to obain inner voriticity values and save in matrix v1 
//            
//            LidDrivenCavity::CopyVorticity(vnew,v,Nx*Ny); //Update vorticity so that v=vnew 
//            
//            A=psolver->MatrixPoisson(); //Call function to obtain matrix A
//               
//            cout << "Works here 1" << endl;
//            
//            cout << "Updated vorticity at new timestep: " << endl; 
//            printmat(Ny,Nx,v); 
//            cout << endl << endl; 
//            
//             s1=psolver->matrixSolve(A,var,v,s1); //Time to implement lapack 
//            
//            //s1=LidDrivenCavity::SolveMatrix(A,v1,s1,var); 
//              cout << "Works here 2" << endl; 
//               
//            LidDrivenCavity::UpdateInnerStream(s,s1,Nx,Ny);  //Need to input s1 back in streamfunction 
//
//              cout << "Works here 3" << endl; 
//          t+=dt; 
//      }
//    cout << "Final streamfunction values: " << endl; 
//    printmat(Nx,Ny,v); //Check if final values make sense 
//
//  //delete [] A; 
//  delete [] v1; 
//  delete [] s1; 
//  delete [] v; 
//  delete [] s; 
//     
}
//      
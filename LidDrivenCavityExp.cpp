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

int LidDrivenCavityExp::CheckParallel(){
     //Check if parallelisation successful 
     
     retval_rank = MPI_Comm_rank(MPI_COMM_WORLD, &rank); // zero-based //Initial rank is defined as zero 
     retval_size = MPI_Comm_size(MPI_COMM_WORLD, &nprocs); //Obtain rank and process number 
     if (retval_rank == MPI_ERR_COMM || retval_size == MPI_ERR_COMM) {
     std::cout << "Invalid communicator" << std::endl;
     return 1;
     }  
     else {
         return 0;
     }
}

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

void LidDrivenCavityExp::SubDomainInfo(){
     
     row_loc=Py; //No. local rows and columns expected due to partitioning in process 
     col_loc=Px; 
     
     Cblacs_pinfo(&rank,&nprocs); 
     assert(nprocs>=row_loc*col_loc); //Assertion statement to validate user input np for chosen array size nx and ny 
     
     Cblacs_get(0,0,&ctxt); //Attains system context 
     Cblacs_gridinit(&ctxt,"Col-major",row_loc,col_loc); //Creates process grid for a given context 
     
     Cblacs_gridinfo(ctxt,&row_loc,&col_loc,&myrow,&mycol);
//     
//     cout << "Myrow is: " << myrow << endl; 
//     cout << "Mycol is: " << mycol << endl; 
//     cout << "Mstrix rank is: "<<rank << endl; 
//     cout << endl; 
//     
    
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
      
//     cout << "Local row length: " << Nyloc << endl; //Verify parallelisation split by checking local subdomain size 
//     cout << "Local Column length: " << Nxloc << endl;  
     
}


void LidDrivenCavityExp::Initialise() //Use initialise to initialise all pointer values 
{    
      
      v=new double[Nxloc*Nyloc]; //Initialize vorticity values 
      vnew=new double[Nxloc*Nyloc];//Voriticity at new timestep value 
      s=new double[Nxloc*Nyloc]; //Initialize streamfunction value 
//      v1=new double[var]; //Initialize inner vorticity and streamfunction values 
//      s1=new double[var]; 
//      
      v_left=new double[Nyloc]; 
      s_left=new double[Nyloc];
      y_left=new double[Nyloc]; //Column vector storing streamfunction values sent from left column of adjacent subdomain at boundary
      x_left=new double[Nyloc]; //Column vector storing vorticity values sent from left column of adjacent subdomain at boundary
      
      v_right=new double[Nyloc]; 
      s_right=new double[Nyloc];
      y_right=new double[Nyloc]; //Column vector storing streamfunction values sent from right column of adjacent subdomain at boundary
      x_right=new double[Nyloc];  //Column vector storing vorticity values sent from right column of adjacent subdomain at boundary
      
      v_top=new double[Nxloc]; 
      s_top=new double[Nxloc]; 
      y_top=new double[Nxloc]; //Row vector storing streamfunction values sent from top row of adjacent subdomain at boundary
      x_top=new double[Nxloc]; //Row vector storing vorticity values sent from top row of adjacent subdomain at boundary
      
      v_bot=new double[Nxloc]; 
      s_bot=new double[Nxloc]; 
      y_bot=new double[Nxloc]; //Row vector storing streamfunction values sent from top row of adjacent subdomain at boundary
      x_bot=new double[Nxloc]; //Row vector storing vorticity values sent from bottom row of adjacent subdomain at boundary
 
 
      for (int i=0;i<Nxloc*Nyloc;i++){
              v[i]=1; 
              s[i]=1; 
      }
//      for (int i=0;i<Nxloc;i++){
//          for (int j=0;j<Nyloc;j++){
//              s[i*Nxloc+j]=
//          }
//      }
      
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

void LidDrivenCavityExp::BoundaryVectors(double* v, double*s, int Nxloc, int Nyloc){ //Inputs local array arr and arr dimensions //Rewrites left, right, top and bottom vorticity and streamfunction boundary values 
//For each subdomain 
   for (int i=0;i<Nyloc;i++){
    v_left[i]=v[i];
    s_left[i]=s[i]; 
   }
    
    for (int i=0;i<Nyloc;i++){
    v_right[i]=v[Nyloc*(Nxloc-1)+i]; 
    s_right[i]=s[Nyloc*(Nxloc-1)+i];
    }
    
    for (int i=0;i<Nxloc;i++){
    v_bot[i]=v[i*Nyloc];
    s_bot[i]=s[i*Nyloc];
    }
    
    for (int i=0;i<Nxloc;i++){
    v_top[i]=v[(i+1)*(Nyloc-1)]; 
    s_top[i]=s[(i+1)*(Nyloc-1)]; 
    }

//    if (rank==1){  //Check value by printing 
//        cout << "Vright for this subdomain is: " << endl; 
//        printmat(Nyloc,1,v_right); 
//        }
}

void LidDrivenCavityExp::CommunicateBound(){ //As boundary values at each subdomain are required for calculating vorticity at current and next timesteps, function
//is created to communicate and store values between different subdomains
  if (Px>1){
    if (mycol<Px-1){

        MPI_Send(&s_left,Nyloc,MPI_DOUBLE,rank+Py,0,MPI_COMM_WORLD);
        MPI_Send(&v_left,Nyloc,MPI_DOUBLE,rank+Py,2,MPI_COMM_WORLD);
        cout << "Sends first" << endl; 
        MPI_Recv(&y_right,Nyloc,MPI_DOUBLE,rank+Py,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
        MPI_Recv(&x_right,Nyloc,MPI_DOUBLE,rank+Py,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
        cout << "Receives second" << endl; 
     }
   if (mycol>0){

    //       int tag=; 
    //       int tag1=1; 
           MPI_Recv(&y_left,Nyloc,MPI_DOUBLE,rank-Py,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
           MPI_Recv(&x_left,Nyloc,MPI_DOUBLE,rank-Py,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
                   cout << "Receives first" << endl; 
           MPI_Send(&s_right,Nyloc,MPI_DOUBLE,rank-Py,1,MPI_COMM_WORLD); 
           MPI_Send(&v_right,Nyloc,MPI_DOUBLE,rank-Py,3,MPI_COMM_WORLD); 
                   cout << "Sends second" << endl; 
       }
 }
//  
  if (Py>1){
   if (myrow<Py-1){
       MPI_Send(&s_bot,Nxloc,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD); 
       MPI_Send(&v_bot,Nxloc,MPI_DOUBLE,rank+1,2,MPI_COMM_WORLD); 
       cout << "Sends third" << endl; 
//        cout << "Rank+1 is: " << rank+1 << endl; 
       MPI_Recv(&y_top,Nxloc,MPI_DOUBLE,rank+1,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
       MPI_Recv(&x_top,Nxloc,MPI_DOUBLE,rank+1,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
             cout << "Receives fourth" << endl; 
    
      
   }
   if (myrow>0){
       MPI_Recv(&y_bot,Nxloc,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
       MPI_Recv(&x_bot,Nxloc,MPI_DOUBLE,rank-1,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
       cout << "Receives third" << endl; 
       MPI_Send(&s_top,Nxloc,MPI_DOUBLE,rank-1,1,MPI_COMM_WORLD); 
       MPI_Send(&v_top,Nxloc,MPI_DOUBLE,rank-1,3,MPI_COMM_WORLD); 
       cout << "Sends Fourth" << endl;
//       cout << rank-1<< endl; 
     
      }

  }
  
    
}

 void LidDrivenCavityExp::BoundaryConditions(int rank, int Nxloc, int Nyloc, double* v, double dx, double dy, double dt, double U){//){
      
      if (rank % Py==0) { //Bottom BC points 
         if (Nyloc==1){ //For the case when communication with another subdomain is necessary 
//                int src=rank+1; 
//                int dest=rank;
//                int tag=0; 
//                cout << "Nxloc is: " << Nxloc << endl; 
//                double* y=new double[Nxloc]; //Create variable y to which the message can be copied 
//                y=LidDrivenCavityExp::CommunicateBound(s_bot,y,Nxloc,src,dest,tag); 
//                cout << "Retrieved message is: " << endl; 
//                printmat(1,Nxloc,y); //For some reason, y isn't received correctly 
//                    for (int i=0;i<Nxloc;i++){
//                                v[i]=(v[i]-y[i])*(2.0/pow(dy,2)); 
//                            }
//              printmat(1,Nxloc,y); 
//              printmat(1,Nxloc,v); 
//              
//             delete[] y; 
         
           }
            
         else {
              for (int i=0;i<Nxloc;i++){
                         v[i*Nyloc]=(2.0/pow(dy,2))*(s[i*Nyloc]-s[i*Nyloc+1]); //Bottom 
                }
                    
              }
     }
//
       if (rank % (Py-1)==0) { //Top BC points 
             if (Nyloc==1){ //Communication necessary
//                int src=rank-1; 
//                int dest=rank;
//                int tag=0; 
//                double* y=new double[Nxloc];
//                y=LidDrivenCavityExp::CommunicateBound(s_top,y,Nxloc,src,dest,tag); 
//                        for (int i=0;i<Nxloc;i++){
//                            v[i]=(v[i]-y[i])*(2.0/pow(dy,2)); 
//                        }
//              delete[] y; 
              
            }
             else {
                  for (int i=0;i<Nxloc;i++){
                            v[(i+1)*(Nyloc-1)+i]=(2/pow(dy,2))*(s[(i+1)*(Nyloc-1)]-s[(i+1)*(Nyloc-1)-1])-(2*U/dy); //Top
                        }
                        
                  }
         }
//////         
//////         
        if  (rank<Py){ //left 
             if (Nxloc==1){ //Communication necessary, call communicator class function
//                int src=rank+Py; 
//                int dest=rank;
//                int tag=0; 
//                double* y=new double[Nyloc];
//                y=LidDrivenCavityExp::CommunicateBound(s_left,y,Nyloc,src,dest,tag); 
//                     for (int i=0;i<Nyloc;i++){
//                        v[i]=(v[i]-y[i])*(2.0/pow(dx,2)); 
//                    }
//          
//                 delete[] y; 
                }
                
              else {
                   for (int i=0;i<Nyloc;i++){
                    v[i]=(2/(pow(dx,2)))*(s[i]-s[i+Nyloc]); 
                    }
              }
         }
////         
//     
      if  (rank<(Px*Py) && rank>=(Px-1)*Py){
             if (Nxloc==1){ //Communication necessary, call communicator class functionccccccc
//                int src=rank-Py; 
//                int dest=rank;
//                int tag=0; 
//                double* y=new double[Nyloc];
//                y=LidDrivenCavityExp::CommunicateBound(s_right,y,Nyloc,src,dest,tag); 
//                     for (int i=0;i<Nyloc;i++){
//                        v[i]=(v[i]-y[i])*(2.0/pow(dx,2)); 
//                    }
//                   delete[] y;
                }

              else {
               for (int i=0;i<Nyloc;i++){
                v[Nyloc*(Nxloc-1)+i]=(2/pow(dx,2))*(s[(Nxloc-1)*Nyloc+i]-s[(Nxloc-2)*Nyloc+i]); //Right 
                }
              }
         }



//    if (rank==2){
//    cout << "Rank is: " << rank << endl; 
//    printmat(Nyloc,Nxloc,v);    
//    cout << "Streamfunctio: " <<endl; 
//    printmat(Nyloc,Nxloc,s);  
    
 
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



void LidDrivenCavityExp::Integrate()
{    
  //PoissonSolver* psolver=new PoissonSolver(); //Create new instance psolver of class Poisson Solver to implement solver calculations  

  //psolver->Initialise(Nx,Ny,dx,dy); //Initialises values within PoissonSolver class 
  double t=0; //First time step value  
  
   // while (t<T){
//            cout << "Time step is: " << t << endl << endl; 
             LidDrivenCavityExp::BoundaryVectors(v,s,Nxloc,Nyloc); //Obtain rows and columns located at the boundaries of subdomains 
            
             LidDrivenCavityExp::BoundaryConditions(rank,Nxloc,Nyloc,v,dx,dy,dt,U); //Update with BCs //
             
             LidDrivenCavityExp::CommunicateBound(); 
             
             
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
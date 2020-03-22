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
              v[i]=0; 
              s[i]=0; 
              vnew[i]=0; 
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

void LidDrivenCavityExp::BoundaryVectors(){ //Inputs local array arr and arr dimensions //Rewrites left, right, top and bottom vorticity and streamfunction boundary values 
//For each subdomain 
   for (int i=0;i<Nyloc;i++){
    v_left[i]=v[i];
    s_left[i]=s[i]; 
   }
    
    for (int i=0;i<Nyloc;i++){
    v_right[i]=v[Nyloc*(Nxloc-1)]; 
    s_right[i]=s[Nyloc*(Nxloc-1)];
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
//y_left stores s_left from adjacent subdomain 
//x_left stores v_left from adjacent subdomain 
//y_top stores s_bot from subdomain right below //bottom values are stored on top row 
//y_bot stores s_top from subomain right above  //top values are stored on bottom row 

  if (Px>1){
    if (mycol<Px-1){

        MPI_Send(s_right,Nyloc,MPI_DOUBLE,rank+Py,0,MPI_COMM_WORLD);
        MPI_Send(v_right,Nyloc,MPI_DOUBLE,rank+Py,2,MPI_COMM_WORLD);
       // cout << "Sends first" << endl; 
        MPI_Recv(y_left,Nyloc,MPI_DOUBLE,rank+Py,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
        //printmat(Nyloc,1,y_right);
        MPI_Recv(x_left,Nyloc,MPI_DOUBLE,rank+Py,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
        //printmat(Nyloc,1,x_right); 
        //cout << "Receives second" << endl; 
     }
   if (mycol>0){
           MPI_Recv(y_right,Nyloc,MPI_DOUBLE,rank-Py,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
           MPI_Recv(x_right,Nyloc,MPI_DOUBLE,rank-Py,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
              //     cout << "Receives first" << endl; 
           MPI_Send(s_left,Nyloc,MPI_DOUBLE,rank-Py,1,MPI_COMM_WORLD); 
           //printmat(Nyloc,1,s_right); 
           MPI_Send(v_left,Nyloc,MPI_DOUBLE,rank-Py,3,MPI_COMM_WORLD); 
                //   cout << "Sends second" << endl; 
       }
 }
//  
  if (Py>1){
   if (myrow<Py-1){
       MPI_Send(s_top,Nxloc,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD); 
       MPI_Send(s_top,Nxloc,MPI_DOUBLE,rank+1,2,MPI_COMM_WORLD); 
       //cout << "Sends third" << endl; 
//        cout << "Rank+1 is: " << rank+1 << endl; 
       MPI_Recv(y_bot,Nxloc,MPI_DOUBLE,rank+1,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
      MPI_Recv(x_bot,Nxloc,MPI_DOUBLE,rank+1,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
            // cout << "Receives fourth" << endl; 
    
//      
   }
   if (myrow>0){
       MPI_Recv(y_top,Nxloc,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
      MPI_Recv(x_top,Nxloc,MPI_DOUBLE,rank-1,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
//          if (rank==0){
//          cout << "Nxloc is: " << Nxloc << endl;
//          cout <<"received y_top is: " << endl; 
//          printmat(1,Nxloc,y_top); }
      // cout << "Receives third" << endl; 
       MPI_Send(s_bot,Nxloc,MPI_DOUBLE,rank-1,1,MPI_COMM_WORLD); 
       MPI_Send(v_bot,Nxloc,MPI_DOUBLE,rank-1,3,MPI_COMM_WORLD); 
      // cout << "Sends Fourth" << endl;
//       cout << rank-1<< endl; 
     
      }

  }
  
  MPI_Barrier(MPI_COMM_WORLD); 
  
    
}

 void LidDrivenCavityExp::BoundaryConditions(){
      
      if (myrow==0) { //Bottom BC points 
         if (Nyloc==1){ //For the case when communication with another subdomain is necessary 
                    for (int i=0;i<Nxloc;i++){
                                v[i]=(s[i]-y_bot[i])*(2.0/pow(dy,2)); //y_top stores s_top of row beneath 
                    }
           }
            
         else {
              for (int i=0;i<Nxloc;i++){
                         v[i*Nyloc]=(2.0/pow(dy,2))*(s[i*Nyloc]-s[i*Nyloc+1]); //Bottom 
                }
               
                    
              }
     }

       if (myrow==Py-1) { //Top BC points 
             if (Nyloc==1){ //Communication necessary
                        for (int i=0;i<Nxloc;i++){
                            v[i]=(s[i]-y_top[i])*(2.0/pow(dy,2))-(2*U/dy); //y_bot stores s_top 
                        }
            }
             else {
                  for (int i=0;i<Nxloc;i++){
                            v[(i+1)*(Nyloc-1)+i]=(2/pow(dy,2))*(s[(i+1)*(Nyloc-1)]-s[(i+1)*(Nyloc-1)-1])-(2*U/dy); //Top
                        }
                        
                  }
         }
        
        if  (mycol==0){ //left 
             if (Nxloc==1){ //Communication necessary, call communicator class function
                     for (int i=0;i<Nyloc;i++){
                        v[i]=(s[i]-y_right[i])*(2.0/pow(dx,2)); //Stores s_left from adjacent array
                    }
                }
//                
              else {
                   for (int i=0;i<Nyloc;i++){
                    v[i]=(2/(pow(dx,2)))*(s[i]-s[i+Nyloc]); 
                    }
              }
         }
//////         
////     
      if  (mycol==Px-1){
             if (Nxloc==1){ //Communication necessary, call communicator class function
                     for (int i=0;i<Nyloc;i++){
                        v[i]=(s[i]-y_left[i])*(2.0/pow(dx,2)); 
                    }
//                    cout << "v is: " << endl; 
//                               printmat(Nxloc,1,v);
                }

              else {
               for (int i=0;i<Nyloc;i++){
                v[Nyloc*(Nxloc-1)+i]=(2/pow(dx,2))*(s[(Nxloc-1)*Nyloc+i]-s[(Nxloc-2)*Nyloc+i]); //Right 
                }
              }
         }
    
 
 }

int* LidDrivenCavityExp::GetIndex(){  //Returns an array of indices consisting of i,imax, etc. to be used for calculating inner vorticity values 
  int* ind = new int[4]; 
         if (myrow==0 && mycol==0){
             ind[0]=1; 
             ind[1]=Nxloc; 
             ind[2]=1; 
             ind[3]=Nyloc;
         }
         else if (myrow==Py-1&&mycol==0){
             ind[0]=1; 
             ind[1]=Nxloc; 
             ind[2]=0; 
             ind[3]=Nyloc-1;
         }
         else if (myrow==Py-1&&mycol==Px-1){
             ind[0]=0; 
             ind[1]=Nxloc-1; 
             ind[2]=0; 
             ind[3]=Nyloc-1;
         }
         else if (mycol==Px-1&&myrow==0){
             ind[0]=0; 
             ind[1]=Nxloc-1; 
             ind[2]=1; 
             ind[3]=Nyloc;
         }
         else if (mycol==0) {
            ind[0]=1; 
             ind[1]=Nxloc; 
             ind[2]=0; 
             ind[3]=Nyloc; 
         }
         else if (myrow==0){
             ind[0]=0; 
             ind[1]=Nxloc; 
             ind[2]=1; 
             ind[3]=Nyloc;
         }
         else if (myrow==Py-1){
             ind[0]=0; 
             ind[1]=Nxloc; 
             ind[2]=0; 
             ind[3]=Nyloc-1;
         }
         else if (mycol==Px-1){
             ind[0]=0; 
             ind[1]=Nxloc-1; 
             ind[2]=0; 
             ind[3]=Nyloc;
         }
         else {
             ind[0]=0; 
             ind[1]=Nxloc; 
             ind[2]=0; 
             ind[3]=Nyloc;
             
         }
         return ind; 

}

//void LidDrivenCavityExp::Updatetemps(int i,int imax, int j, int jmax, int size1, int size2){
//
//for (i;i<imax;i++){
//        for (j;j<jmax;j++){
//            if (i==0){
//                stemp[i*size2+j]=y_right[j]; 
//                vtemp[i*size2+j]=x_right[j]; 
//            }
//            else {
//                stemp[i*size2+j]=s[(i-1)*Nyloc+j]; 
//                vtemp[i*size2+j]=v[(i-1)*Nyloc+j];
//            }
//            if (i==Nxloc-1){
//                stemp1[i*size2+j]=y_left[j]; 
//                vtemp1[i*size2+j]=x_left[j];
//            }
//            else {
//                stemp1[i*size2+j]=s[(i+1)*Nyloc+j];
//                vtemp1[i*size2+j]=v[(i+1)*Nyloc+j];
//            }
//            if (j==0) {
//                stemp2[i*size2+j]=y_top[i];
//                vtemp2[i*size2+j]=x_top[i];
//            }
//            else {
//                stemp2[i*size2+j]=s[i*Nyloc+j-1];
//                vtemp2[i*size2+j]=v[i*Nyloc+j-1];
//            }
//            if (j==Nyloc-1){
//                stemp3[i*size2+j]=y_bot[i];
//                vtemp3[i*size2+j]=x_bot[i]; 
//            }
//            else {
//                stemp3[i*size2+j]=s[i*Nyloc+j+1]; //After checking all the conditions, v is finally computed 
//                vtemp3[i*size2+j]=v[i*Nyloc+j+1];
//            }
//            temp1=(stemp3[i*size2+j]-stemp2[i*size2+j])*(vtemp1[i*size2+j]-vtemp[i*size2+j])/(4*dy*dx);
//            temp2=(stemp1[i*size2+j]-stemp[i*size2+j])*(vtemp3[i*size2+j]-vtemp2[i*size2+j])/(4*dy*dx); 
//            temp3=(1/Re)*((vtemp2[i*size2+j]-2*v[i*Nyloc+j]+vtemp3)/(pow(dy,2))+(vtemp-2*v[i*Nyloc+j]+vtemp1)/(pow(dx,2))); 
//            vnew[i*Nyloc+j]=v[i*Nyloc+j]+dt*(temp2-temp1+temp3); 
////            if (rank==4){
////                cout << "temp3 is: " << vtemp3 <<  endl; 
////            }
//            
//        }
//    }
//    
//}

//void LidDrivenCavityExp::InnerVorticity(int i, int imax, int j, int jmax,int size1, int size2){
//    for (i;i<imax;i++){
//        for (j;j<jmax;j++){
//            
//             v[i*size2+j]=-((stemp1[i*size2+j]-2*v[i*size2+j]+stemp[i*size2+j])/(pow(dx,2))-(stemp3[i*size2+j]-2*v[i*size2+j]+stemp2[i*size2+j])/(pow(dy,2))); 
//                    cout << "v is: " <<  v[i*Nyloc+j] << endl; 
//                        
//        }
//    }
//}
//
//void LidDrivenCavityExp::NextInnerVorticity(int i,int imax, int j,int jmax, int size1, int size2){
//    for (i;i<imax;i++){
//         for (j;j<jmax;j++){
//            temp1=(stemp3[i*size2+j]-stemp2[i*size2+j])*(vtemp1[i*size2+j]-vtemp[i*size2+j])/(4*dy*dx);
//            temp2=(stemp1[i*size2+j]-stemp[i*size2+j])*(vtemp3[i*size2+j]-vtemp2[i*size2+j])/(4*dy*dx); 
//            temp3=(1/Re)*((vtemp2[i*size2+j]-2*v[i*Nyloc+j]+vtemp3)/(pow(dy,2))+(vtemp-2*v[i*Nyloc+j]+vtemp1)/(pow(dx,2))); 
//            vnew[i*Nyloc+j]=v[i*Nyloc+j]+dt*(temp2-temp1+temp3); 
//         }
//    }
//    
//}



void LidDrivenCavityExp::InnerVorticity(int* ind, int size1, int size2){
//    int i, imax, j, jmax; 
//    
  int i, j; 
  int imax=ind[1];  
  int jmax=ind[3];  
    
    
    double stemp, stemp1, stemp2, stemp3; 
   
    
   
//        cout << "size1 is: " << size1 << endl; 
//        cout << "size2 is: " << size2 << endl; 
//        cout << "Index matrix, storing i,imax,j and jmax: " << endl; 
//        cout << i << " " << imax << " " << j << " " << jmax << endl << endl;   
//        if (rank==0){
//        cout << "Nxloc and Nyloc are: " << Nxloc << ", " << Nyloc << endl;
//        cout << "i, imax, j and jmax are" << ind[0] << ", " << ind[1] << ", " << ind[2] << ", " << ind[3] << endl; 
//        cout << "stored y_left is" << endl; 
//                    printmat(1,Nxloc,y_left); 
//        }
      for (i=ind[0];i<imax;i++){
//        cout << "i works here" << endl; 
            for (j=ind[2];j<jmax;j++){
                if (rank==0){
                cout << "i and j are: " << i << ", " << j << endl; 
                }
//                cout << "j works here" << endl; 
                if (i==0){
                    stemp=y_right[j]; 
                }
                else{
                    stemp=s[(i-1)*Nyloc+j]; 
                }
                if (i==Nxloc-1){
                    stemp1=y_left[j]; 
                }
                else {
                    stemp1=s[(i+1)*Nyloc+j];
                }
                if (j==0) {
                    stemp2=y_top[i];
                }
                else {
                    stemp2=s[i*Nyloc+j-1];
                }
                if (j==Nyloc-1){
                    stemp3=y_bot[i];
                }
                else {
                    stemp3=s[i*Nyloc+j+1]; //After checking all the conditions, v is finally computed 
                }
//                if (rank==0){
//                    cout << "stemp, stemp1, stemp 2 and stemp 3 are: " << stemp << ", " << stemp1 << ", " << stemp2 << ", " <<stemp3 << endl;
//                    
//                }
                v[i*Nyloc+j]=-((stemp1-2*s[i*Nyloc+j]+stemp)/(pow(dx,2))+(stemp3-2*s[i*Nyloc+j]+stemp2)/(pow(dy,2))); 
                  //  cout << "v is: " <<  v[i*size2+j] << endl; 
                    
                    
                }
        }
//    
    if (rank==5){
        cout << "Inner vorticity matrix is: " << endl; 
        printmat(Nyloc,Nxloc,v); 
//        cout << "Streamfunction is" << endl; 
//        printmat(Nyloc,Nxloc,s); 
    }
 }
// 
//
//
 void LidDrivenCavityExp::NextInnerVorticity(int* ind, int size1, int size2){
     
 LidDrivenCavityExp::BoundaryVectors();     
 LidDrivenCavityExp::CommunicateBound(); 
     
    int i, imax, j, jmax; 
    

    imax=ind[1];
    jmax=ind[3]; 
    
    
    double stemp, stemp1, stemp2, stemp3; 
    double vtemp,vtemp1,vtemp2,vtemp3; 
    double temp1,temp2,temp3; 
    for (i=ind[0];i<imax;i++){
        for (j=ind[2];j<jmax;j++){
            if (i==0){
                stemp=y_right[j]; 
                vtemp=x_right[j]; 
            }
            else {
                stemp=s[(i-1)*Nyloc+j]; 
                vtemp=v[(i-1)*Nyloc+j];
            }
            if (i==Nxloc-1){
                stemp1=y_left[j]; 
                vtemp1=x_left[j];
            }
            else {
                stemp1=s[(i+1)*Nyloc+j];
                vtemp1=v[(i+1)*Nyloc+j];
            }
            if (j==0) {
                stemp2=y_top[i];
                vtemp2=x_top[i];
            }
            else {
                stemp2=s[i*Nyloc+j-1];
                vtemp2=v[i*Nyloc+j-1];
            }
            if (j==Nyloc-1){
                stemp3=y_bot[i];
                vtemp3=x_bot[i]; 
            }
            else {
                stemp3=s[i*Nyloc+j+1]; //After checking all the conditions, v is finally computed 
                vtemp3=v[i*Nyloc+j+1];
            }
            temp1=(stemp3-stemp2)*(vtemp1-vtemp)/(4*dy*dx);
            temp2=(stemp1-stemp)*(vtemp3-vtemp2)/(4*dy*dx); 
            temp3=(1/Re)*((vtemp2-2*v[i*Nyloc+j]+vtemp3)/(pow(dy,2))+(vtemp-2*v[i*Nyloc+j]+vtemp1)/(pow(dx,2))); 
            vnew[i*Nyloc+j]=v[i*Nyloc+j]+dt*(temp2-temp1+temp3); 
//            if (rank==4){
//                cout << "temp3 is: " << vtemp3 <<  endl; 
//            }
            
        }
    }
//    if (rank==5){
////        cout << "Temp1 is: " << temp1 << endl; 
////        cout << "Temp2 is: " << << endl; 
////        cout << "Temp3 is: " << endl;
//        cout << "Next timestep Inner vorticity matrix is: " << endl; 
//        printmat(Nyloc,Nxloc,vnew); 
    }
 }
 
 void LidDrivenCavityExp::UpdateVorticity(int* ind,int size1,int size2){
     
     for (i=ind[0];i<ind[1];i++){
        for (j=ind[2];j<ind[3];j++){
            v[i*size2+j]=vnew[i*size2+j]; 
            
        }
     }
     
     
 }
 
//Moment of truth...

void LidDrivenCavityExp::Jacobi(){
    
    
}



void LidDrivenCavityExp::Integrate()
{    
  //PoissonSolver* psolver=new PoissonSolver(); //Create new instance psolver of class Poisson Solver to implement solver calculations  

  //psolver->Initialise(Nx,Ny,dx,dy); //Initialises values within PoissonSolver class 
  double t=0; //First time step value  
  
  int* ind=new int[4]; //Required for using GetIndex()
  
  
  ind=LidDrivenCavityExp::GetIndex(); 
  
  int i=ind[0]; 
  int imax=ind[1]; 
  int j=ind[2]; 
  int jmax=ind[3]; 
  
              
  int size1=imax-i; 
  int size2=jmax-j; 
//  double* stemp=new double[size1*size2]; //Initialise temp storage variables 
//  double* stemp1=new double[size1*size2]; 
//  double* stemp2=new double[size1*size2]; 
//  double* stemp3=new double[size1*size2];
//
//  double* vtemp=new double[size1*size2]; 
//  double* vtemp1=new double[size1*size2]; 
//  double* vtemp2=new double[size1*size2];  
//  double* vtemp3=new double[size1*size2];  
  
   // while (t<T){
//            cout << "Time step is: " << t << endl << endl; 
             LidDrivenCavityExp::BoundaryVectors(); //Obtain rows and columns located at the boundaries of subdomains 
             
             LidDrivenCavityExp::CommunicateBound(); 
            
             LidDrivenCavityExp::BoundaryConditions(); //Update with BCs //
             
           //  LidDrivenCavityExp::Updatetemps(i,imax,j,jmax,size1,size2
             
             LidDrivenCavityExp::InnerVorticity(ind,size1,size2); 
             
             
             LidDrivenCavityExp::NextInnerVorticity(ind,size1,size2); 
             
             LidDrivenCavityExp::UpdateVorticity(); //Must be updated for the new timestep 
             
             
            
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
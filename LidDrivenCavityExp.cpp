#include "LidDrivenCavityExp.h"
#include "PoissonSolverExp.h"
#include <iostream>
#include <math.h>
//#include "PrintMat.h"
#include <cblas.h>
#include <mpi.h> 
#include <cassert>
//#include "Communicate.h"


using namespace std;
 
 extern "C" {
  void Cblacs_pinfo(int*, int*);
  void Cblacs_get(int, int, int*);
  void Cblacs_gridinit(int*, char const*, int, int);
  void Cblacs_gridinfo(int, int*, int*, int*, int*);
  
  int numroc_(int const& n, int const& nb, int const& iproc, int const& isproc, int const& nprocs);
}
 
 //Create default destructor and constructor 
LidDrivenCavityExp::LidDrivenCavityExp()=default; //Constructor 


LidDrivenCavityExp::~LidDrivenCavityExp()=default;  //Destructor 

int LidDrivenCavityExp::CheckParallel() { //Check if parallelisation successful 
     
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

void LidDrivenCavityExp::AssignGlobal(string *val){ //Assign global values based on user input info 

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
     
     row_loc=Py; //No. subdomain rows and columns expected due to partitioning in process 
     col_loc=Px; 
     
     Cblacs_pinfo(&rank,&nprocs); 
     assert(nprocs>=row_loc*col_loc); //Assertion statement to validate user input np for chosen partition size px and py 
     
     Cblacs_get(0,0,&ctxt); //Attains system context 
     Cblacs_gridinit(&ctxt,"Col-major",row_loc,col_loc); //Creates process grid for a given context 
     
     Cblacs_gridinfo(ctxt,&row_loc,&col_loc,&myrow,&mycol); 
}


void LidDrivenCavityExp::SetDomainSize()
{  
   Lxloc=Lx/Px; 
   Lyloc=Ly/Py; 
}



void LidDrivenCavityExp::SetGridSize() //Determine no. of rows and columns of each subdomain based on its position on the global grid using numroc 
{
      Nyloc=numroc_(Ny,1,myrow,0,row_loc);
      Nxloc=numroc_(Nx,1,mycol,0,col_loc); 
     
}


void LidDrivenCavityExp::Initialise() //Initialise class variables 
{    
      
      v=new double[Nxloc*Nyloc]; //Initialize vorticity values 
      vnew=new double[Nxloc*Nyloc];//Voriticity at new timestep value 
      s=new double[Nxloc*Nyloc]; //Initialize streamfunction value 

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
      y_top=new double[Nxloc]; //Row vector storing streamfunction values sent from bottom row (storing top wall conditions) of the subdomain above at the boundary
      x_top=new double[Nxloc]; //Row vector storing vorticity values sent from bottom row of subdomain right above. 
      
      v_bot=new double[Nxloc]; 
      s_bot=new double[Nxloc]; 
      y_bot=new double[Nxloc]; //Row vector storing streamfunction values sent from top row (storing bottom wall conditions) of subdomain right below. 
      x_bot=new double[Nxloc]; //Row vector storing vorticity values sent from top row of adjacent subdomain right below. 
 
 
      for (int i=0;i<Nxloc*Nyloc;i++){
              v[i]=0; 
              s[i]=0; 
              vnew[i]=0; //Initialise all values to be zero 
       }

      
        dx=Lx/(Nx-1.0); //Set global variables dx, dy and U
        dy=Lx/(Ny-1.0); 
        U=1.0;
        
        if (dt>=Re*dx*dy/4){ //Ensure timestep value satosfoes restriction
           throw logic_error("Chosen value of dt too large"); 
    
        }
      
}


void LidDrivenCavityExp::BoundaryVectorsGen(double* arr, double* arr_left, double* arr_right, double* arr_bot, double* arr_top){
 
    for (int i=0;i<Nyloc;i++){
    arr_left[i]=arr[i]; 
    arr_right[i]=arr[Nyloc*(Nxloc-1)];
    }
    
    for (int i=0;i<Nxloc;i++){
    arr_bot[i]=arr[i*Nyloc];
    arr_top[i]=arr[(i+1)*(Nyloc-1)]; 
    }
    
    
}



void LidDrivenCavityExp::Communicate(double* arr_top,double* buf_top, double*arr_bot, double* buf_bot, double* arr_left, double* buf_left, double* arr_right, double* buf_right){
    
    if (Px>1){
        if (mycol<Px-1){
                MPI_Send(arr_right,Nyloc,MPI_DOUBLE,rank+Py,0,MPI_COMM_WORLD);
                MPI_Recv(buf_left,Nyloc,MPI_DOUBLE,rank+Py,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 

         }
       if (mycol>0){
               MPI_Recv(buf_right,Nyloc,MPI_DOUBLE,rank-Py,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
               MPI_Send(arr_left,Nyloc,MPI_DOUBLE,rank-Py,1,MPI_COMM_WORLD); 
               
       }
    }
//  
   if (Py>1){
       if (myrow<Py-1){
              MPI_Send(arr_top,Nxloc,MPI_DOUBLE,rank+1,2,MPI_COMM_WORLD); 
              MPI_Recv(buf_bot,Nxloc,MPI_DOUBLE,rank+1,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
    
       }
       if (myrow>0){
              MPI_Recv(buf_top,Nxloc,MPI_DOUBLE,rank-1,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
              MPI_Send(arr_bot,Nxloc,MPI_DOUBLE,rank-1,3,MPI_COMM_WORLD); 
     
      }

  }
  
  MPI_Barrier(MPI_COMM_WORLD); 
  

}
    

 void LidDrivenCavityExp::BoundaryConditions(){
      
     LidDrivenCavityExp::BoundaryVectorsGen(s,s_left,s_right,s_bot,s_top); 
     LidDrivenCavityExp::Communicate(s_top,y_top,s_bot,y_bot,s_left,y_left,s_right,y_right); 
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
                            v[i]=(s[i]-y_top[i])*(2.0/pow(dy,2))-(2*U/dy);  
                        }
            }
             else {
                  for (int i=0;i<Nxloc;i++){
                            v[(i+1)*(Nyloc-1)+i]=(2/pow(dy,2))*(s[(i+1)*(Nyloc-1)+i]-s[(i+1)*(Nyloc-1)+i-1])-(2*U/dy); //Top
                        }
                        
                  }
         }
        
        if  (mycol==0){ //left 
             if (Nxloc==1){ //Communication necessary, call communicator class function
                     for (int i=0;i<Nyloc;i++){
                        v[i]=(s[i]-y_right[i])*(2.0/pow(dx,2)); //Stores s_left from adjacent array
                    }
                }
                
              else {
                   for (int i=0;i<Nyloc;i++){
                    v[i]=(2/(pow(dx,2)))*(s[i]-s[i+Nyloc]); 
                    }
              }
         }

     
      if  (mycol==Px-1){
             if (Nxloc==1){ //Communication necessary, call communicator class function
                     for (int i=0;i<Nyloc;i++){
                        v[i]=(s[i]-y_left[i])*(2.0/pow(dx,2)); 
                    }

                }

              else {
               for (int i=0;i<Nyloc;i++){
                v[Nyloc*(Nxloc-1)+i]=(2/pow(dx,2))*(s[(Nxloc-1)*Nyloc+i]-s[(Nxloc-2)*Nyloc+i]); //Right 
                }
              }
         }
    
 
 }

int* LidDrivenCavityExp::GetIndex(){  //Returns an array of indices consisting of i,imax,j,jmax, etc. to be used for calculating inner vorticity values 
  int* ind = new int[4]; 
         if (myrow==0 && mycol==0){ //Leftmost corner, bottom wall
             ind[0]=1; 
             ind[1]=Nxloc; 
             ind[2]=1; 
             ind[3]=Nyloc;
         }
         else if (myrow==Py-1&&mycol==0){ //Leftmost corner, top wall 
             ind[0]=1; 
             ind[1]=Nxloc; 
             ind[2]=0; 
             ind[3]=Nyloc-1;
         }
         else if (myrow==Py-1&&mycol==Px-1){ //Rightmost corner, top wall 
             ind[0]=0; 
             ind[1]=Nxloc-1; 
             ind[2]=0; 
             ind[3]=Nyloc-1;
         }
         else if (mycol==Px-1&&myrow==0){ //Rightmost corner, bottom wall 
             ind[0]=0; 
             ind[1]=Nxloc-1; 
             ind[2]=1; 
             ind[3]=Nyloc;
         }
         else if (mycol==0) { //General Left wall 
            ind[0]=1; 
             ind[1]=Nxloc; 
             ind[2]=0; 
             ind[3]=Nyloc; 
         }
         else if (myrow==0){ //General bottom wall 
             ind[0]=0; 
             ind[1]=Nxloc; 
             ind[2]=1; 
             ind[3]=Nyloc;
         }
         else if (myrow==Py-1){ //General top wall 
             ind[0]=0; 
             ind[1]=Nxloc; 
             ind[2]=0; 
             ind[3]=Nyloc-1;
         }
         else if (mycol==Px-1){ //General Right wall 
             ind[0]=0; 
             ind[1]=Nxloc-1; 
             ind[2]=0; 
             ind[3]=Nyloc;
         }
         else {              //Interior subdomain points 
             ind[0]=0; 
             ind[1]=Nxloc; 
             ind[2]=0; 
             ind[3]=Nyloc;
             
         }
         return ind; 
         
         delete[] ind; 

}





void LidDrivenCavityExp::InnerVorticity(int* ind, int size1, int size2){
  
  int i, j; 
  int imax=ind[1];  
  int jmax=ind[3];  
    
    
    double stemp, stemp1, stemp2, stemp3; //temp variables used for storage 
   
      for (i=ind[0];i<imax;i++){
            for (j=ind[2];j<jmax;j++){
                if (i==0){
                    stemp=y_right[j]; //Conditional statements introduced to ensure communicated values are used at boundaries //Communication for s already conducted in BoundaryConditions
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
                v[i*Nyloc+j]=-((stemp1-2*s[i*Nyloc+j]+stemp)/(pow(dx,2))+(stemp3-2*s[i*Nyloc+j]+stemp2)/(pow(dy,2))); 
                    
                    
                }
        }

 }

void LidDrivenCavityExp::NextInnerVorticity(int* ind, int size1, int size2){
         
     LidDrivenCavityExp::BoundaryVectorsGen(v,v_left,v_right,v_bot,v_top); //Ensure boundary conditions at timestep t are assigned to boundary vectors    
     LidDrivenCavityExp::Communicate(v_top,x_top,v_bot,x_bot,v_left,x_left,v_right,x_right);  //Ensure boundary vectors are communicated to buffer arrays 
    // Communicate(v_top,x_top,v_bot,x_bot,v_left,x_left,v_right,x_right, Px, Py, myrow, mycol);
    
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
                stemp3=s[i*Nyloc+j+1]; 
                vtemp3=v[i*Nyloc+j+1];
            }
            temp1=(stemp3-stemp2)*(vtemp1-vtemp)/(4*dy*dx);
            temp2=(stemp1-stemp)*(vtemp3-vtemp2)/(4*dy*dx); 
            temp3=(1/Re)*((vtemp2-2*v[i*Nyloc+j]+vtemp3)/(pow(dy,2))+(vtemp-2*v[i*Nyloc+j]+vtemp1)/(pow(dx,2))); 
            vnew[i*Nyloc+j]=v[i*Nyloc+j]+dt*(temp2-temp1+temp3); 
            
        }
    }
 }
 
 void LidDrivenCavityExp::UpdateVorticity(int* ind,int size1,int size2){ //Ensures interior vorticity values are updated so that v=vnew for next timestep 
     
     for (int i=ind[0];i<ind[1];i++){
        for (int j=ind[2];j<ind[3];j++){
            v[i*Nyloc+j]=vnew[i*Nyloc+j]; 
            
        }
     }
 }


double* LidDrivenCavityExp::Jacobi(double*v, double* s, int* ind,int size1, int size2, int Nyloc, int Nxloc){
    double* snew=new double[Nyloc*Nxloc];
    double* snew1=new double[Nyloc*Nxloc]; //Dummy variable 
    double* diff=new double[Nyloc*Nxloc]; //To calculate difference between snew and snew1 
        for (int i=0;i<Nxloc*Nyloc;i++){
            snew[i]=s[i];  //Initialize before iterating, initial guessed value selected as zero at each timestep 
            snew1[i]=0; 
        }

    double stemp, stemp1, stemp2, stemp3; 
    double alph=(pow(dx,2)*pow(dy,2))/(2*pow(dy,2)+2*pow(dx,2));
    
    double* snew_left=new double[Nyloc]; 
    double* snew_right=new double[Nyloc]; 
    double* snew_top=new double[Nxloc]; 
    double* snew_bot=new double[Nxloc]; 
    
    double* ynew_left=new double[Nyloc]; 
    double* ynew_right=new double[Nyloc]; 
    double* ynew_top=new double[Nxloc]; 
    double* ynew_bot=new double[Nxloc]; 
    
    int count=1; 
    
    double err; 

    while(count<100) {
          LidDrivenCavityExp::BoundaryVectorsGen(snew,snew_left,snew_right,snew_bot,snew_top); //Will assign values to boundary vectors. snew changes with every iteration 
          LidDrivenCavityExp::Communicate(snew_top,ynew_top,snew_bot,ynew_bot,snew_left,ynew_left,snew_right,ynew_right); //Communicate new values
            for (int i=ind[0];i<ind[1];i++){
                for (int j=ind[2];j<ind[3];j++){
                        if (i==0){
                            stemp=ynew_right[j]; 
                        }
                        else{
                            stemp=snew[(i-1)*Nyloc+j]; 
                        }
                        if (i==Nxloc-1){
                            stemp1=ynew_left[j]; 
                        }
                        else {
                            stemp1=snew[(i+1)*Nyloc+j];
                        }
                        if (j==0) {
                            stemp2=ynew_top[i];
                        }
                        else {
                            stemp2=snew[i*Nyloc+j-1];
                        }
                        if (j==Nyloc-1){
                            stemp3=ynew_bot[i];
                        }
                        else {
                            stemp3=snew[i*Nyloc+j+1]; //After checking all the conditions, v is finally computed 
                        }
                        snew1[i*Nyloc+j]=alph*(v[i*Nyloc+j]+((stemp1+stemp)/pow(dx,2))+((stemp3+stemp2)/pow(dy,2)));
                }
            }

        cblas_dcopy(Nxloc*Nyloc,snew1,1,diff,1); //Copy so that diff=snew1
        cblas_daxpy(Nxloc*Nyloc,-1,snew,1,diff,1); //diff=diff-snew or snew1-snew
         err=cblas_dnrm2(Nxloc*Nyloc,diff,1); //Calculate norm error
        cblas_dcopy(Nxloc*Nyloc,snew1,1,snew,1); //Update so that snew=snew1 
//        if (rank==5){
//        cout << "snew1 is: " << endl; 
//        printmat(Nyloc,Nxloc,snew1); 
//        cout << "err is: " << err << endl; 
//        }
        count++;
    }

    
    return snew; 
    
    
    delete[] snew;
    delete[] snew1; 
    delete[] diff;  
    delete[] ynew_left;
    delete[] ynew_right; 
    delete[] ynew_top; 
    delete[] ynew_bot; 
    delete[] snew_left;
    delete[] snew_right; 
    delete[] snew_top; 
    delete[] snew_bot; 
}

void LidDrivenCavityExp::OutputValues(){ //Also used for parallelisation 

    int var=Nxloc*Nyloc; 
    double* sglobal=nullptr; //recvbuf, global array storing values from all processes 
    int* temp=nullptr; //temp array, used for storing subdomain array size for each process 
    int* displs=nullptr; //array showing placement of recevied values with respect to sglobal in MPI_Gatherv

      if (rank==0){
        sglobal=new double[Nx*Ny]; //Initialise array size, values will be gathered at root process 
        temp=new int[nprocs]; 
        displs=new int[nprocs]; 
      }


    MPI_Gather(&var,1,MPI_INT,temp,1,MPI_INT,0,MPI_COMM_WORLD); //Obtain values for array temp 
//
    if (rank==0){
        displs[0]=0; 
        for (int i=1;i<nprocs;i++){
            displs[i]=displs[i-1]+temp[i-1];  
        }
      
    }

    MPI_Gatherv(s,var,MPI_DOUBLE,sglobal,temp,displs,MPI_DOUBLE,0,MPI_COMM_WORLD); //Used gatherv as subdomain size varies for each process 
    //    
        if (rank==0){
            cout << "sglobal is: " << endl; 
            printmat(Ny*Nx,1,sglobal);
        }
    }



void LidDrivenCavityExp::Integrate()
{    
      PoissonSolverExp* psolver=new PoissonSolverExp(); //Create new instance psolver of class Poisson Solver to implement solver calculations  
 
      psolver->Initialise(Nxloc,Nyloc,Px,Py,dx,dy,myrow,mycol,rank); 
      
      double* snew=new double[Nxloc*Nyloc];
      double t=0; //First time step value  
      
      int* ind=new int[4]; //Required for using GetIndex()
      
      ind=LidDrivenCavityExp::GetIndex(); //Obtain indices for each process 
      
      int i=ind[0]; 
      int imax=ind[1]; 
      int j=ind[2]; 
      int jmax=ind[3]; 
      
                  
      int size1=imax-i; 
      int size2=jmax-j;  
  
       while (t<T){
//                 if (rank==0){
//                cout << "Time step is: " << t << endl << endl; 
//                
//                 }
                 //LidDrivenCavityExp::BoundaryVectors(); //Obtain rows and columns located at the boundaries of subdomains 
                 
                // LidDrivenCavityExp::BoundaryVectorsGen(s,s_left,s_right,s_bot,s_top); 
                
                 LidDrivenCavityExp::BoundaryConditions(); //Update with BCs // Communicator function already defined within BCs 
                 
                 LidDrivenCavityExp::InnerVorticity(ind,size1,size2); 
                 
                 LidDrivenCavityExp::NextInnerVorticity(ind,size1,size2); //Communicator function already called 
                 
                 LidDrivenCavityExp::UpdateVorticity(ind,size1,size2); //Must be updated for the new timestep 
                 
                 for (int i=0;i<Nxloc*Nyloc;i++){
                     snew[i]=s[i]; 
                 }
                 
                 s=LidDrivenCavityExp::Jacobi(v,s,ind,size1,size2,Nyloc,Nxloc); //Determine new s 
                 //s=psolver->Jacobi(ind,s,v,snew); //Results in a segmentation fault somewhere 
                            
                 t+=dt;
                 
           }
             
    }
    


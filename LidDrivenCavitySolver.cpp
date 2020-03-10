#include <iostream>
using namespace std;

#include "LidDrivenCavity.h"

int main(int argc, char **argv)
{
    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();
    string inputs[argc-1]={"Lx","Ly","Nx","Ny","Px","Py","dt","T","Re"}; 
    
     int n=argc; //Record no. of inputs 
     double Ly=stod(argv[2]); 
     int Nx=stoi(argv[4]); 
     int Ny=stoi(argv[6]); 
     int Px=stoi(argv[8]); 
     int Py=stoi(argv[10]); 
     double dt=stod(argv[12]);
     double T=stoi(argv[14]); 
     double Re=stod(argv[16]); 
     
     for(int i=0;i<argc/2;i++){
         cout << inputs[i] << ": " << argv[2*(i+1)] << endl;
     }
     
     
//    double T; 
//    double Re;  
//    for(int i=1; i<argc; i++){



    // Configure the solver here...
    // ...
//    solver->Initialise();
//
//    // Run the solver
//    solver->Integrate();

	return 0;
}


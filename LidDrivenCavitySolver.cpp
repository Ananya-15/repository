#include <iostream>
#include <iomanip> 
using namespace std;

#include "LidDrivenCavity.h"
//#include "PrintMat.h"
//#include "PoissonSolver.h"


int main(int argc, char **argv)
{

   
    
    string inputs[(argc-1)/2]={"Lx","Ly","Nx","Ny","Px","Py","dt","T","Re"};

    cout << "Argument length is: " << argc << endl; 
    int n=argc/2; //Record no. of inputs 
    string *val=new string [n]; //Create array values, storing all the input values, to be passed as inputs to class function 
    for (int i=0;i<n;i++){
        val[i]=argv[2*(i+1)]; 
    }
     LidDrivenCavity* solver = new LidDrivenCavity();     // Create a new instance of the LidDrivenCavity class
     solver->Initialise(val); 


    // Configure the solver here...
    // ...
//    solver->Initialise();
//
//    // Run the solver
    solver->Integrate();

	return 0;
}


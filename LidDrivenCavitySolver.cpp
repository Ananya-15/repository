#include <iostream>
using namespace std;

#include "LidDrivenCavity.h"

void printmat(int N, int M, string* A){ //matrix printing algorithm 
for (int i = 0; i < N; i++){ //Array printing algo 
           for (int j = 0; j < M; j++){
            cout << A[i*M+j] << " ";
            }
          cout << endl;
         }
}

int main(int argc, char **argv)
{
    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();
    
    string inputs[(argc-1)/2]={"Lx","Ly","Nx","Ny","Px","Py","dt","T","Re"};

    
    int n=argc/2; //Record no. of inputs 
    string *val=new string [n]; //Create array values, storing all the input values, to be passed as inputs to class function 
    for (int i=0;i<n;i++){
        val[i]=argv[2*(i+1)]; 
    }
     
     solver->Initialise(val); 
     
     


    // Configure the solver here...
    // ...
//    solver->Initialise();
//
//    // Run the solver
//    solver->Integrate();

	return 0;
}


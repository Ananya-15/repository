#include <iostream>
#include <iomanip>
using namespace std; 


void printmat(int N, int M, double* A){ //matrix printing algorithm 
for (int i = 0; i < N; i++){ //Array printing algo 
           for (int j = 0; j < M; j++){
            cout << setw(7) << setprecision(5) << A[j*N+i] << " ";
            }
          cout << endl;
         }
}
default: myprog 

LidDrivenCavitySolver.o: LidDrivenCavitySolver.cpp
	g++ -Wall -o LidDrivenCavitySolver.o -c LidDrivenCavitySolver.cpp

LidDrivenCavity.o: LidDrivenCavity.cpp LidDrivenCavity.h PrintMat.h PoissonSolver.h
	g++ -Wall -o LidDrivenCavity.o -c LidDrivenCavity.cpp

myprog: LidDrivenCavitySolver.o LidDrivenCavity.o 
	g++ -o myprog LidDrivenCavitySolver.o LidDrivenCavity.o -llapack -lblas

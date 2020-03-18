default: myprog 

LidDrivenCavitySolver.o: LidDrivenCavitySolver.cpp
	g++ -Wall -o LidDrivenCavitySolver.o -c LidDrivenCavitySolver.cpp

LidDrivenCavity.o: test.cpp LidDrivenCavity.h PrintMat.h PoissonSolver.h
	g++ -Wall -o test.o -c test.cpp

myprog: LidDrivenCavitySolver.o test.o 
	g++ -o myprog LidDrivenCavitySolver.o test.o -llapack -lblas

.PHONY: clean 
	target
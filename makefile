default: myprogpar

Paralleltry2.o: Paralleltry2.cpp
	mpicxx -Wall -o Paralleltry2.o -c Paralleltry2.cpp

LidDrivenCavityExp.o: LidDrivenCavityExp.cpp LidDrivenCavity.h PrintMat.h 
	mpicxx -Wall -o LidDrivenCavityExp.o -c LidDrivenCavityExp.cpp

myprogpar: Paralleltry2.o LidDrivenCavityExp.o 
	mpicxx -o myprogpar Paralleltry2.o LidDrivenCavityExp.o -lscalapack-openmpi -lblas

.PHONY: clean 
	target
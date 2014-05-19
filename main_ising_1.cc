#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <omp.h>

using namespace std;

/* compile with "g++ -fopenmp main_ising_1.cc lattice.cc -o ising_ex -lgsl -lgslcblas" */

#include "lattice.h"	// lattice class

int main(int argc, char *argv[])
{	
	if(argc != 5){
		cout<<"usage : ising_ex <length> <dimension> <B-field> <iterations>"<<endl;
		return 1;
	}
	
	int Len = atoi(argv[1]);
	int dim = atoi(argv[2]);
	double Bfield = atof(argv[3]);
	int iterations = atoi(argv[4]);
	

	lattice l1(Len,dim,Bfield,iterations);		// lattice(L,d,B,iter);
	
	l1.hot_start();

	l1.betarun();
	
	return 0;
}


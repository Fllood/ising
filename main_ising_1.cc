#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <fstream>


using namespace std;


#include "lattice.h"	// lattice class

int main(int argc, char *argv[])
{	
	if(argc != 7){
		cout<<"usage : ising_ex <length> <dimension> <B-field> <iterations> <T> <eq_time>"<<endl;
		return 1;
	}
	
	int Len = atoi(argv[1]);
	int dim = atoi(argv[2]);
	double Bfield = atof(argv[3]);
	int iterations = atoi(argv[4]);
	double Temp = atof(argv[5]);
	int eq_time = atoi(argv[6]);

	lattice l1(Len,dim,Bfield,iterations,Temp,eq_time);		// lattice(L,d,B,iter);
	
	
	l1.hot_start();

	l1.run();
	
	return 0;
}


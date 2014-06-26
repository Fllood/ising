#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <string>


using namespace std;


#include "lattice.h"	// lattice class


int main(int argc, char *argv[])
{	
	
	if(argc != 9){
		cout<<"usage : ising_ex <length> <dimension> <B-field> <iterations> <T> <eq_time> <algorithm> <output>"<<endl;
		return 1;
	}
	
	int Len = atoi(argv[1]);
	int dim = atoi(argv[2]);
	double Bfield = atof(argv[3]);
	int iterations = atoi(argv[4]);
	double Temp = atof(argv[5]);
	int eq_time = atoi(argv[6]);
	string mode_for_sweep = argv[7];
	string output_mode = argv[8];

	lattice l1(Len,dim,Bfield,iterations,Temp,eq_time, mode_for_sweep, output_mode);	
	
	
	if(Temp){		//skip when Temp = 0
		//l1.one_temp();
		cout<<"distance: "<<l1.dist(0,7)<<endl;
	}
	else{
		l1.scan_t();
		}
	
	return 0;
}


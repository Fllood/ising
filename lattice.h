#ifndef LATTICE_H
#define LATTICE_H

#include <vector>

/* GNU Scientific Library random number generator (rng) */
#include <gsl/gsl_rng.h>

using namespace std;

class lattice{
	private:
		vector<short> spins;		// Spins contained in one dim vector
		
		int L;						// Length of axis		
		
		int d;						// Dimensions	
		
		int V;						// Volume
		
		double b;					// Inverse temperature beta
		
		double T;					// Temperature
		
		double B;					// Magnetic field
		
		int iter;
		
		
		vector<double> lookup_J;	// Lookup table for (2*d+1) exponential values containing J
		
		vector<double> lookup_B;	// Lookup table for 2 exponential values containing B
		
		/* pointer to the rng */
		gsl_rng *rng;
		
		ofstream file;
		
		vector<double> betas;
	
	public:
		
		lattice(int length, int dim, double Bfield, int iter);
		
		void update_lookups();		
		
		void cold_start();
		
		void hot_start();
		
		void display();	
		
		int get_V();
		
		int get_nn_sum(int pos);
		
		double get_mag();
		
		double get_eng();
		
		void sweep();
		
		void run();
		
		void betarun();
	};

#endif
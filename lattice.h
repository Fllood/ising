#ifndef LATTICE_H
#define LATTICE_H

#include <vector>

/* GNU Scientific Library random number generator (rng) */
#include <gsl/gsl_rng.h>

using namespace std;

class lattice{
	private:
		vector<short> spins;		// Spins contained in one dim vector
		
		vector<double> mag;		// Magnetization measurements
		
		vector<double> eng;		// Energy measurements
		
		vector<double> corr_mag // Magnetic correlation function
		
		vector<double> corr_eng // Energy correlation function
		
		
		int L;						// Length of axis		
		
		int d;						// Dimensions	
		
		int V;						// Volume
		
		double b;					// Inverse temperature beta
		
		double T;					// Temperature
		
		double B;					// Magnetic field
		
		int iter;					// Number of steps
		
		int t_eq;					// equilibration time
		
		
		vector<double> lookup_J;	// Lookup table for (2*d+1) exponential values containing J
		
		vector<double> lookup_B;	// Lookup table for 2 exponential values containing B
		
		/* pointer to the rng */
		gsl_rng *rng;
		
		ofstream file;
		
		vector<double> betas;
	
	public:
		
		lattice(int length, int dim, double Bfield, int iter, double Temp, int eq_time);
		
		void update_lookups();		
		
		void cold_start();
		
		void hot_start();
		
		void set_T(double Temp);
		
		void display();	
		
		int get_V();
		
		int get_nn_sum(int pos);
		
		double get_mag();
		
		vector<double> get_mag_vec();
		
		double get_eng();
		
		vector<double> get_eng_vec();
		
		double corr_func(int t_c, const vector<double>& vec);
		
		void calc_corr_t(const vector<double>& vec, vector<double>& corr);
		
		void rem_equilib(vector<double>& vec);
		
		void sweep();
		
		void run();
		
		void betarun();
	};

#endif

#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <string>

/* GNU Scientific Library random number generator (rng) */
#include <gsl/gsl_rng.h>

using namespace std;

class lattice{
	private:
		vector<short> spins;		// Spins contained in one dim vector
		
		vector<double> mag;			// Magnetization measurements
		
		vector<double> mag_2;
		
		vector<double> eng;			// Energy measurements
		
		vector<double> eng_2;
		
		vector<double> cov_mag; 	// Magnetic covariance function
		
		vector<double> corr_mag;	// Magnetic correlation function
		
		vector<double> cov_eng; 	// Energy covariance function
		
		vector<double> corr_eng; 	// Energy correlation function
		
		double avg_mag;				// Average magnetization per spin
		
		double avg_mag_2;
		
		double avg_eng;				// Average energy per spin
		
		double avg_eng_2;
		
		int L;						// Length of axis		
		
		int d;						// Dimensions	
		
		int V;						// Volume
		
		double b;					// Inverse temperature beta
		
		double T;					// Temperature
		
		double B;					// Magnetic field
		
		int iter;					// Number of steps
		
		int t_eq;					// equilibration time
		
		
		vector<double> lookup_met_J;	// Lookup table for (2*d+1) exponential values containing J
		
		vector<double> lookup_met_B;	// Lookup table for 2 exponential values containing B
		
		vector<double> lookup_heat_J;	// Lookup table for (2*d+1) exponential values containing J (heat bath)
		
		string mode;		// metropolis or heat bath
		
		string output;		// choice of output
		
		/* pointer to the rng */
		gsl_rng *rng;
		
		ofstream file;
		
		vector<double> betas;
	
	public:
		
		lattice(int length, int dim, double Bfield, int iter, double Temp, int eq_time, string mode_for_sweep, string output_mode);
		
		void update_lookups();		
		
		void cold_start();
		
		void hot_start();
		
		void set_T(double Temp);
		
		void display();	
		
		int get_V();
		
		int get_nn_sum(int pos);
		
		double get_mag();
		
		double get_eng();
		
		double cov_func(int t_c, const vector<double>& vec);
		
		void calc_cov_t(const vector<double>& vec, vector<double>& corr);
		
		void calc_mag_cov();
		
		void calc_eng_cov();
		
		void calc_mag_corr();
		
		void calc_eng_corr();
		
		double calc_tau(const vector<double>& corr);
		
		double get_avg(const vector<double>& vec);
		
		double get_std_err(const vector<double>& cov, const vector<double>& corr);
		
		double get_spec_heat();
		
		double get_mag_sus();
		
		vector<double> get_vec(string choice);
		
		double get_val(string choice);
		
		void equilibrate();
		
		void sweep_met();
		
		void sweep_heat();
		
		void run();
		
	};

#endif

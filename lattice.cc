#include <stdlib.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>

#include "lattice.h"

#include <gsl/gsl_rng.h>

using namespace std;

lattice::lattice(int length, int dim, double Bfield, int iterations, double Temp, int eq_time){
	
	L = length;
	d = dim;
	V = pow(L,d);
	
	t_eq = eq_time;
	
	T = Temp;
	b = 1/T;

	B = Bfield;
	iter = iterations;
	
	spins.reserve(V);
	
	mag.reserve(iter);
	
	eng.reserve(iter);
	
	lookup_J.reserve(2*d+1);
	lookup_B.reserve(2);
	
	this->update_lookups();	
	
	// gsl function to initialize the rng
	rng = gsl_rng_alloc(gsl_rng_mt19937);
	
	// set seed based on current time
	long seed = time(NULL);	
	
	gsl_rng_set(rng,seed);	
	
	
	}

void lattice::update_lookups(){
	lookup_J.clear();
	lookup_B.clear();
	mag.clear();
	eng.clear();
	
	for(int q = -2*d; q <= 2*d; q += 2){
		lookup_J.push_back(exp(2*q*b));				// J = 1
	}
	
	for(int i = -1; i < 2; i += 2){					// i = +- 1
		lookup_B.push_back(exp(i*b*B));
	}
	}
	
void lattice::cold_start(){			//all spins pointing up
	for(int i = 0; i < V; i++){
		spins.push_back(1);		
		}	
	}

void lattice::hot_start(){			//spins randomly set
	spins.clear();
	for(int i = 0; i < V; i++){
		double rn = gsl_rng_uniform(rng);
		if(rn < 0.5){
			spins.push_back(1);
			}
		else{
			spins.push_back(-1);	
			}	
		}	
	}

int lattice::get_V(){
	return V;	
	}

int lattice::get_nn_sum(int pos){		// works for d dimensions
	int nn,sum = 0,div;
	
	for(int i = d; i > 0; i--){
		div = V/(pow(L,i));
		
		if((nn=pos+div)>=V)nn -= V;			// helical boundary conditions
		sum += spins[nn];
		if((nn=pos-div)<0)nn += V;
		sum += spins[nn];
	}
	
	return sum;
	}

double lattice::get_mag(){
	double sum = 0;
	
	for(int i = 0; i < V; i++){
		sum += spins[i];
		}
	return fabs(sum/((double)V));
	}

vector<double> lattice::get_mag_vec(){
	return mag;	
	}

double lattice::get_eng(){
	double sum = 0;
	for(int i = 0; i < V; i++){
		sum	-= this->get_nn_sum(i)*spins[i];
		sum -= B*spins[i];
		}
	return sum/((double)V);
	}

vector<double> lattice::get_eng_vec(){
	return eng;	
	}

void lattice::set_T(double Temp){
	T = Temp;
	b = 1/T;
	}

void lattice::sweep(){
	int s_j_i, q, s_j;							// Index of proposed spin flip, sum of nearest neighbors * value(s_j), value of flipped spin, sum of nearest neighbors
	double dE, rho, rn;									// Energy difference, acceptance probab. rho, random number
	for(int i = 0; i < V; i++){
		s_j_i = floor(V * gsl_rng_uniform(rng));
		s_j = spins.at(s_j_i);
		
		q = this->get_nn_sum(s_j_i);
		dE = 2*s_j*q + 2*B*s_j;	
		//cout<<"dE = "<<dE<<endl;							// J = 1;
		if(dE <= 0) spins.at(s_j_i) = -s_j;				// If dE <= 0 accept flip
		else{
			//cout<<"q = "<<q<<endl;
			//cout<<"lookup_index_J = "<<(-s_j*q/2)+d<<endl;
			
			rho = lookup_J.at(((-s_j*q)/2)+d);					// Build up rho
			if(s_j < 0) rho *= lookup_B.at(0);
			else rho *= lookup_B.at(1);
			//cout<<"rho = "<<rho<<endl;
			
			rn = gsl_rng_uniform(rng);
			if(rn <= rho) spins.at(s_j_i) = -s_j;		// Accept flip with prob. rho
		}
	}
}

void lattice::run(){
	ostringstream fs;	
	
	this->update_lookups();
	
	
	fs<<"data/ising_"<<T<<".dat";
		
	file.open(fs.str().c_str());
	file<<"#t mag eng# at L="<<L<<" d="<<d<<endl;
	file.precision(10);
	cout<<"T = "<<T<<endl;
	file<<"#T = "<<T<<endl;
		
	for(int t = 0; t<iter; t++){
		this->sweep();
		file<<t<<" "<<this->get_mag()<<" "<<this->get_eng()<<endl;	// Write measurements to file
		
		mag.push_back(this->get_mag());										// Save measurements in vector
		eng.push_back(this->get_eng());	
		}
	
	file.close();	
	
	}

void lattice::betarun(){
	for(int i = 0; i < betas.size(); i++){
		b = betas.at(i);
		T = 1/betas.at(i);
			
		}	
	}

double lattice::corr_func(int t_c, const vector<double>& y){
	double sum1 = 0, sum2 = 0, sum3 = 0;
	for(int i = 1; i <= (iter - t); i++){
		sum1 += 	y.at(i)*y.at(i+t);
		sum2 +=  y.at(i);
		sum3 +=  y.at(i+t);	
		}
	double norm = 1/(double(iter-t));
	sum1 *= norm;
	sum2 *= norm;
	sum3 *= norm;
	
	return sum1 - sum2 * sum3;
	}

void lattice::display(){
	for(int i = 0; i < V; i++){
		if(spins.at(i)>0)cout<<"  ";
		else cout<<" ";
		cout<<spins.at(i);
		if( (i+1) % L == 0){
			cout<<endl;			
			}
		if( (i+1) % (L*L) == 0){
			cout<<endl;			
			}
		}	
	}

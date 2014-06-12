#include <stdlib.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include "lattice.h"

#include <gsl/gsl_rng.h>

using namespace std;

lattice::lattice(int length, int dim, double Bfield, int iterations, double Temp, int eq_time, string mode_for_sweep, string output_mode){
	
	L = length;
	d = dim;
	V = pow(L,d);
	
	avg_mag=0;
	avg_eng=0;
	
	t_eq = eq_time;
	
	T = Temp;
	b = 1/T;

	B = Bfield;
	iter = iterations;
	
	M = 1000;
	
	mode = mode_for_sweep;
	
	output = output_mode;
	
	cout<< "selected algorithm: "<<mode<<endl;
	
	spins.reserve(V);
	
	mag.reserve(iter);
	
	eng.reserve(iter);
	
	boot_samples.reserve(iter);
	
	boot_values.reserve(M);
	
	this->update_lookups();	
	
	// gsl function to initialize the rng
	rng = gsl_rng_alloc(gsl_rng_mt19937);
	
	// set seed based on current time
	long seed = time(NULL);	
	
	gsl_rng_set(rng,seed);	
	
	
	}

void lattice::update_lookups(){
	lookup_met_J.clear();
	lookup_met_B.clear();
	
	lookup_heat_J.clear();
	
	mag.clear();
	eng.clear();
	
	avg_mag = 0;
	avg_eng = 0;
	
	for(int q = -2*d; q <= 2*d; q += 2){
		lookup_met_J.push_back(exp(2*q*b));				// J = 1
	}
	
	for(int i = -1; i < 2; i += 2){					// i = +- 1
		lookup_met_B.push_back(exp(i*b*B));
	}
	
	
	for(int q = -2*d; q <= 2*d; q += 2){
		lookup_heat_J.push_back(1/(1+exp(-2*b*q)));				// J = 1
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


double lattice::get_eng(){
	double sum = 0;
	for(int i = 0; i < V; i++){
		sum	-= this->get_nn_sum(i)*spins[i];
		sum -= B*spins[i];
		}
	return sum/((double)V);
	}


void lattice::set_T(double Temp){
	T = Temp;
	b = 1/T;
	}

void lattice::sweep_met(){
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
			
			rho = lookup_met_J.at(((-s_j*q)/2)+d);					// Build up rho
			if(s_j < 0) rho *= lookup_met_B.at(0);
			else rho *= lookup_met_B.at(1);
			//cout<<"rho = "<<rho<<endl;
			
			rn = gsl_rng_uniform(rng);
			if(rn <= rho) spins.at(s_j_i) = -s_j;		// Accept flip with prob. rho
		}
	}
}

void lattice::sweep_heat(){
	int s_j_i,q;
	double rho,rn;
	for(int i = 0; i < V; i++){
		s_j_i = floor(V * gsl_rng_uniform(rng));
		q = this->get_nn_sum(s_j_i);
		rho = lookup_heat_J.at((q/2)+d);					// Build up rho
		rn = gsl_rng_uniform(rng);
		if(rho<=rn) spins.at(s_j_i)=-1;
		else spins.at(s_j_i)=1;
		}
	}

void lattice::run(){
	ostringstream fs;	
	ofstream file;
	
	this->update_lookups();
	
	this->equilibrate();
	
	fs<<"data/ising_"<<T<<".dat";
		
	file.open(fs.str().c_str());
	file<<"#t mag eng# at L="<<L<<" d="<<d<<endl;
	file.precision(10);
	cout<<"T = "<<T<<endl;
	file<<"#T = "<<T<<endl;
	
	int time_s = time(NULL);	
	
	double magn, engy;
	
	for(int t = 0; t<iter; t++){
		if(mode=="heatbath") this->sweep_heat();
		else this->sweep_met();
		
		magn = this->get_mag();

		engy = this->get_eng();
		
		if(output=="file")file<<t<<" "<<magn<<" "<<engy<<endl;	// Write measurements to file
		
		mag.push_back(magn);										// Save measurements in vectors
		
		eng.push_back(engy);	
		
		avg_mag += magn;
		avg_eng += engy;
		
		
		if((t%100) == 0){		// feedback every 100 sweeps
				int time_el = time(NULL)-time_s;
				int time_e = round((time_el)*(iter/double(t)-1));
			
				int min_e = floor(time_e/double(60));
				int min_el = floor(time_el/double(60));
				int s_e = time_e - min_e*60;
				int s_el = time_el - min_el*60;
				cout <<"\r";
				cout.flush();
				cout<<"Elapsed time: "<<min_el<<" min "<<s_el<<" sec; Estimated time: "<<min_e<<" min "<<s_e<<" sec  Progress: "<< round(100*t/iter)<<"%     ";
			}
		}
	
	avg_mag 	/= double(iter); 	
	avg_eng	/= double(iter);
	
	file.close();	
	
	}


void lattice::scan_t(){
	vector<double> t_vec;
	int num = 20;
	int t_c_num = 10;
	double start_t = 1;
	double end_t = 3;
	double margin = 0.25;
	double t_c = 2.3;
	double step = (end_t-start_t-2*margin)/double(num);
	double step_t_c = 2*margin/double(t_c_num);
	
	double cur_t= start_t;
	for (int i = 0; i <= num+t_c_num; i++){
		if(cur_t < t_c-margin || cur_t > t_c+margin){
			t_vec.push_back(cur_t);
			cur_t += step;
			}
		else{
			t_vec.push_back(cur_t);
			cur_t += step_t_c;
			}
		}
	cout<<"Scanning over the following range of temperatures:"<<endl;
	for (int i = 0; i < t_vec.size(); i++)cout<<t_vec.at(i)<<" ";
	cout<<endl;
	
	this->cold_start();
	ofstream file;
	file.open("data/t_scan.dat");
	file<<"# Lattice: "<<L<<"^"<<d<<endl;
	file<<"# T mag magerr sus suserr eng engerr heat heaterr"<<endl;
	file.close();
	for (int i = 0; i<t_vec.size(); i++){
		ofstream appfile;							// open(ios::app) files
		appfile.open("data/t_scan.dat",ios::app);
		appfile<<T<<" ";
		cout<<"T = "<<t_vec.at(i)<<endl;
		this->set_T(t_vec.at(i));
		this->update_lookups();
		this->run();
		appfile<<this->get_val("avg_mag")<<" ";
		this->calc_mag_corr();
		appfile<<this->get_std_err(this->get_vec("cov_mag"),this->get_vec("corr_mag"))<<" ";
		appfile<<this->get_mag_sus()<<" "<<this->get_mag_sus_err()<<" ";
		
		appfile<<this->get_val("avg_eng")<<" ";
		this->calc_eng_corr();
		appfile<<this->get_std_err(this->get_vec("cov_eng"),this->get_vec("corr_eng"))<<" ";
		appfile<<this-> get_spec_heat()<<" "<<this->get_spec_heat_err()<<endl;
		appfile.close();
		}
	}

double lattice::cov_func(int t_c, const vector<double>& y){		//compute covariance of a vector at time t_c
	double sum1 = 0, sum2 = 0, sum3 = 0;
	int size = int(y.size());
	for(int i = 0; i < (size - t_c); i++){
		sum1 +=  y.at(i)*y.at(i+t_c);
		sum2 +=  y.at(i);
		sum3 +=  y.at(i+t_c);	
		}
	double norm = 1/(double(size-t_c));
	
	return norm*(sum1 - norm * sum2 * sum3);
	}

void lattice::calc_cov_t(const vector<double>& vec, vector<double>& cov){
	for(int t_c = 0; t_c<int(vec.size())/10; t_c++){
		cov.push_back(cov_func(t_c,vec));
		}
	}

void lattice::calc_mag_cov(){
	calc_cov_t(mag,cov_mag);
	}

void lattice::calc_eng_cov(){
	calc_cov_t(eng,cov_eng);
	}

void lattice::calc_mag_corr(){
	this->calc_mag_cov();
	for(int i = 0; i<cov_mag.size(); i++){
		corr_mag.push_back(cov_mag.at(i)/cov_mag.at(0));
		}
	}

void lattice::calc_eng_corr(){
	this->calc_eng_cov();
	for(int i = 0; i<cov_eng.size(); i++){
		corr_eng.push_back(cov_eng.at(i)/cov_eng.at(0));
		}
	}

double lattice::calc_tau(const vector<double>& corr){
	double tau = 0.5;
	int flag = 0;
	int t = 1;
	while(flag == 0){
		tau += corr.at(t);
		if(corr.at(t)<0) flag = 1;
		t++;		
		}	
	return tau;
	}

double lattice::get_avg(const vector<double>& vec){
	double sum = 0;
	for(int i = 0; i<vec.size(); i++ ){
		sum+=vec.at(i);
		}
	return sum/double(vec.size());
	}

void lattice::equilibrate(){
	for(int i = 0; i<20*t_eq; i++){
		if(mode=="heatbath") this->sweep_heat();
		else this->sweep_met();		
		}
	}

double lattice::get_std_err(const vector<double>& cov,const vector<double>& corr){
	double variance = cov.at(0);
	int N = iter;
	double tau = this->calc_tau(corr);	
	return sqrt((variance/double(N))*2*tau);
	}

vector<double> lattice::get_vec(string choice){
	if(choice == "mag") return mag;
	else if(choice == "eng") return eng;
	else if(choice == "corr_mag") return corr_mag;
	else if(choice == "corr_eng") return corr_eng;
	else if(choice == "cov_mag") return cov_mag;
	else if(choice == "cov_eng") return cov_eng;
	
	vector<double> def_vec;
	return def_vec;
	}

double lattice::get_spec_heat(){
	double avg_eng_2 = this->get_pow_2_avg(eng);
	return (avg_eng_2 - pow(avg_eng,2))/pow(T,2);
	}

double lattice::get_spec_heat_err(){
	boot_values.clear();
	for(int i = 0; i<M; i++){
		boot_samples.clear();
		for(int j = 0; j<iter; j++){
			int rn_i = floor(iter*gsl_rng_uniform(rng));		// Random index between 0 and iter-1
			boot_samples.push_back(eng.at(rn_i));
			}
		double avg_boot_eng = this->get_avg(boot_samples);
		double avg_boot_eng_2 = this->get_pow_2_avg(boot_samples);
		boot_values.push_back((avg_boot_eng_2 - pow(avg_boot_eng,2))/pow(T,2));
		}
	double boot_avg = this->get_avg(boot_values);
	double sum = 0;
	for(int i = 0; i<boot_values.size(); i++){
		sum += pow(boot_values.at(i)-boot_avg,2);		
		}
	sum /= boot_values.size();
	return sqrt(sum);
	}

double lattice::get_mag_sus(){
	double avg_mag_2 = this->get_pow_2_avg(mag);
	return (avg_mag_2 - pow(avg_mag,2))/T;
	}

double lattice::get_mag_sus_err(){
	boot_values.clear();
	for(int i = 0; i<M; i++){
		boot_samples.clear();
		for(int j = 0; j<iter; j++){
			int rn_i = floor(iter*gsl_rng_uniform(rng));		// Random index between 0 and iter-1
			boot_samples.push_back(mag.at(rn_i));
			}
		double avg_boot_mag = this->get_avg(boot_samples);
		double avg_boot_mag_2 = this->get_pow_2_avg(boot_samples);
		boot_values.push_back((avg_boot_mag_2 - pow(avg_boot_mag,2))/T);
		}
	double boot_avg = this->get_avg(boot_values);
	double sum = 0;
	for(int i = 0; i<boot_values.size(); i++){
		sum += pow(boot_values.at(i)-boot_avg,2);		
		}
	sum /= boot_values.size();
	return sqrt(sum);
	}

double lattice::get_pow_2_avg(const vector<double>& vec){
	double sum = 0;
	for(int i = 0; i<vec.size(); i++){
		sum += pow(vec.at(i),2);		
		}	
	return sum/(double(vec.size()));
	}

double lattice::get_val(string choice){
	if(choice == "avg_mag") return avg_mag;
	else if(choice == "avg_eng") return avg_eng;
	
	double def;
	return def;
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

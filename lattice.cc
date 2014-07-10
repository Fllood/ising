#include <stdlib.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <ctime>

#include "lattice.h"

#include <gsl/gsl_rng.h>

#include "gnuplot_i.hpp" //Gnuplot class handles POSIX-Pipe-communication with Gnuplot

using namespace std;

lattice::lattice(int length, int dim, double Bfield, int iterations, double Temp, int eq_time, string mode_for_sweep, string output_mode, string start_mode){
	
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
	
	if(mode != "heatbath" && mode != "wolff") mode = "metro";
	
	output = output_mode;
	
	cout<< "selected algorithm: "<<mode<<endl;
	
	start = start_mode;
	
	spins.reserve(V);
	
	mag.reserve(iter);
	
	eng.reserve(iter);
	
	boot_samples.reserve(iter);
	
	cluster_sizes.reserve(iter);
	
	boot_values.reserve(M);
	
	
	s_ij_avg.reserve(V);	
	
	for(int i = 0; i < V; i++) {
		s_ij_avg.push_back(0);
		}
	
	this->update_lookups();	
	
	// gsl function to initialize the rng
	rng = gsl_rng_alloc(gsl_rng_mt19937);
	
	// set seed based on current time
	long seed = time(NULL);	
	
	gsl_rng_set(rng,seed);
	
	s_cl = floor(V * gsl_rng_uniform(rng));	
	
	
	
	}

void lattice::update_lookups(){
	lookup_met_J.clear();
	lookup_met_B.clear();
	
	lookup_heat_J.clear();
	
	mag.clear();
	eng.clear();
	
	cov_mag.clear();
	corr_mag.clear();
	
	cov_eng.clear();
	corr_eng.clear();
	
	cov_clu.clear();
	corr_clu.clear();
	
	boot_samples.clear();
	boot_values.clear();
	
	cluster_sizes.clear();
	
	corr_length_func.clear();
	r_values.clear();
	
	s_ij_avg.clear();

	for(int i = 0; i < V; i++) {
		s_ij_avg.push_back(0);
		}
			
	
	avg_mag = 0;
	avg_eng = 0;
	
	wolff_prob = 1 - exp(-2*b); 
	
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

void lattice::sweep_wolff(){
	int rs,s,nn;
	double rho;
	vector<int> added_spins, current_spins;
	
	rs = floor(V * gsl_rng_uniform(rng));			// choose random seed spin
	
	int old_value = spins.at(rs);
	int new_value = -old_value;
	
	spins.at(rs) = new_value;		// flip seed spin
	added_spins.push_back(rs);		
	
	bool flag = true;
	
	int cluster_size = 1;
	
	while(flag == true){
		
		flag = false;
		
		current_spins = added_spins;
		added_spins.clear();
		
		for(unsigned int j = 0; j < current_spins.size(); j++){
			
			s = current_spins.at(j);
			
			for(int i = d; i > 0; i--){
				int div = V/(pow(L,i));
				
				if((nn=s+div)>=V)nn -= V;			// helical boundary conditions
				if( spins[nn]==old_value){
					rho = gsl_rng_uniform(rng);
					if(rho <= wolff_prob){
						spins.at(nn) = new_value;	// flip spin
						added_spins.push_back(nn);
						flag = true;
						cluster_size++;
						}
					}
				
				if((nn=s-div)<0)nn += V;
				if(spins[nn]==old_value ){
					rho = gsl_rng_uniform(rng);
					if(rho <= wolff_prob){
						spins.at(nn) = new_value;	// flip spin
						added_spins.push_back(nn);
						flag = true;
						cluster_size++;
						}
					}
				
				}
			}
		}
	cluster_sizes.push_back(cluster_size);		// save cluster size
	}

void lattice::run(){
	
	this->update_lookups();
	
	this->equilibrate();
	
	ofstream file;
	string filename ("data/ising_temp_");
	ostringstream convert0;
	convert0<<T;
	filename.append(convert0.str());
	filename.append("_");

	filename.append(get_time_str());
	filename.append(".dat");
	file.open(filename.c_str());
	
	file<<"#t mag eng# at L="<<L<<" d="<<d<<endl;
	file.precision(10);
	file<<"#T = "<<T<<endl;
	
	int time_s = time(NULL);	
	
	double magn, engy;
	
	for(int t = 0; t<iter; t++){
		if(mode=="heatbath") this->sweep_heat();
		else if(mode =="wolff") this->sweep_wolff();
		else  this->sweep_met();
		
		// Calculation of averages
		this->calc_corr_length_avg();		
		
		magn = this->get_mag();

		engy = this->get_eng();
		
		if(output=="file" || output=="fileplot")file<<t<<" "<<magn<<" "<<engy<<endl;	// Write measurements to file
		
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
	int t_c_num = 20;
	double start_t = 1.5;
	double end_t = 3.5;
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
	for (unsigned int i = 0; i < t_vec.size(); i++)cout<<t_vec.at(i)<<" ";
	cout<<endl;
	
	if (start == "hot")this->hot_start();
	else this->cold_start();
	
	ofstream file;
	string filename ("data/t_scan_");
	filename.append(get_time_str());
	filename.append(".dat");
	if(!fexists(filename)){
		file.open(filename.c_str());
		}
	file<<"# Lattice: "<<L<<"^"<<d<<" B = "<<B<<" iterations = "<<iter<<" Algorithm: "<<mode<<endl;
	file<<"# T mag magerr sus suserr eng engerr heat heaterr";
	if(mode == "wolff"){
		file<<" clustersize clustersizeerr";		
		}
	file<<" engcorrtime magcorrtime";
	file<<endl;
	file.close();
	for (unsigned int i = 0; i<t_vec.size(); i++){
		ofstream appfile;							// open(ios::app) files
		appfile.open(filename.c_str(),ios::app);
		
		this->set_T(t_vec.at(i));
		this->update_lookups();
		
		appfile<<T<<" ";
		cout<<endl<<"T = "<<t_vec.at(i)<<endl;
		this->run();
		appfile<<this->get_val("avg_mag")<<" ";
		this->calc_mag_corr();
		appfile<<this->get_std_err(this->get_vec("cov_mag"),this->get_vec("corr_mag"))<<" ";
		appfile<<this->get_mag_sus()<<" "<<this->get_mag_sus_err()<<" ";
		
		appfile<<this->get_val("avg_eng")<<" ";
		this->calc_eng_corr();
		appfile<<this->get_std_err(this->get_vec("cov_eng"),this->get_vec("corr_eng"))<<" ";
		appfile<<this-> get_spec_heat()<<" "<<this->get_spec_heat_err();
		double avg_clu_size = 0;
		if(mode == "wolff"){
			avg_clu_size = this->get_avg(cluster_sizes);
			appfile<<" "<<avg_clu_size<<" ";
			this->calc_clu_corr();
			double err = this->get_std_err(this->get_vec("cov_clu"),this->get_vec("corr_clu"));
			appfile<<err;		
			}
		double tau_int = this->calc_tau(this->get_vec("corr_mag"));
		if(mode == "wolff")appfile<<" "<<tau_int*avg_clu_size/double(V)<<" ";
		else appfile<<" "<<tau_int<<" ";
		tau_int = this->calc_tau(this->get_vec("corr_eng"));
		if(mode == "wolff")appfile<<" "<<tau_int*avg_clu_size/double(V)<<" ";
		else appfile<<" "<<tau_int<<" ";
		appfile<<endl;
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

void lattice::calc_clu_cov(){
	calc_cov_t(cluster_sizes,cov_clu);
	}

void lattice::calc_mag_corr(){
	this->calc_mag_cov();
	for(unsigned int i = 0; i<cov_mag.size(); i++){
		corr_mag.push_back(cov_mag.at(i)/cov_mag.at(0));
		}
	}

void lattice::calc_eng_corr(){
	this->calc_eng_cov();
	for(unsigned int i = 0; i<cov_eng.size(); i++){
		corr_eng.push_back(cov_eng.at(i)/cov_eng.at(0));
		}
	}

void lattice::calc_clu_corr(){
	this->calc_clu_cov();
	for(unsigned int i = 0; i<cov_clu.size(); i++){
		corr_clu.push_back(cov_clu.at(i)/cov_clu.at(0));
		}
	}

double lattice::calc_tau(const vector<double>& corr){
	double tau = 0.5;
	int flag = 0;
	int t = 1;
	while(flag == 0 && t < int(corr.size()-1)){
		tau += corr.at(t);
		if(corr.at(t)<0) flag = 1;
		t++;		
		}	
	return tau;
	}

void lattice::calc_corr_length_avg(){
	int i = s_cl;
	//s_i_avg += spins[i]/double(iter);
	for(int j = 0; j < V; j++){
		s_ij_avg.at(j) += spins[i]*spins[j]/double(iter);
		//s_j_avg.at(j) += spins[j]/double(iter);		
		}	
	}

void lattice::calc_corr_length_func(){
	int i = s_cl;
	for(int j = 0; j<V; j++){
		double r = dist(i,j);
		if(find(r_values.begin(), r_values.end(), r)==r_values.end())r_values.push_back(r);
	}
	sort(r_values.begin(),r_values.end());
	for(unsigned int k = 0; k<r_values.size(); k++){
		int n_r = 0;
		double sum = 0;
		for(int j = 0; j<V; j++){
			double distance = dist(i,j);
			
			if(r_values.at(k) == distance ){
				n_r++;
				sum += s_ij_avg.at(j)-avg_mag*avg_mag;				
				}			
		}	
		sum /= double(n_r);
		
		corr_length_func.push_back(sum);
		
		}
	
	}

double lattice::get_avg(const vector<double>& vec){
	double sum = 0;
	for(unsigned int i = 0; i<vec.size(); i++ ){
		sum+=vec.at(i);
		}
	return sum/double(vec.size());
	}

double lattice::get_avg(const vector<int>& vec){
	double sum = 0;
	for(unsigned int i = 0; i<vec.size(); i++ ){
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
	else if(choice == "corr_clu") return corr_clu;
	else if(choice == "cov_mag") return cov_mag;
	else if(choice == "cov_eng") return cov_eng;
	else if(choice == "cov_clu") return cov_clu;
	
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
	for(unsigned int i = 0; i<boot_values.size(); i++){
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
	for(unsigned int i = 0; i<boot_values.size(); i++){
		sum += pow(boot_values.at(i)-boot_avg,2);		
		}
	sum /= boot_values.size();
	return sqrt(sum);
	}

double lattice::get_pow_2_avg(const vector<double>& vec){
	double sum = 0;
	for(unsigned int i = 0; i<vec.size(); i++){
		sum += pow(vec.at(i),2);		
		}	
	return sum/(double(vec.size()));
	}

double lattice::get_val(string choice){
	if(choice == "avg_mag") return avg_mag;
	else if(choice == "avg_eng") return avg_eng;
	
	double def = 0.0;
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

bool lattice::fexists(string filename){
	ifstream ifile(filename.c_str());
	return ifile;
	}

string lattice::get_time_str(){
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];
	
	time (&rawtime);
	timeinfo = localtime(&rawtime);
	
	strftime(buffer,80,"%d-%m-%Y_%I-%M-%S",timeinfo);
	string str(buffer);
	return str;
	}

bool lattice::in_vec(const vector<int>& vec, int a){
	bool found = false;
	for (unsigned int i = 0; i<vec.size(); i++){
		if(a == vec.at(i)) found = true;
		}
	return found;
	}

void lattice::one_temp(){
	if (start == "hot")this->hot_start();
	else this->cold_start();
	
	this->run();
	
	// fill vectors with G(r) and corresponding r values
	this->calc_corr_length_func();
	
	//write to file:
	ofstream file;
	string filename ("data/spatial_corr_");
	filename.append(get_time_str());
	filename.append(".dat");
	file.open(filename.c_str());
	file<<"# Lattice: "<<L<<"^"<<d<<" B = "<<B<<" iterations = "<<iter<<" Algorithm: "<<mode<<endl;
	file<<"# radius correlation"<<endl;
	file.precision(15);
	for(unsigned int i = 0; i < r_values.size(); i++ ){
		file<<r_values.at(i)<<" "<<corr_length_func.at(i)<<endl;			
		}
	file.close();
	
	// Magnetization measurement
	
	cout<<endl<<"avg mag: "<<avg_mag;
	
	int tau_eq_mag = 0;
	
	if(start == "hot"){
		while(fabs(mag.at(tau_eq_mag))<= fabs(avg_mag) || tau_eq_mag > iter) tau_eq_mag++;
		}
	else while(fabs(mag.at(tau_eq_mag))>= fabs(avg_mag) || tau_eq_mag > iter) tau_eq_mag++;	
	
	
	this->calc_mag_corr();	
	
	cout<<"(+/-)"<<this->get_std_err(cov_mag,corr_mag)<<endl;
	
	double mag_sus = 	this->get_mag_sus();
	double mag_sus_err = this->get_mag_sus_err();
	cout<<"mag suscep: "<<mag_sus<<"(+/-)"<<mag_sus_err<<endl;
	
	double tau_int_mag = this->calc_tau(corr_mag);
	
	cout<<"The integrated autocorrelation time of the magnetization is: "<<tau_int_mag<<endl;	
	
	// Energy measurement
	cout<<endl<<"avg eng: "<<avg_eng;	
	
	int tau_eq_eng = 0;
	
	if(start == "hot"){
		while(fabs(eng.at(tau_eq_eng))<= fabs(avg_eng) || tau_eq_eng > iter) tau_eq_eng++;
		}
	else while(fabs(eng.at(tau_eq_eng))>= fabs(avg_eng) || tau_eq_eng > iter) tau_eq_eng++;
	
	this->calc_eng_corr();	
	
	cout<<"(+/-)"<<this->get_std_err(cov_eng,corr_eng)<<endl;
	
	cout<<"specific heat: "<<this->get_spec_heat()<<"(+/-)"<<this->get_spec_heat_err()<<endl;
	
	double tau_int_eng = this->calc_tau(corr_eng);
	
	cout<<"The integrated autocorrelation time of the energy is: "<<tau_int_eng<<endl;
	
	if(mode == "wolff"){
		double avg_size = this->get_avg(cluster_sizes);
		this->calc_clu_corr();
		double err = this->get_std_err(this->get_vec("cov_clu"),this->get_vec("corr_clu"));
		double percent = 100*avg_size/double(V);
		cout<<endl<<"The average cluster size is: "<<avg_size<<"(+/-)"<<err<<" ("<<percent<<"% of lattice)"<<endl;		
		}	
	
	if(output == "plot" || output == "fileplot"){
	try
	{
		Gnuplot g1("lines");
		stringstream ss1;
		ostringstream convert1,convert2;
		convert1<<tau_int_mag;
		convert2<<tau_eq_mag;
		ss1<<"set termopt enhanced; set xlabel 'MC time ({/Symbol t}_{int} = "<<convert1.str()<<" {/Symbol t}_{eq} = "<<convert2.str()<<")'";
		g1<<ss1.str();
		g1<<"set ylabel 'Magnetization per spin'";
		g1<<"set xrange [0:3000]";
		
		g1.plot_x(mag,"Mag per spin versus MC time");
		
		/*
		string gnucmd = "";
		while(gnucmd != "exit" && gnucmd != "quit" && gnucmd != "stay"){
			cout<<"gnuplot>";
			std::getline (cin,gnucmd);
			g1.cmd(gnucmd);
		}
		*/
		
		g1<<"set term pdfcairo";
		stringstream ss2;
		ostringstream convert3;
		convert3<<T;
		ss2<<"set termopt enhanced; set output 'output/mag_plot_"<<mode<<"_"<<convert3.str()<<"_"<<get_time_str()<<".pdf'";
		g1<<ss2.str();
		g1<<"replot";
		g1<<"set term pop";
		
		
		Gnuplot g3("lines");
			
		g3.plot_x(corr_mag,"Correlation of mag per spin versus MC time");
		
		
		
		Gnuplot g2("lines");
		
		stringstream ss4;
		ostringstream convert4,convert5;
		convert4<<tau_int_eng;
		convert5<<tau_eq_eng;
		ss4<<"set termopt enhanced; set xlabel 'MC time ({/Symbol t}_{int} = "<<convert4.str()<<" {/Symbol t}_{eq} = "<<convert5.str()<<")'";
		g2<<ss4.str();
		g2<<"set ylabel 'Energy per spin'";
		g2<<"set xrange [0:3000]";
		
		g2.plot_x(eng,"Eng per spin versus MC time");
		
		/*
		gnucmd = "";
		while(gnucmd != "exit" && gnucmd != "quit" && gnucmd != "stay"){
			cout<<"gnuplot>";
			std::getline (cin,gnucmd);
			g2.cmd(gnucmd);
		}
		*/
		
		g2<<"set term pdfcairo";
		stringstream ss3;
		ostringstream convert6;
		convert6<<T;
		ss3<<"set termopt enhanced; set output 'output/eng_plot_"<<mode<<"_"<<convert6.str()<<"_"<<get_time_str()<<".pdf'";
		g2<<ss3.str();
		g2<<"replot";
		g2<<"set term pop";
		
		
		Gnuplot g4("lines");
			
		g4.plot_x(corr_eng,"Correlation of eng per spin versus MC time");
		
		
		this->wait_for_key();	
		}
	catch (GnuplotException ge){
        cout << ge.what() << endl;
		}
	}
	}
	
void lattice::wait_for_key(){
	std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
	#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)  // every keypress registered, also arrow keys
    cout << endl << "Press any key to continue..." << endl;

    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
    _getch();
	#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    cout << endl << "Press ENTER to continue..." << endl;

    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get();
	#endif
    return;
	}



double lattice::dist(int i, int j){
	int p;
	if(i >= V) i %= V;
	if(j >= V) j %= V;
	if(i<j) p = j-i;
	else p = i-j;
	
	int sum = 0;
	for(int dim = d; dim > 0; dim-- ){
		int fac = p / (pow(L,dim-1));
		int rest = p % int(pow(L,dim-1));
		p = rest;
		int taxlen = fac;
		if(fac >= (L-1)/sqrt(2)) taxlen = L - fac;
		sum += pow(taxlen,2);
		}
			
	return sqrt(sum);
	}
	

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <string>


using namespace std;


#include "lattice.h"	// lattice class

#include "gnuplot_i.hpp" //Gnuplot class handles POSIX-Pipe-communication with Gnuplot

void wait_for_key(); // Programm halts until keypress

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
	
	l1.hot_start();

	l1.run();
	
	
	// Magnetization measurement
	cout<<endl<<"avg mag: "<<l1.get_val("avg_mag");	
	
	
	l1.calc_mag_corr();	
	
	cout<<"(+/-)"<<l1.get_std_err(l1.get_vec("cov_mag"),l1.get_vec("corr_mag"))<<endl;
	
	double mag_sus = 	l1.get_mag_sus();
	double mag_sus_err = l1.get_mag_sus_err();
	cout<<"mag suscep: "<<mag_sus<<"(+/-)"<<mag_sus_err<<endl;
	
	double tau_int = l1.calc_tau(l1.get_vec("corr_mag"));
	
	cout<<"The integrated autocorrelation time of the magnetization is: "<<tau_int<<endl;	
	
	// Energy measurement
	cout<<endl<<"avg eng: "<<l1.get_val("avg_eng");	
	
	l1.calc_eng_corr();	
	
	cout<<"(+/-)"<<l1.get_std_err(l1.get_vec("cov_eng"),l1.get_vec("corr_eng"))<<endl;
	
	cout<<"std dev eng: "<<sqrt(l1.get_vec("cov_eng").at(0))<<endl;
	
	cout<<"specific heat: "<<l1.get_spec_heat()<<"(+/-)"<<l1.get_spec_heat_err()<<endl;
	
	tau_int = l1.calc_tau(l1.get_vec("corr_eng"));
	
	cout<<"The integrated autocorrelation time of the energy is: "<<tau_int<<endl;
	
	try
	{
		Gnuplot g1("lines");
		
		g1.plot_x(l1.get_vec("mag"),"Mag per spin versus MC time");
		
		
		
		Gnuplot g3("lines");
			
		g3.plot_x(l1.get_vec("corr_mag"),"Correlation of mag per spin versus MC time");
		

		
		Gnuplot g2("lines");
		
		g2.plot_x(l1.get_vec("eng"),"Eng per spin versus MC time");
		
		
		
		Gnuplot g4("lines");
			
		g4.plot_x(l1.get_vec("corr_eng"),"Correlation of eng per spin versus MC time");
		
		wait_for_key();	
		}
	catch (GnuplotException ge)
    {
        cout << ge.what() << endl;
    }
	
	
	return 0;
}

void wait_for_key ()
{
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


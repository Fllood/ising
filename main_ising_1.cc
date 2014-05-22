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
	
	l1.calc_mag_corr();
	
	cout<<"std dev: "<<sqrt(l1.get_vec("corr_mag").at(0))<<endl;
	
	try
	{
		Gnuplot g1("lines");
		
		g1.plot_x(l1.get_vec("mag"),"Mag per spin versus MC time");
		
		wait_for_key();
			
		g1.plot_x(l1.get_vec("corr_mag"),"Correlation of mag per spin versus MC time");
		
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


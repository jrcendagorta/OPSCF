#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <string>
#include <cstring>
#include <sstream>
#include "System.h"
#include "Initial.h"
#include "Integrators.h"
#include "Force.h"

#define TWO_PI 6.2831853071795864769252866
#define PI  3.14159265358979323846

using namespace std;


/**
* Currently the program requires a few command line arguments. Further implementations
* will read in an input file with the arguments of the simulation. 
* Specifically, the code requires 6 arguments. 
* 1) # path integral beads
* 2) Beta = 1/kbT
* 3) Time for the correlation function
* 4) Number of simulation steps
* 5) timestep of simulation  (dt)
* 6) Potential type (1=HO, 2=Mildly Anharmonic, 3= Quartic, 4= DW)
**/

int main(int argc, char* argv[]) {

    System s;
    srand (100);


    int P = atoi(argv[1]);
    double beta = atof(argv[2]);
    double time = atof(argv[3]);  
    int nsteps = atoi(argv[4]);
    double dt = atof(argv[5]);
    int pot_typ = atoi(argv[6]);
    
    int nwrite = 100; // Write out to files every 100 steps, this can be altered.    
    
    double mass = 1.0;
    
    // Allocate memory 
    Allocate(&s,P);
    double r;
       
    // Initialize variables and perform initial force calculation.
    Initialize(&s,beta,mass,dt,time,pot_typ);
	
    	
	
	// Output Files for averages and instantaneous quantities
	string file1,file2,file3;
	stringstream ss1, ss2, ss3;
	ss1 << P;
	ss2 << beta;
	ss3 << time;
	string sss1 = ss1.str();
	string sss2 = ss2.str();
	string sss3 = ss3.str();
	file1 = "average_P" + sss1 + "_Beta-"+ sss2 + "_T-"+sss3 + ".out";
	file2 = "instant_P" + sss1 + "_Beta-"+ sss2 + "_T-"+sss3 + ".out";
	char * fd1 = new char[file1.length() + 1];
	char * fd2 = new char[file2.length() + 1];
	strcpy(fd1,file1.c_str());
	strcpy(fd2,file2.c_str());
    ofstream aver,inst;
    aver.open(fd1);
    inst.open(fd2);
        

    double Gt = 0.0;
    double Gt_avg = 0.0;
    double Gt_ste = 0.0;

    double Gt_2 = 0.0;
    double Gt_avg_2 = 0.0;
    double Gt_ste_2 = 0.0;
    
    double temp_avg = 0.0;
    double temp_avg_2 = 0.0;
    
    int count = 1;
    for(int m = 0; m < nsteps; m++) {            
        Integrate_Nose(&s, dt);
        Integrate_Velocity(&s, dt);
        
        Integrate_Position(&s,dt/2.0);
        Integrate_OU(&s, dt);
        Integrate_Position(&s,dt/2.0);

        TransformP_TtoP(&s); 
        Force_Control(&s,pot_typ);
        Integrate_Velocity(&s,dt);

        Integrate_Nose(&s, dt);

        Gt = s.position_p[0].x*s.position_p[P].x;        
        Gt_2 = Gt*Gt;

        if((m+1)%nwrite == 0) {
            inst << fixed <<right<< setprecision(3) << setw(15) << m 
                << fixed << setprecision(10) << setw(25) << Gt 
                << fixed << setprecision(10) << setw(25) << Gt_2 << endl;
        }
        if(count == 1) {
            Gt_avg = Gt;
            Gt_ste = 0.0;
            Gt_avg_2 = Gt_2;
            Gt_ste_2 = 0.0;

            }
        else {
            temp_avg = Gt_avg;
            temp_avg_2 = Gt_avg_2;
            Gt_avg = Gt_avg + (Gt-Gt_avg)/((double)count);
            Gt_ste = Gt_ste + (Gt-temp_avg)*(Gt-Gt_avg);
            Gt_avg_2 = Gt_avg_2 + (Gt_2-Gt_avg_2)/((double)count);
            Gt_ste_2 = Gt_ste_2 + (Gt_2-temp_avg_2)*(Gt_2-Gt_avg_2);
        }
            
        count++;
        if((m+1)%nwrite == 0) {
            aver << fixed <<right<< setprecision(3) << setw(15) << (m+1)*dt
                << fixed << setprecision(10) << setw(20) << Gt_avg
                << fixed << setprecision(10) << setw(20) << sqrt(Gt_ste/((double)count)/((double)count-1))
                << fixed << setprecision(10) << setw(20) << Gt_avg_2 
                << fixed << setprecision(10) << setw(20) << sqrt(Gt_ste_2/((double)count)/((double)count-1)) << endl;
         }
       
    }
    Deallocate(&s);
    aver.close();
    inst.close();
    delete [] fd1;
    delete [] fd2;

	return 0; 

}

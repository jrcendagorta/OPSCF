#include "System.h"
#include "Force.h"
#include "Initial.h"
#include <cstdlib>
#include <cmath>
#include <iostream>

#define TWO_PI 6.2831853071795864769252866
#define PI  3.14159265358979323846

double randGaus(double sigma) {
        static bool haveSpare = false;
        static double rand1, rand2;
        rand2 = rand()/ ((double) RAND_MAX);
        rand2 = 1.0-rand2;
        if(rand2 < 1e-100) rand2 = 1e-100;
        rand2 = -2*log(rand2);
        rand1 = (rand() / ((double) RAND_MAX)) * TWO_PI;
        return sigma*sqrt(rand2)*cos(rand1);
}

void Allocate(System *s, int& P) {
    s->P = P;
    s->position_p = new pnt [P+1];
    s->velocity_p = new pnt [P+1];
    s->force_p    = new pnt [P+1];

    s->mass_harm  = new double [P+1];
    s->lambda     = new double [P+1];

    s->position_t = new pnt [P+1];
    s->velocity_t = new pnt [P+1];
    s->force_t    = new pnt [P+1];
    
    s->force_harm = new pnt [P+1];
    s->mass_vel   = new double [P+1];
    
    s->force_matrix_t = new pnt [P+1];
    s->force_matrix_p = new pnt [P+1];
    
    s->forces_total = new pnt [P+1];
    
    s->m.M_di = new double [P-1];
    s->m.M_di_inv = new double [P-1];
    int num_off = (P*P - P)/2;
    s->m.M_off_inv = new double [num_off];

    s->m.K = new double [P-1];
    s->m.KM = new double [P-1];
    s->U1=new double [P+1];
    s->U2=new double [P+1];
    s->U3=new double [P+1];
    

    
    for(int i = 0; i <= P; i++) {
        s->position_p[i] = 0.0;
        s->velocity_p[i] = 0.0;
        s->force_p[i] = 0.0;
        s->mass_harm[i] = 0.0;
        s->lambda[i] = 0.0;
        s->force_harm[i] = 0.0;
	    s->position_t[i] = 0.0;
	    s->velocity_t[i] = 0.0;
	    s->force_t[i] = 0.0;
	    s->force_matrix_t[i] = 0.0;	
	    s->force_matrix_p[i] = 0.0;
	    s->forces_total[i] = 0.0;	
	    s->mass_vel[i] = 0.0; 
	    s->U1[i] = 0.0;
        s->U2[i] = 0.0;
        s->U3[i] = 0.0;   
	}
	for(int i = 0; i < P-1; i++) {
	    s->m.K[i] = 0.0;
	    s->m.KM[i] = 0.0;
        s->m.M_di[i] = 0.0;
        s->m.M_di_inv[i] = 0.0; 

	}
	for(int i = 0; i < num_off; i++) {
	    s->m.M_off_inv[i] = 0.0;
	}
}

void Deallocate(System *s) {
    delete [] s->position_p;
    delete [] s->velocity_p;
    delete [] s->force_p;    
    delete [] s->m.M_di;
    delete [] s->m.M_di_inv;
    delete [] s->m.M_off_inv;

    delete [] s->m.K;
    delete [] s->m.KM;
    
    delete [] s->lambda; 

    delete [] s->position_t;
    delete [] s->velocity_t;
    delete [] s->force_t;
    
    delete [] s->force_harm;
    delete [] s->force_matrix_p;
    delete [] s->force_matrix_t;
    
    delete [] s->forces_total;
   
    delete [] s->mass_vel;
    delete [] s->mass_harm; 
   

    delete [] s->wdt;
    delete [] s->U1;
    delete [] s->U2;
    delete [] s->U3;

    
    delete [] s->isok.v1;
    delete [] s->isok.v2;
    delete [] s->isok.Q1;
    delete [] s->isok.Q2;
    delete [] s->isok.fric;


}

void Initialize_Mass(System *s, double mass) {
    
    int P = s->P;
    for(int i = 0; i <= P; i++) {
        s->mass_harm[i] = mass;   
	    if(i == 0) {
    	    s->lambda[i] = 0;
		    s->mass_vel[i] = mass;
	    }
	    else if(i == P) {
	        s->lambda[i] = 1.0/((double)P);;
	        s->mass_vel[i] = mass*s->lambda[i];
	    }
	    else {
	        s->lambda[i] = ((double) i+1)/((double)i);
		    s->mass_vel[i] = mass*s->lambda[i];
	    }
	    

    } 
}

void Initialize_Position(System *s) {
    int P = s->P;
    double sigma = 0.0;
    
    // Initilization of Positions
    // This need to be better developed
    // Currently the initial bead is set to zero and all subsequent beads
    // are displaced randomly from the previous bead
	for(int i = 0; i <= P; i++) {
		if(i == 0) {
			s->position_p[i].x = 0.0;
		}
		else {
			sigma = 0.01;
			s->position_p[i].x = s->position_p[i-1].x+randGaus(sigma);
		
		}
		
	}
	
	// Transformation into transformed staging coordinates
    s->position_t[0].x = 0.5*(s->position_p[0].x + s->position_p[P].x);
    
    s->position_t[P].x = (s->position_p[0].x - s->position_p[P].x);
    
    for(int i = P-1; i > 0; i--) {
        s->position_t[i].x = s->position_p[i].x - (((double)i)*s->position_p[i+1].x + s->position_p[0].x)/((double) i+1); 
    }

  

}

void Initialize_Velocity(System *s) {
    double sigma = 0.0;
    int P = s->P;
    
    // Velocity are randomly samples from Gaussian distribution
    // THere is no need to transform to primitive vleocities. 
    for(int i = 0; i <= P; i++) {
	    sigma = 0.0;
	    sigma = sqrt(s->mass_vel[i]/s->beta);
	    s->velocity_t[i].x = randGaus(sigma);
    }  

}

void Initialize_Isokinetic(System* s, int L) {
    int P = s->P;

    s->isok.v1 = new pnt [(P+1)*L];
    s->isok.v2 = new pnt [(P+1)*L];
    s->isok.fric = new double [P+1];
    s->isok.Q1 = new double [P+1];
    s->isok.Q2 = new double [P+1];
    
// Initialization of these values can be optimized and have not been explored in depth  
    for(int i = 0; i <= P; i++) {
        if( i == 0 || i == P) {
            s->isok.fric[i] = s->wp;
        }
        else {
            s->isok.fric[i] = s->wp;
        }
        s->isok.Q1[i] = 1.0/s->beta;
        s->isok.Q2[i] = 1.0/(s->wp*s->wp*s->beta);
 
        for(int j = 0; j < s->L; j++) {
            s->isok.v2[j*(P+1)+i].x = 0.1;
            s->isok.v1[j*(P+1)+i].x = sqrt(2.0*(1.0/s->beta)/s->isok.Q1[i]);
        }
    }
    
}


void Initialize(System *s, double beta, double mass, double timestep, double time, int pot_typ) {
    int P = s->P;
    s->hbar = 1.0;
    s->beta = beta;
    
    s->dt = timestep;
    s->L = 1;
    
    s->t = time;
    s->tau_c = sqrt(beta*beta*s->hbar*s->hbar/4.0 + s->t*s->t);
    s->gamma = mass*((double)P)*s->t/s->hbar/s->tau_c/s->tau_c;
    s->alpha = mass*((double)P)*beta/4.0/s->tau_c/s->tau_c;
    s->wp = sqrt((double)P)/(s->tau_c);

    Initialize_Mass(s, mass);    
    Initialize_Position(s);
    Initialize_Velocity(s);

    
    Initialize_Isokinetic(s,s->L);
    s->wdt = new double [3];
    s->nsy = 3;
    s->nres = 5;
   
    double pp = 1.0/(2.0 - pow(2.0,(1.0/3.0)));
    s->wdt[0] = pp*s->dt/((double) s->nres);
    s->wdt[1] = -pow(2.0,(1.0/3.0))*pp*s->dt/((double) s->nres);
    s->wdt[2] = pp*s->dt/((double) s->nres);
    
    
    Force_Control(s,pot_typ);
    

}

/**
* Transformation of Position Coordinates from Staging to Primitive that 
* can be used for force evauation. 
**/

void TransformP_TtoP(System *s) {
    int P = s->P;
    double xsum = 0.0;
    s->position_p[0].x = s->position_t[0].x + 0.5*s->position_t[P].x;

    s->position_p[P].x = s->position_t[0].x - 0.5*s->position_t[P].x;
    
    for(int i = P-1; i > 0; i--) {
        xsum = 0.0;
        s->position_p[i].x = s->position_t[0].x + (((double) P)/2.0-((double)i))/((double)P)*s->position_t[P].x;
        for(int j = i; j < P; j++) {
            s->position_p[i].x += ((double)i)/((double)j)*s->position_t[j].x;
        }  
    } 

}


/** 
* Transformation of Force Cooridinates from Primitive to Staging Coordinates
*
**/
void TransformF_PtoT(System *s) {
    int P = s->P;
    double fsum = 0.0;
    for(int j = 0; j <= P; j++) {
        fsum = 0.0;
        if(j == 0) {   
            for(int i = 0; i <= P; i++) {
                fsum += s->force_p[i].x;
            }
            s->force_t[j].x = fsum;
        }
        else if(j == P) {
            s->force_t[j].x = 0.5*(s->force_p[0].x - s->force_p[P].x);
            for(int i = 1; i < P; i++) {
                s->force_t[j].x += (((double)P)/2.0 - ((double) i))/((double)P)*s->force_p[i].x;
            }
        }
        else {
            s->force_t[j].x = s->force_p[j].x + ((double)j-1)/((double)j)*s->force_t[j-1].x;
        }
    }    
}





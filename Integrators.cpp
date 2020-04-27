#include "System.h"
#include "Initial.h"
#include "Integrators.h"
#include <cstdlib>
#include <cmath>
#include <iostream>



void Integrate_Position(System* s, double dt) {
    int P = s->P;
    for(int i = 0; i <= P; i++) {
        s->position_t[i].x += dt*s->velocity_t[i].x;      
    }
    
}

void Integrate_Velocity(System* s, double dt) {
    int P = s->P;
    pnt* forces;
    forces = s->forces_total; //allocate new force
    
    for(int i = 0; i <= P; i++) {
        forces[i].x = 0.0;        
        forces[i].x += s->force_matrix_t[i].x;
        forces[i].x += s->force_t[i].x;
        forces[i].x += s->force_harm[i].x;
    }
    
    double arg = 0.0;
    double b = 0.0;
    double st = 0.0;
    double rb = 0.0;
    double sdot = 0.0;
    double a = 0.0;
    int j_ind;
    double dt2 = dt/2.0;
    
    for(int i = 0; i <= P; i++) {
        a = forces[i].x*s->velocity_t[i].x*s->beta/((double) s->L);
        b = forces[i].x*forces[i].x/s->mass_vel[i]*s->beta/((double) s->L);

        rb = sqrt(b);
        arg = dt*rb*0.5;
        
        if(arg > 0.00001) {            
            st = (1.0/rb)*sinh(arg) + (a/b)*(cosh(arg)-1.0);
            sdot = cosh(arg) + (a/rb)*sinh(arg);
        }
        else {            
            st = ((((b*a/24.0)*dt2 + b/6.0)*dt2 + 0.5*a)*dt2 + 1.0)*dt2;
            sdot = (((b*a/6.0)*dt2 + 0.5*b)*dt2 + a)*dt2 + 1.0;           
        }
        
        s->velocity_t[i].x = (s->velocity_t[i].x + forces[i].x/s->mass_vel[i]*st)/sdot;
        
        for(int j = 0; j < s->L; j++) {
            j_ind = i+j*(P+1);
            s->isok.v1[j_ind].x /= sdot;
        }
    
    }
}

void Integrate_Nose(System* s, double dt) {
    int P = s->P;
    double G = 0.0;
    double H = 0.0;
    double H_sum = 0.0;
    double ev2 = 0.0;
    int j_ind;
    double dti = 0.0;

    for(int i = 0; i <= P; i++) {
        for(int idt = 0; idt < s->nsy; idt++){
            dti = s->wdt[idt];
            
	        for(int ires = 0; ires < s->nres; ires++){
                H_sum = 0.0;
                for(int j = 0; j < s->L; j++) {
                    j_ind = i+j*(P+1);
                    G = (s->isok.Q1[i]*s->isok.v1[j_ind].x*s->isok.v1[j_ind].x - 1.0/s->beta)/s->isok.Q2[i];
                    s->isok.v2[j_ind].x += 0.25*dti*G;
                    ev2 = exp(-dti*0.5*s->isok.v2[j_ind].x);
                    H_sum += s->isok.Q1[i]*s->isok.v1[j_ind].x*s->isok.v1[j_ind].x*ev2*ev2;
                }
                H_sum = H_sum*((double) s->L)/((double) s->L+1);
                H = sqrt(s->L/s->beta/( s->mass_vel[i]*s->velocity_t[i].x*s->velocity_t[i].x + H_sum));
      
                s->velocity_t[i].x *= H;

                for(int j = 0; j < s->L; j++) {
                    j_ind = i+j*(P+1);
                    s->isok.v1[j_ind].x *= (H*exp(-dti*0.5*s->isok.v2[j_ind].x));
                }
                for(int j = 0; j < s->L; j++) {
                    j_ind = i+j*(P+1);
                    G = (s->isok.Q1[i]*s->isok.v1[j_ind].x*s->isok.v1[j_ind].x - 1/s->beta)/s->isok.Q2[i];
                    s->isok.v2[j_ind].x += 0.25*dti*G;
                }
            }
        } 
    }
} 

void Integrate_OU(System*s, double dt) {
    int P = s->P;
    for(int i = 0; i <= P; i++) {
        s->isok.v2[i].x = exp(-s->isok.fric[i]*dt)*s->isok.v2[i].x + sqrt((1.0-exp(-s->isok.fric[i]*dt*2.0))/(s->beta*s->isok.Q2[i]))*randGaus();
    }

}  




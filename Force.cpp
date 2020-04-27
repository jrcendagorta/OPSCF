#include "System.h"
#include "Force.h"
#include "Initial.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>


/**
* This is a mtrix inversion scheme specific to a symmetric trdiagonal matrix.
* This is not a general matrix inversion scheme. This is the most important part 
* as it solves for the inversion of the matrix. Specifically, the function parameters are
* double * d : diagonal elements of the original matrix M
* int n      : dimension of the matrix (P-1)
* double * md: diagonal elements of the inverted matrix (return value)
* double * mo: off-diagonal elements of inverted matrix (return value)
* double a   : parameter alpha that correpsond to the off-diagonal element of the original matrix M
* double det : determinant of the matrix M (return value)
*
**/

double InverseMatrix(double* d, int n, double* md, double *mo, double a, double* det) {
    double * theta;
    double * psi; 
    theta = new double [n+1];
    psi = new double [n+1];
    
    for(int i = 0; i <= n; i++) {
        theta[i] = 0.0;
        psi[i] = 0.0;
    }
    
    theta[0] = 1.0;
    theta[1] = d[0];
    psi[n-1] = d[n-1];
    psi[n] = 1.0;
    for(int i = 2; i <= n; i++) {
        theta[i] = d[i-1]*theta[i-1]-theta[i-2];
    }
    for(int i = n-2; i >= 0; i--) {
        psi[i] = d[i]*psi[i+1] - psi[i+2];
    }
    int count = 0;
    for(int i = 0; i < n; i++) {
        for(int j = i; j < n; j++) {
            if(i == j) {
                md[i] = theta[i]*psi[j+1]/theta[n]/a;
            }
            else {
                mo[count] = theta[i]*psi[j+1]/theta[n]/a; 
                count++; 
            }
        }
    }
    det[0]=theta[n]*pow(a,n);

    delete [] theta;
    delete [] psi;
}


void Force_Control(System *s, int pot_typ) {
    Derivatives(s,pot_typ);
    Matrix_Calc(s);
    Harm_Calc(s);
    Force_Calc(s);    
}

void Derivatives(System *s, int pot_typ) {
    int P = s->P;
    for(int i = 0; i <= P; i++) {
        if(pot_typ == 1) {// Harmonic Oscillator           
            s->U1[i] = 1.0*s->position_p[i].x;
            s->U2[i] = 1.0;
            s->U3[i] = 0.0; 
        }
        else if(pot_typ == 2) {// Anaharmonic Oscillator
            s->U1[i] = s->position_p[i].x + 0.3*s->position_p[i].x*s->position_p[i].x + 0.04*s->position_p[i].x*s->position_p[i].x*s->position_p[i].x;
            s->U2[i] = 1.0 + 0.6*s->position_p[i].x+ 0.12*s->position_p[i].x*s->position_p[i].x;
            s->U3[i] = 0.6 + 0.24*s->position_p[i].x;
        }
        else if(pot_typ == 3) {// Quartic Oscillator
            s->U1[i] = s->position_p[i].x*s->position_p[i].x*s->position_p[i].x;
            s->U2[i] = 3.0*s->position_p[i].x*s->position_p[i].x;
            s->U3[i] = 6.0*s->position_p[i].x;
        }
        else if(pot_typ == 4) {
            double a = 2.5;
            double D = 0.1;
            s->U1[i] = 4.0*D*s->position_p[i].x*(s->position_p[i].x*s->position_p[i].x - a);
            s->U2[i] = 4.0*D*(3.0*s->position_p[i].x*s->position_p[i].x - a);
            s->U3[i] = 24.0*D*s->position_p[i].x;
        }
    }

}


// Calculate Primitive force from Potential
void Force_Calc(System *s) {
    int P = s->P ;
    for(int i = 0; i <= P; i++) {
       if(i == 0 || i == P) {
            s->force_p[i].x = -s->U1[i]/((double) 2*P);
        }
        else {
            s->force_p[i].x = -s->U1[i]/((double) P);
        }
    
    }
    TransformF_PtoT(s);

}

// Calculate Harmonic Forces
void Harm_Calc(System *s) {
    int P = s->P;
    for(int i = 0; i <= P; i++) {
        s->force_harm[i].x = -s->mass_harm[i]*s->lambda[i]*s->wp*s->wp*s->position_t[i].x;
    }
    
}


// This is the Matrix Force Calculation, (i.e. the heart of OPSCF method)
void Matrix_Calc(System *s) {
    int P = s->P;
    double t = s->t;
    double hbar = s->hbar;
    double gamma = s->gamma;
    double alpha = s->alpha;
    double *K = s->m.K;
    double * KM = s->m.KM;
    double *M_di = s->m.M_di;
    double *M_di_inv = s->m.M_di_inv;
    double *M_off_inv = s->m.M_off_inv;
    pnt* fm_p = s->force_matrix_p;
    pnt* fm_t = s->force_matrix_t;    
    pnt* p = s->position_p;   
    double det, KMK;
    double K1, K2, dln;

    
    
   
    K1 = 0.0; K2 = 0.0; dln = 0.0;
    
    for(int i = 0; i < P-1; i++) {
        KM[i] = 0.0;
    }


//  Define the matrix diagonal elements, off diagonal elements are constant
//  Also define the vector K elements 
//  Matrix diagonal elements and K-elements are defined in the original JCP paper.
    for(int i = 0; i < P-1; i++) {     
        K[i] =gamma*(2*p[i+1].x - p[i].x - p[i+2].x)-t/((double) P)/hbar*s->U1[i+1];
        M_di[i] = 2.0 + s->beta/4.0/((double)P)*s->U2[i+1]/alpha;   
    }

    InverseMatrix(M_di, P-1, M_di_inv, M_off_inv,alpha,&det);
  
  
// Evaluate the K^T*M^{-1} terms 
// Matrix multiplication 
    int num_off = (P*P - P)/2;
    int count = 0;
    
    for(int i = 0; i < P-1; i++) {
        KM[i] = K[i]*M_di_inv[i];
        for(int j = i+1; j < P-1; j++) {
            KM[i] +=K[j]*M_off_inv[count];
            count++;
        } 
    } 
    
    count = 0;
    
    for(int i = 0; i < P-2; i++) {
        for(int j = i+1; j < P-1; j++) {
            KM[j] += K[i]*M_off_inv[count];
            count++;
        }
        
    }
    

    KMK = 0.0;
    for(int i = 0; i < P-1; i++ ){
                 KMK += K[i]*KM[i];
             
    }
    

// Evaluate primitive forces  
    for(int m = 0; m <= P; m++) {
       K2 = 0.0;
       K1 = 0.0;
       dln = 0.0;
        if(m == 0) {
            K1 = -(KM[0])*2.0*gamma;
            
        }
        else if(m == P) {
            K1 = -(KM[P-2])*2.0*gamma;
        }
        else {

            K2 = (KM[m-1])*(KM[m-1])*s->beta/4.0/((double)P)*s->U3[m];
            dln = M_di_inv[m-1]*s->beta/4.0/((double)P)*s->U3[m];
            if(m == 1) {
                K1 = (2.0*gamma - t/((double)P)/hbar*s->U2[m])*(KM[0])*2.0;
                K1 += -gamma*2.0*(KM[1]);
                
                  
            } 
            else if(m == P-1) {
                K1 = (2.0*gamma - t/((double)P)/hbar*s->U2[m])*(KM[P-2])*2.0;
                K1 += -gamma*2.0*(KM[P-3]);                
            }
            else {
                K1 = (2.0*gamma - t/((double)P)/hbar*s->U2[m])*(KM[m-1])*2.0;
                K1 += -gamma*2.0*((KM[m])+(KM[m-2]));               
            }
        }    
        
        fm_p[m].x = -(K1 - K2 + dln)/(2.0*s->beta);
    }

// Transform primitive matrix forces to transformed staging forces   
    double fsum = 0.0;
    for(int j = 0; j <= P; j++) {
        fsum = 0.0;
        if(j == 0) {   
            for(int i = 0; i <= P; i++) {
                fsum += fm_p[i].x;
            }
            fm_t[j].x = fsum;
        }
        else if(j == P) {
            fm_t[j].x = 0.5*(fm_p[0].x - fm_p[P].x);
            for(int i = 1; i < P; i++) {
                fm_t[j].x += (((double)P)/2.0 - ((double) i))/((double)P)*fm_p[i].x;
            }
        }
        else {
            fm_t[j].x = fm_p[j].x + ((double)j-1)/((double)j)*fm_t[j-1].x;
        }
        
        
    }   
   

    
    
}

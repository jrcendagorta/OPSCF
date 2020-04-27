#ifndef SYSTEM_H
#define SYSTEM_H

#include <cmath>

/**
* This is the system header file. Even though these simulations are performed in 
* one dimension, the entire architecture is based on the a point structure called 
* pnt. The point structure has some basic operations.

**/



struct pnt {
    double x;
    double y;
    double z;
    
    pnt& operator=(const double& scalar) {
            x = scalar;
            y = scalar;
            z = scalar;
            return* this;
    } 
    
    pnt& operator=(const pnt& rhs) {
            x = rhs.x;
            y = rhs.y;
            z = rhs.z;
            return *this;
    }
    
    pnt& operator*=(const pnt& rhs) {
            x *= rhs.x;
            y *= rhs.y;
            z *= rhs.z;
            return* this;
    } 
    
    pnt& operator*=(const double& scalar) {
            x *= scalar;
            y *= scalar;
            z *= scalar;
            return *this;
    }
    
    friend pnt operator*(pnt lhs, const double& scalar) {
            lhs.x *= scalar;
            lhs.y *= scalar;
            lhs.z *= scalar;
            return lhs;
    } 
    
    friend pnt operator*( const double& scalar,pnt rhs) {
            rhs.x *= scalar;
            rhs.y *= scalar;
            rhs.z *= scalar;
            return rhs;
    } 
    
    friend pnt operator*( const pnt& lhs,pnt rhs) {
            rhs.x *= lhs.x;
            rhs.y *= lhs.y;
            rhs.z *= lhs.z;
            return rhs;
    } 
    
    
    pnt& operator/=(const pnt& rhs) {
            x /= rhs.x;
            y /= rhs.y;
            z /= rhs.z;
            return* this;
    } 
    
    pnt& operator/=(const double& scalar) {
            x /= scalar;
            y /= scalar;
            z /= scalar;
            return *this;
    }
    
    friend pnt operator/(pnt lhs, const double& scalar) {
            lhs.x /= scalar;
            lhs.y /= scalar;
            lhs.z /= scalar;
            return lhs;
    } 
    
    friend pnt operator/(const double& scalar, pnt rhs) {
            rhs.x /= scalar;
            rhs.y /= scalar;
            rhs.z /= scalar;
            return rhs;
    } 
    
    friend pnt operator-(pnt lhs, const double& scalar) {
            lhs.x -= scalar;
            lhs.y -= scalar;
            lhs.z -= scalar;
            return lhs;
    } 
    
    friend pnt operator-(const double& scalar, pnt rhs) {
            rhs.x -= scalar;
            rhs.y -= scalar;
            rhs.z -= scalar;
            return rhs;
    } 
    
    
    pnt& operator+=(const double& scalar) {
            x += scalar;
            y += scalar;
            z += scalar;
            return* this;
    } 
    
    pnt& operator+=(const pnt& rhs) {
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;
            return* this;
    } 

    
    
};


/**
* This is the structure to define the isokinetic thermostat system for the system.
*
**/

struct Isokinetic {
    pnt* v1;
    pnt* v2;

    double* fric;
    double* Q1;
    double* Q2;
};

/**
* This is the structure for the various matrix elements. These specific variables
* will be discussed in the Force.cpp file.
*
**/

struct Matrix {
    double * M_di;
    double * M_di_inv;
    double * M_off_inv;
    double *K;
    double *KM;
};

/**
*   This is the main structure based on the system with a variety of variables. 
*
**/

struct System {
    int P;                          // # of PI beads                                                             
    pnt* position_p;                // Primitive position coordinates
    pnt* velocity_p;                // Primitive velocity coordinates
    pnt* force_p;                   // Primitive force coordinates
      
    pnt* position_t;                // Transformed position coordinates
    pnt* velocity_t;                // Transformed velocity coordinates
    pnt* force_t;                   // Transformed force coordinates 
    
    pnt* force_harm;                // Forces associated with bead spring              
    pnt* force_matrix_p;            // Primitive forces for external potential
    pnt* force_matrix_t;                 // Transformed coordinates for external pot.
    pnt* forces_total;
    
	double* mass_harm;				// Mass related to harmonic coupling
	double* mass_vel;               // Mass associated with velocities
	double* lambda;				    // Eigenvalue for staging transformation
	double beta;                    // 1/kbT
	double hbar;                    
	double sigma;                   // Dummy variable used for RNG
	double * U1, *U2, *U3;          // Derivatives of the PES


	double *energy;                 
	
	double t;                       // Time of the correlation function
	
	// See original paper for definition of these variables
	double tau_c;
	double gamma;
	double alpha;	
	double wp;
	
	double dt;                      // Timestep of simulation
	
	int nsteps;                     // Number of simulation steps
    Isokinetic isok;                // Definition of the Isokinetic thermostat
    
    // Variables needed for Isokinetic Thermostat
    int L;
    int nsy;
    int nres;
    double* wdt;
        
    Matrix m;                       // Declaration of the matrix structure
    
    
    
};





#endif //SYSTEM_H



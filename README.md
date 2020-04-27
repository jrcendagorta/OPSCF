To compile code:

make

To run code, requires 6 arguments

1) Number of path integral beads
2) Beta = 1/kbT
3) Time for the correlation function
4) Number of simulation steps
5) timestep of simulation  (dt)
6) Potential type (1=HO, 2=Mildly Anharmonic, 3= Quartic, 4= DW)

For example a simulation of the Quartic well, at beta=8, for t=3 a.u. with 64 beads would be

./OPSCF 64 8.0 3.0 100000000 0.01 3


Typically you will need around 10^7-10^8 steps with a time step less than 0.01 depending on the system. 
All the units are in these dummy natural units where k_b (boltzmann) = 1, hbar = 1. 

Two files are written out, one that writes out the instantaneous value of the correlation functions 
(x(0)*x(t) and x^2(0)*x^2(t)) and another file that calculates the running average for all points. 

It is recommended to not use all points for average but the average file gives an appropriate indication of if 
the correlation functions are converging as a function of the number of simulation steps. 


Example output files are provided for a simulation:

./OPSCF 8 1.0 3.0 100000 0.01 3

This is nowhere near enough simulation steps, but this will give you an idea of how to run the program

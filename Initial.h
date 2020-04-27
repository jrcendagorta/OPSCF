#ifndef INITIAL_H
#define INITIAL_H

double randGaus(double sigma = 1);
void Allocate(System *s, int& P);
void Deallocate(System *s);
void Initialize_Mass(System *s, double mass);
void Initialize_Position(System *s);
void Initialize_Velocity(System *s);
void Initialize(System *s, double beta, double mass, double timestep, double time, int pot_typ);
void TransformP_TtoP(System *s);
void TransformF_PtoT(System *s);
void Initialize_Isokinetic(System* s, int L);

#endif //INITIAL_H

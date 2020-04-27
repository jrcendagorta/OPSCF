#ifndef INTEGRATORS_H
#define INTEGRATORS_H

void Integrate_Velocity(System *s, double dt);
void Integrate_Position(System* s, double dt);
void Integrate_Nose(System* s, double dt);
void Integrate_OU(System* s, double dt);
#endif //INTEGRATORS_H

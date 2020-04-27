#ifndef FORCE_H
#define FORCE_H

void Force_Calc(System *s);
void Force_Control(System* s, int pot_typ);
void Harm_Calc(System* s);
void Matrix_Calc(System *s);
double InverseMarix(double **m, int n, double** mi, double mo, double a, double* det);
void Derivatives(System *s, int pot_type);

#endif //FORCE_H

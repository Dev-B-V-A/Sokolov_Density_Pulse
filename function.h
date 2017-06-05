#ifndef FUNCTION_H
#define FUNCTION_H

double pulse_1 (double t, double x, double y);
double pulse_2 (double t, double x, double y);
double density_ (double t, double x, double y);

double f0 (double t, double x, double y);
double f1 (double t, double x, double y, double mu);
double f2 (double t, double x, double y, double mu);

void print_array (double *a, int m, int n);

#endif  //FUNCTION_H

#ifndef __FUNCTIONS_H
#define __FUNCTIONS_H

#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>

#include "Types.h"

const double g = 1.4;   // Adiabatic constant
const int NGC = 2;       // Number of ghost cells per side

// Functions written for approximate solvers: 
void initialConditions(Primitive *W, Conserved *U, 
        double x_0, int N, double dx,
        double rho_L, double u_L, double p_L, 
        double rho_R, double u_R, double p_R);
double maxWaveSpeed(Primitive *W, double a, int N);
void PrimitiveToConserved(Primitive W, Conserved &U);
void ConservedToPrimitive(Conserved U, Primitive &W);
void boundaryConditions(Primitive *W, Conserved *U, int N);
void flux(Primitive W, Conserved U, Conserved &F);
void laxFriedrich(Primitive W_L, Primitive W_R, 
        Conserved U_L, Conserved U_R, double dt, double dx,
        Conserved &f);
void richtmeyer(Primitive W_L, Primitive W_R, 
        Conserved U_L, Conserved U_R, double dt, double dx,
        Conserved &f);
void force(Primitive W_L, Primitive W_R, 
        Conserved U_L, Conserved U_R, double dt, double dx,
        Conserved &f);
double minbee(double r); 
double vanleer(double r, double a_max);
double superbee(double r, double a_max);
double fluxLimiter(double E_min, double E_0, 
        double E_plus, double E_2plus, char *limitFunc);
void flic(Primitive W_L, Primitive W_R, 
        Conserved U_L, Conserved U_R, 
        Conserved U_2L, Conserved U_2R, 
        double dt, double dx, Conserved &f);
void hllc(Primitive W_L, Primitive W_R, Conserved &f);
void godunov(Primitive W_L, Primitive W_R, Conserved &f);
double sminbee(double r);
double svanleer(double r);
double ssuperbee(double r);
double slopeLimiter(double E_min, double E_0, 
        double E_plus, char *limitFunc);
void muscl(Conserved *U, 
        double dt, double dx, double omega, 
        int N, char *limitFunc, char *RP, 
        Conserved *f);
void slic(Conserved *U, 
        double dt, double dx, double omega, 
        int N, char *limitFunc, Conserved *f);
void advance(Primitive *W, Conserved *U, Conserved *U_old, 
        double dt, double dx, int N, 
        char *scheme, char *limitFunc, char *RP);

// Functions written for exact solver: 
void initiateTestCase(int testCase, double &tStop, double &x_0,
        double &rho_L, double &u_L, double &p_L, 
        double &rho_R, double &u_R, double &p_R);
void starRegionPressureVelocity(double &p, double &u,
        double rho_L, double u_L, double p_L, double a_L,
        double rho_R, double u_R, double p_R, double a_R);
void pressureFunctions(double &f, double &f_Diff, double p_Old, 
        double rho, double u, double p, double a);
void sample(double &rho, double &u, double &p,
        double p_S, double u_S, double S,
        double rho_L, double u_L, double p_L, double a_L,
        double rho_R, double u_R, double p_R, double a_R);

#endif

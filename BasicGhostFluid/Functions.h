#ifndef __FUNCTIONS_H
#define __FUNCTIONS_H

#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>

#include "Types.h"

const int NGC = 2;       // Number of ghost cells per side

void updateGhostCells(Primitive *W_A, Conserved *U_A, 
        Primitive *W_B, Conserved *U_B, double *phi, 
        double gamma_L, double gamma_R);

void advanceLevelSet(double *phi, Primitive *W_A, Primitive *W_B, 
        double dt, double dx, int N);

int interfacePosition(double *phi);

void reinitialize(double *phi, int N, double dx);

void initialConditions(Primitive *W_A, Conserved *U_A, 
        Primitive *W_B, Conserved *U_B, double *phi,
        double x_0, int N, double dx,
        double rho_L, double u_L, double p_L, double gamma_L,
        double rho_R, double u_R, double p_R, double gamma_R);

double maxWaveSpeed(Primitive *W, double a, int N);

void PrimitiveToConserved(Primitive W, Conserved &U,double gamma);

void ConservedToPrimitive(Conserved U, Primitive &W,double gamma);

void boundaryConditions(Primitive *W, Conserved *U, double gamma, 
        int N);

void flux(Primitive W, Conserved U, Conserved &F);

void laxFriedrich(Primitive W_L, Primitive W_R, 
        Conserved U_L, Conserved U_R, double dt, double dx,
        Conserved &f);

void richtmeyer(Primitive W_L, Primitive W_R, 
        Conserved U_L, Conserved U_R, double dt, double dx,
        double gamma, Conserved &f);

void force(Primitive W_L, Primitive W_R, 
        Conserved U_L, Conserved U_R, double dt, double dx,
        double gamma, Conserved &f);

double minbee(double r); 
double vanleer(double r, double a_max);
double superbee(double r, double a_max);
double fluxLimiter(double E_min, double E_0, 
        double E_plus, double E_2plus, char *limitFunc);

void flic(Primitive W_L, Primitive W_R, 
        Conserved U_L, Conserved U_R, 
        Conserved U_2L, Conserved U_2R, 
        double dt, double dx, double gamma, Conserved &f);

double sminbee(double r);
double svanleer(double r);
double ssuperbee(double r);
double slopeLimiter(double E_min, double E_0, 
        double E_plus, char *limitFunc);

void slic(Conserved *U, 
        double dt, double dx, double omega, 
        int N, char *limitFunc, double gamma, Conserved *f);

void advance(Primitive *W, Conserved *U, Conserved *U_old, 
        double dt, double dx, int N, double gamma, 
        char *scheme, char *limitFunc);

void initiateTestCase(int testCase, double &tStop, double &x_0,
        double &rho_L, double &u_L, double &p_L, double &gamma_L,
        double &rho_R, double &u_R, double &p_R, double &gamma_R);

void setScheme(char *scheme, char *limitFunc, 
        int schemeChoice, int limitChoice);
#endif

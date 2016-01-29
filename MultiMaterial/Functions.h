#ifndef __FUNCTIONS_H
#define __FUNCTIONS_H

#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>

#include "Types.h"

const int NGC = 2;       // Number of ghost cells per side

void updateGhostCells(Primitive *W_L, Conserved *U_L, 
        Primitive *W_R, Conserved *U_R, int pos, double gamma_L, 
        double gamma_R, double p_Inf_L, double p_Inf_R, 
        int method);

void hllcStarStates(Primitive W_L, Primitive W_R, Conserved U_L, 
        Conserved U_R, double gamma_L, double gamma_R, 
        double p_Inf_L, double p_Inf_R, Conserved &U_L_star, 
        Conserved &U_R_star);

void advanceLevelSet(double *phi, Primitive *W_A, Primitive *W_B, 
        double dt, double dx, int N);

int* interfacePosition(double *phi, int nMaterialInterfaces);

void reinitialize(double *phi, int N, double dx, 
        int nMaterialInterfaces);

void initialConditions(Primitive *W_A, Conserved *U_A, 
        Primitive *W_B, Conserved *U_B, double *phi, 
        int N, double dx, int nInterfaces, double *x, 
        int nMaterialInterfaces, double *x_M, double *rho, 
        double *u, double *p, double *gamma, 
        double *p_Inf);

double maxWaveSpeed(Primitive *W, double a, int N);

void PrimitiveToConserved(Primitive W, Conserved &U,double gamma,
        double p_Inf);

void ConservedToPrimitive(Conserved U, Primitive &W,double gamma, 
        double p_Inf);

void boundaryConditions(Primitive *W, Conserved *U, double gamma, 
        double p_Inf, int N);

void flux(Primitive W, Conserved U, Conserved &F);

void laxFriedrich(Primitive W_L, Primitive W_R, Conserved U_L, 
        Conserved U_R, double dt, double dx, Conserved &f);

void richtmeyer(Primitive W_L, Primitive W_R, Conserved U_L, 
        Conserved U_R, double dt, double dx, double gamma, 
        double p_Inf, Conserved &f);

void force(Primitive W_L, Primitive W_R, Conserved U_L, 
        Conserved U_R, double dt, double dx, double gamma, 
        double p_Inf, Conserved &f);

double minbee(double r); 
double vanleer(double r, double a_max);
double superbee(double r, double a_max);
double fluxLimiter(double E_min, double E_0, 
        double E_plus, double E_2plus, char *limitFunc);

void flic(Primitive W_L, Primitive W_R, Conserved U_L, 
        Conserved U_R, Conserved U_2L, Conserved U_2R, 
        double dt, double dx, char *limitFunc, double gamma, 
        double p_Inf, Conserved &f);

double sminbee(double r);
double svanleer(double r);
double ssuperbee(double r);
double slopeLimiter(double E_min, double E_0, double E_plus, 
        char *limitFunc);

void slic2(Conserved U_L, Conserved U_0, Conserved U_R, 
        Conserved U_2R, double dt, double dx, char *limitFunc, 
        double gamma, double p_Inf, Conserved &f);

void advance(Primitive *W, Conserved *U, Conserved *U_old, 
        double dt, double dx, int N, double gamma, double p_Inf, 
        char *scheme, char *limitFunc);

void initiateTestCase(int testCase, int &nInterfaces, 
        double x[4], double rho[5], double u[5], double p[5], 
        int mat[5], double gamma[2], double p_Inf[2], 
        double &tStop);

void setScheme(char *scheme, char *limitFunc, 
        int schemeChoice, int limitChoice);

void starRegionPressureVelocity(double &p, double &u,
        double rho_L, double u_L, double p_L, double a_L, 
        double gamma_L, double p_Inf_L, double rho_R, double u_R, 
        double p_R, double a_R, double gamma_R, double p_Inf_R);

void pressureFunctions(double &f, double &f_Diff, double p_Old, 
        double rho, double p, double a, double g, 
        double p_Inf);

void starRegionRho(double &rho_L_S, double &rho_R_S, double p_S, 
        double rho_L, double p_L, double g_L, double p_Inf_L,
        double rho_R, double p_R, double g_R, double p_Inf_R);

void sample(double &rho, double &u, double &p, double &e,
        double p_S, double u_S, double S,
        double rho_L, double u_L, double p_L, double a_L, 
        double g_L, double p_Inf_L, double rho_R, double u_R, 
        double p_R, double a_R, double g_R, double p_Inf_R);

#endif

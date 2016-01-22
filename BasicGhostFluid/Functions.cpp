#ifndef __FUNCTIONS_CPP
#define __FUNCTIONS_CPP

#include "Functions.h"

void updateGhostCells(Primitive *W_A, Conserved *U_A, 
        Primitive *W_B, Conserved *U_B, double *phi, 
        double gamma_L, double gamma_R){
    int pos = interfacePosition(phi);
    // Isobaric fix: 
/*    W_A[pos-1].rho = pow(W_A[pos-1].p/W_A[pos-2].p, 1/gamma_L)
            *W_A[pos-2].rho;
    W_B[pos].rho = pow(W_B[pos].p/W_A[pos+1].p, 1/gamma_R)
            *W_A[pos+1].rho;
    PrimitiveToConserved(W_A[pos-1], U_A[pos-1], gamma_L);
    PrimitiveToConserved(W_B[pos], U_B[pos], gamma_R);*/
    // Ghost fluid update, constant entropy: 
    for(int i=0; i<2; i++){
        W_A[pos+i].rho = pow(W_B[pos+i].p/W_A[pos-1].p, 1/gamma_L)
            *W_A[pos-1].rho;
        W_A[pos+i].u = W_B[pos+i].u; 
        W_A[pos+i].p = W_B[pos+i].p;
        W_B[pos-1-i].rho = 
            pow(W_A[pos-1-i].p/W_B[pos].p, 1/gamma_R)
            *W_B[pos].rho;
        W_B[pos-1-i].u = W_A[pos-1-i].u; 
        W_B[pos-1-i].p = W_A[pos-1-i].p;
        PrimitiveToConserved(W_A[pos+i], U_A[pos+i], gamma_L);
        PrimitiveToConserved(W_B[pos-1-i], U_B[pos-1-i], gamma_R); 
    }
}

void advanceLevelSet(double *phi, Primitive *W_A, Primitive *W_B, 
        double dt, double dx, int N){
    int pos = interfacePosition(phi); 
    double tmp[N+2*NGC]; 
    for(int i=0; i<N+2*NGC; i++){
        tmp[i] = phi[i]; 
    }
    for(int i=NGC; i<NGC+N; i++){
        if(i<pos){
            if(W_A[i].u>0.0)
                phi[i]=tmp[i]-dt/dx*W_A[i].u*(tmp[i]-tmp[i-1]);
            else
                phi[i]=tmp[i]-dt/dx*W_A[i].u*(tmp[i+1]-tmp[i]);
        }else{
            if(W_B[i].u>0.0)
                phi[i]=tmp[i]-dt/dx*W_B[i].u*(tmp[i]-tmp[i-1]);
            else
                phi[i]=tmp[i]-dt/dx*W_B[i].u*(tmp[i+1]-tmp[i]);
        }
    }
}

int interfacePosition(double *phi){
    int pos = 0;
    while(phi[pos]<0.0) pos++;
    return pos;
}

void reinitialize(double *phi, int N, double dx){
    int pos = interfacePosition(phi);
    for(int i=0; i<pos-1; i++){
        phi[i] = phi[pos-1]-(pos-1-i)*dx;
    }
    for(int i=pos+1; i<N+2*NGC; i++){
        phi[i] = phi[pos]+(i-pos)*dx;
    }
}

void initialConditions(Primitive *W_A, Conserved *U_A, 
        Primitive *W_B, Conserved *U_B, double *phi, 
        double x_0, int N, double dx, 
        double rho_L, double u_L, double p_L, double gamma_L,
        double rho_R, double u_R, double p_R, double gamma_R){
   for(int i=0; i<N+2*NGC; i++){
       phi[i] = (i-NGC+0.5)*dx - x_0;
//       if((i-NGC+0.5)*dx <= x_0){
           W_A[i].rho = rho_L; 
           W_A[i].u = u_L; 
           W_A[i].p = p_L; 
//       }else{
           W_B[i].rho = rho_R; 
           W_B[i].u = u_R; 
           W_B[i].p = p_R; 
//       }
       PrimitiveToConserved(W_A[i], U_A[i], gamma_L);
       PrimitiveToConserved(W_B[i], U_B[i], gamma_R);
   } 
   reinitialize(phi, N, dx);
}

double maxWaveSpeed(Primitive *W, double a, int N){
    // Smax = max_i {abs(S_L_i+0.5), abs(S_R_i+0.5)}
    double u_max = 0;
    for(int i=0; i<N+2*NGC; i++){
        if(u_max < fabs(W[i].u)) u_max = fabs(W[i].u);
    }
    return u_max+a; 
}

void PrimitiveToConserved(Primitive W, Conserved &U, double gamma){
    U.rho = W.rho;
    U.rho_u = W.rho*W.u; 
    U.E = 0.5*W.rho*W.u*W.u + W.p/(gamma-1);
}

void ConservedToPrimitive(Conserved U, Primitive &W, double gamma){
    W.rho = U.rho;
    W.u = U.rho_u/U.rho; 
    W.p = (U.E-0.5*U.rho_u*U.rho_u/U.rho)*(gamma-1);
}

void boundaryConditions(Primitive *W, Conserved *U, double gamma, 
        int N){
    // Transmissive boundary conditions
    for(int i=0; i<NGC; i++){
        W[i] = W[i+1];
        W[N+NGC-i] = W[N+NGC-i-1];
        PrimitiveToConserved(W[i], U[i], gamma);
        PrimitiveToConserved(W[N+NGC-i], U[N+NGC-i], gamma); 
    }
}

void flux(Primitive W, Conserved U, Conserved &F){
    F.rho = U.rho_u;
    F.rho_u = U.rho_u*W.u+W.p;
    F.E = W.u*(U.E+W.p);
}

void laxFriedrich(Primitive W_L, Primitive W_R, 
        Conserved U_L, Conserved U_R, double dt, double dx,
        Conserved &f){
    Conserved F_L, F_R, f_LF;
    flux(W_L, U_L, F_L); 
    flux(W_R, U_R, F_R); 
    f = 0.5*(F_L+F_R
            +dx/dt*(U_L-U_R));
}

void richtmeyer(Primitive W_L, Primitive W_R, 
        Conserved U_L, Conserved U_R, double dt, double dx,
        double gamma, Conserved &f){
    Conserved F_L, F_R, U_half;
    Primitive W_half;
    flux(W_L, U_L, F_L); 
    flux(W_R, U_R, F_R); 
    U_half = 0.5*((U_L+U_R)
            +dt/dx*(F_L-F_R));
    ConservedToPrimitive(U_half, W_half, gamma);
    flux(W_half, U_half, f);
}

void force(Primitive W_L, Primitive W_R, 
        Conserved U_L, Conserved U_R, double dt, double dx,
        double gamma, Conserved &f){
    Conserved f_LF, f_RI;
    laxFriedrich(W_L, W_R, U_L, U_R, dt,dx, f_LF);
    richtmeyer(W_L, W_R, U_L, U_R, dt, dx, gamma, f_RI);
    f = 0.5*(f_LF+f_RI);
}

double minbee(double r){
    if(r<=0.0) return 0.0;
    else if(r<=1.0) return r;
    else return 1.0;
}

double vanleer(double r, double a_max){
    double phi_g = (1.0-a_max)/(1.0+a_max);
    if(r<=0.0) return 0.0;
    else if(r<=1.0) return 2.0*r/(1.0+r);
    else return phi_g+2*(1.0-phi_g)*r/(1.0+r);
}

double superbee(double r, double a_max){
    double phi_g = (1.0-a_max)/(1.0+a_max);
    if(r<=0.0) return 0.0;
    else if(r<=0.5) return 2.0*r;
    else if(r<=1.0) return 1.0; 
    else return (2.0 < phi_g+(1.0-phi_g)*r) ? 
        2.0 : phi_g+(1.0-phi_g)*r;
}

double fluxLimiter(double q_min, double q_0, 
        double q_plus, double q_2plus, char *limitFunc){
    double Delta_0, Delta_L, Delta_R, tol = 1e-12;
    if(fabs(q_plus-q_0)<tol){
        Delta_0 = copysign(tol, q_plus-q_0);
    }else{
        Delta_0 = q_plus-q_0;
    }
    if(fabs(q_0-q_min)<tol){
        Delta_L = copysign(tol, q_0-q_min);
    }else{
        Delta_L = q_0-q_min;
    }
    if(fabs(q_2plus-q_plus)<tol){
        Delta_R = copysign(tol, q_2plus-q_plus);
    }else{
        Delta_R = q_2plus-q_plus;
    }
    double r_L = Delta_L/Delta_0;
    double r_R = Delta_R/Delta_0;
    double phi_L, phi_R;
    if(!strcmp(limitFunc, "minbee")){
        phi_L = minbee(r_L);
        phi_R = minbee(r_R);
    }else if(!strcmp(limitFunc, "vanleer")){
        phi_L = vanleer(r_L, 0.9);
        phi_R = vanleer(r_R, 0.9);
    }else{ 
        phi_L = superbee(r_L, 0.9);
        phi_R = superbee(r_R, 0.9);
    }
    return (phi_L < phi_R) ? phi_L : phi_R;
}

void flic(Primitive W_L, Primitive W_R, 
        Conserved U_L, Conserved U_R, 
        Conserved U_2L, Conserved U_2R, 
        double dt, double dx, char *limitFunc, 
        double gamma, Conserved &f){
    Conserved f_RI, f_FORCE;
    double phi = fluxLimiter(U_2L.rho, U_L.rho, 
            U_R.rho, U_2R.rho, limitFunc);
    richtmeyer(W_L, W_R, U_L, U_R, dt, dx, gamma, f_RI);
    force(W_L, W_R, U_L, U_R, dt, dx, gamma, f_FORCE);
    f = f_FORCE+phi*(f_RI-f_FORCE);
}

double sminbee(double r){
    if(r <= 0.0) return 0.0;
    else if(r<=1.0) return r;
    else return (1.0<2.0/(1.0+r)) ? 1.0 : 2.0/(1.0+r);
}
double svanleer(double r){
    if(r<=0.0) return 0.0;
    else if(r<=1.0) return 2.0*r/(1.0+r);
    else return 2.0/(1.0+r);
}

double ssuperbee(double r){
    if(r<=0.0) return 0.0;
    else if(r<=0.5) return 2.0*r;
    else if(r<=1.0) return 1.0; 
    else return (2.0 < r) ? 
        (2.0<2.0/(1.0+r) ? 2.0 : 2.0/(1.0+r)) : 
            (r<2.0/(1.0+r) ? r : 2.0/(1.0+r));
}

double slopeLimiter(double q_min, double q_0, 
        double q_plus, char *limitFunc){
    double Delta_L, Delta_R, tol = 1e-12;
    if(fabs(q_0-q_min)<tol){
        Delta_L = copysign(tol, q_0-q_min);
    }else{
        Delta_L = q_0-q_min;
    }
    if(fabs(q_plus-q_0)<tol){
        Delta_R = copysign(tol, q_plus-q_0);
    }else{
        Delta_R = q_plus-q_0;
    }
    double r = Delta_L/Delta_R;
    if(!strcmp(limitFunc, "minbee")){
        return sminbee(r);
    }else if(!strcmp(limitFunc, "vanleer")){
        return svanleer(r);
    }else if(!strcmp(limitFunc, "superbee")){
        return ssuperbee(r);
    }
    printf("fail!");
    return 1337;
}

void slic(Conserved *U, 
        double dt, double dx, double omega, 
        int N, char *limitFunc, double gamma, Conserved *f){
    Primitive *W_L = new Primitive[N+2*NGC],
              *W_R = new Primitive[N+2*NGC],
              *W_L_bar = new Primitive[N+2*NGC],
              *W_R_bar = new Primitive[N+2*NGC];
    Conserved *Delta = new Conserved[N+2*NGC],
              *U_L = new Conserved[N+2*NGC],
              *U_R = new Conserved[N+2*NGC],
              *F_L = new Conserved[N+2*NGC],
              *F_R = new Conserved[N+2*NGC],
              *U_L_bar = new Conserved[N+2*NGC],
              *U_R_bar = new Conserved[N+2*NGC];
    int L, R;
    double xi;
    for(int i=NGC-1; i<N+NGC+1; i++){
        Delta[i] = 0.5*(1+omega)*(U[i]-U[i-1])
                +0.5*(1-omega)*(U[i+1]-U[i]);
        xi = slopeLimiter(U[i-1].rho, U[i].rho, 
                U[i+1].rho, limitFunc); 
        Delta[i] = xi*Delta[i]; 
        U_L[i] = U[i]-0.5*Delta[i];
        U_R[i] = U[i]+0.5*Delta[i]; 
        ConservedToPrimitive(U_L[i], W_L[i], gamma);
        ConservedToPrimitive(U_R[i], W_R[i], gamma);
        flux(W_L[i], U_L[i], F_L[i]);
        flux(W_R[i], U_R[i], F_R[i]);
        U_L_bar[i] = U_L[i]+0.5*dt/dx*(F_L[i]-F_R[i]);
        U_R_bar[i] = U_R[i]+0.5*dt/dx*(F_L[i]-F_R[i]);
        ConservedToPrimitive(U_L_bar[i], W_L_bar[i], gamma);
        ConservedToPrimitive(U_R_bar[i], W_R_bar[i], gamma);
    }
    for(int i=0; i<N+1; i++){
        L = i+NGC-1; R = i+NGC;
        force(W_R_bar[L], W_L_bar[R], U_R_bar[L], U_L_bar[R], 
                dt, dx, gamma, f[i]); 
    }
    delete []W_L; W_L=NULL;
    delete []W_R; W_R=NULL;
    delete []W_L_bar; W_L_bar=NULL;
    delete []W_R_bar; W_R_bar=NULL;
    delete []U_L; U_L=NULL;
    delete []U_R; U_R=NULL;
    delete []F_L; F_L=NULL;
    delete []F_R; F_R=NULL;
    delete []U_L_bar; U_L_bar=NULL;
    delete []U_R_bar; U_R_bar=NULL;
}

void advance(Primitive *W, Conserved *U, Conserved *U_old, 
        double dt, double dx, int N, double gamma, 
        char *scheme, char *limitFunc){
    Conserved *f = new Conserved[N+1];
    int L, R;
    if(!strcmp(scheme, "SLIC"))
        slic(U_old, dt, dx, 0.0, N, limitFunc, gamma, f); 
    else{
        for(int i=0; i<N+1; i++){
            L = i+NGC-1; R = i+NGC;
            if(!strcmp(scheme, "FORCE"))
                force(W[L], W[R], U_old[L], U_old[R], 
                        dt, dx, gamma, f[i]);
            else if(!strcmp(scheme, "FLIC"))
                flic(W[L], W[R], U_old[L], U_old[R], U_old[L-1], 
                        U_old[R+1], dt, dx, limitFunc, gamma,f[i]);
        }
    }
    for(int i=NGC; i<N+NGC; i++){
        L = i-NGC; R = i-NGC+1;
        U[i] = U_old[i]-dt/dx*(f[R]-f[L]); 
        ConservedToPrimitive(U[i], W[i], gamma); 
    }
    for(int i=NGC; i<N+NGC; i++){
        U_old[i] = U[i];
    }
    delete []f; f = NULL;
}

void initiateTestCase(int testCase, double &tStop, double &x_0,
        double &rho_L, double &u_L, double &p_L, double &gamma_L,
        double &rho_R, double &u_R, double &p_R, double &gamma_R){
    if(testCase==1){
        rho_L = 1.0;
        u_L = 0.75;
        p_L = 1.0;
        gamma_L = 1.4;
        rho_R = 0.125;
        u_R = 0.0;
        p_R = 0.1;
        gamma_R = 1.4; 
        x_0 = 0.3;
        tStop = 0.2; 
    }else if(testCase==2){
        rho_L = 1.0;
        u_L = -2.0;
        p_L = 0.4;
        gamma_L = 1.4;
        rho_R = 1.0;
        u_R = 2.0;
        p_R = 0.4;
        gamma_R = 1.4; 
        x_0 = 0.5;
        tStop = 0.15; 
    }else if(testCase==3){
        rho_L = 1.0;
        u_L = 0.0;
        p_L = 1000;
        gamma_L = 1.4;
        rho_R = 1.0;
        u_R = 0.0;
        p_R = 0.01;
        gamma_R = 1.4; 
        x_0 = 0.5;
        tStop = 0.012; 
    }else if(testCase==4){
        rho_L = 5.99924;
        u_L = 19.5975;
        p_L = 460.894;
        gamma_L = 1.4;
        rho_R = 5.9942;
        u_R = -6.19633;
        p_R = 46.0950;
        gamma_R = 1.4; 
        x_0 = 0.4;
        tStop = 0.035; 
    }else if(testCase==5){
        rho_L = 1.0;
        u_L = -19.59745;
        p_L = 1000;
        gamma_L = 1.4;
        rho_R = 1.0;
        u_R = -19.59745;
        p_R = 0.01;
        gamma_R = 1.4; 
        x_0 = 0.8;
        tStop = 0.012; 
    }else if(testCase==6){
        rho_L = 1.0;
        u_L = 0.5;
        p_L = 1.0;
        gamma_L = 1.4;
        rho_R = 0.138;
        u_R = 0.5;
        p_R = 1.0;
        gamma_R = 1.4; 
        x_0 = 0.25;
        tStop = 1.0; 
    }else if(testCase==7){
        rho_L = 1.0;
        u_L = 0.5;
        p_L = 1.0;
        gamma_L = 1.4;
        rho_R = 0.138;
        u_R = 0.5;
        p_R = 1.0;
        gamma_R = 1.67; 
        x_0 = 0.25;
        tStop = 1.0; 
    }
}

void setScheme(char *scheme, char *limitFunc, 
        int schemeChoice, int limitChoice){
    if(schemeChoice == 0)
        strcpy(scheme, "FORCE");
    else if(schemeChoice == 1){
        strcpy(scheme, "FLIC");
        if(limitChoice == 0)
            strcpy(limitFunc, "minbee");
        else if(limitChoice == 1)
            strcpy(limitFunc, "vanleer");
        else
            strcpy(limitFunc, "superbee");
    }else{
        strcpy(scheme, "SLIC"); 
        if(limitChoice == 0)
            strcpy(limitFunc, "minbee");
        else if(limitChoice == 1)
            strcpy(limitFunc, "vanleer");
        else
            strcpy(limitFunc, "superbee");
    }
}

#endif

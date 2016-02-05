#ifndef __FUNCTIONS_CPP
#define __FUNCTIONS_CPP

#include "Functions.h"
#include <QDebug>

void updateGhostCells(Primitive *W_L, Conserved *U_L, 
        Primitive *W_R, Conserved *U_R, int pos, double gamma_L, 
        double gamma_R, double p_Inf_L, double p_Inf_R, 
        int method)
{
    if(method>0)
    {
        // Basic ghost fluid method
        if(method>1)
        {
            // Isobaric fix: 
            W_L[pos-1].rho = 
                pow(W_L[pos-1].p/W_L[pos-2].p, 1/gamma_L)
                *W_L[pos-2].rho;
            W_R[pos].rho = pow(W_R[pos].p/W_R[pos+1].p, 1/gamma_R)
                    *W_R[pos+1].rho;
            PrimitiveToConserved(W_L[pos-1], U_L[pos-1], gamma_L, 
                    p_Inf_L);
            PrimitiveToConserved(W_R[pos], U_R[pos], gamma_R, 
                    p_Inf_R);
        }
        // Ghost fluid update, constant entropy: 
        for(int i=0; i<2; i++)
        {
            W_L[pos+i].rho = 
                pow(W_R[pos+i].p/W_L[pos-1].p, 1/gamma_L)
                *W_L[pos-1].rho;
            W_L[pos+i].u = W_R[pos+i].u; 
            W_L[pos+i].p = W_R[pos+i].p;
            W_R[pos-1-i].rho = 
                pow(W_L[pos-1-i].p/W_R[pos].p, 1/gamma_R)
                *W_R[pos].rho;
            W_R[pos-1-i].u = W_L[pos-1-i].u; 
            W_R[pos-1-i].p = W_L[pos-1-i].p;
            PrimitiveToConserved(W_L[pos+i], U_L[pos+i], gamma_L, 
                    p_Inf_L);
            PrimitiveToConserved(
                    W_R[pos-1-i], U_R[pos-1-i], gamma_R, p_Inf_R);
        }
    }
    else
    {
        // Riemann ghost fluid method, Exact solver
        double a_L = 
            sqrt(gamma_L*(W_L[pos-2].p+p_Inf_L)/W_L[pos-2].rho);
        double a_R = 
            sqrt(gamma_R*(W_R[pos+1].p+p_Inf_R)/W_R[pos+1].rho);

        double p_Star, u_Star; 
        starRegionPressureVelocity(p_Star, u_Star, 
              W_L[pos-2].rho, W_L[pos-2].u, W_L[pos-2].p, a_L, 
              gamma_L, p_Inf_L, W_R[pos+1].rho, W_R[pos+1].u, 
              W_R[pos+1].p, a_R, gamma_R, p_Inf_R);

        double rho_L_Star, rho_R_Star;
        starRegionRho(rho_L_Star, rho_R_Star, p_Star, 
                W_L[pos-2].rho, W_L[pos-2].p, gamma_L, p_Inf_L,
                W_R[pos+1].rho, W_R[pos+1].p, gamma_R, p_Inf_R); 

        Primitive W_L_star(rho_L_Star, u_Star, p_Star);
        Primitive W_R_star(rho_R_Star, u_Star, p_Star);

        for(int i=0; i<3; i++)
        {
            W_L[pos-1+i] = W_L_star; 
            W_R[pos-i] = W_R_star; 
            PrimitiveToConserved(
                    W_L[pos-1+i], U_L[pos-1+i], gamma_L, p_Inf_L);
            PrimitiveToConserved(
                    W_R[pos-i], U_R[pos-i], gamma_R, p_Inf_R);
        }
    }
}

void advanceLevelSet(double *phi, Primitive *W_A, Primitive *W_B, 
        double dt, double dx, int N)
{
    double tmp[N+2*NGC]; 
    for(int i=0; i<N+2*NGC; i++)
        tmp[i] = phi[i]; 

    for(int i=NGC; i<NGC+N; i++)
    {
        if(phi[i]<0.0)
        {
            if(W_A[i].u>0.0)
                phi[i]=tmp[i]-dt/dx*W_A[i].u*(tmp[i]-tmp[i-1]);
            else
                phi[i]=tmp[i]-dt/dx*W_A[i].u*(tmp[i+1]-tmp[i]);
        }
        else
        {
            if(W_B[i].u>0.0)
                phi[i]=tmp[i]-dt/dx*W_B[i].u*(tmp[i]-tmp[i-1]);
            else
                phi[i]=tmp[i]-dt/dx*W_B[i].u*(tmp[i+1]-tmp[i]);
        }
    }
}

int* interfacePosition(double *phi, int nMaterialInterfaces)
{
    int* pos = new int[nMaterialInterfaces];
    pos[0]=NGC;
    for(int i=0; i<nMaterialInterfaces; i++)
    {
        if(i>0) pos[i]=pos[i-1]+1;
        while((phi[pos[i]]>0.0) == (pow(-1.0,i)<0.0))
                pos[i]++;
    }
    return pos;
}

void reinitialize(double *phi, int N, double dx, 
        int nMaterialInterfaces)
{
    int *pos = interfacePosition(phi, nMaterialInterfaces);
    //TODO: MAKE SURE phi does not become zero. 

    for(int i=0; i<pos[0]-1; i++)
        phi[i] = phi[pos[0]-1]-(pos[0]-1-i)*dx;

    for(int j=0; j<nMaterialInterfaces-1; j++)
    {
        for(int i=pos[j]+1; i<pos[j+1]-1; i++)
        {
            if((i-pos[j]) < (pos[j+1]-i))
                phi[i] = phi[pos[j]]
                    +pow(-1,j)*(i-pos[j])*dx; 
            else
                phi[i] = phi[pos[j+1]-1]
                    +pow(-1,j)*(pos[j+1]-1-i)*dx; 
        }
    }

    for(int i=pos[nMaterialInterfaces-1]+1; i<N+2*NGC; i++)
        phi[i] = phi[pos[nMaterialInterfaces-1]]
            +pow(-1,nMaterialInterfaces-1)
            *(i-pos[nMaterialInterfaces-1])*dx;

    delete[] pos;
}

void initialConditions(Primitive *W_A, Conserved *U_A, 
        Primitive *W_B, Conserved *U_B, double *phi, 
        int N, double dx, int nInterfaces, double *x, 
        int nMaterialInterfaces, double *x_M, double *rho, 
        double *u, double *p, double *gamma, 
        double *p_Inf)
{
    int region=0; int materialRegion=0; 
    for(int i=0; i<N+2*NGC; i++)
    {
        if((i-NGC+0.5)*dx>x[region] 
                && region<nInterfaces)
            region++;
        W_A[i].rho = rho[region]; 
        W_A[i].u = u[region]; 
        W_A[i].p = p[region]; 
        W_B[i].rho = rho[region]; 
        W_B[i].u = u[region]; 
        W_B[i].p = p[region]; 
        PrimitiveToConserved(W_A[i], U_A[i], gamma[0], p_Inf[0]);
        PrimitiveToConserved(W_B[i], U_B[i], gamma[1], p_Inf[1]);

        if((i-NGC+0.5)*dx > 
                (x_M[materialRegion]+x_M[materialRegion+1])/2
                && materialRegion<nMaterialInterfaces-1)
            materialRegion++;
        phi[i] = pow(-1,materialRegion)
            *((i-NGC+0.5)*dx-x_M[materialRegion]);
    } 
    reinitialize(phi, N, dx, nMaterialInterfaces);
}

double maxWaveSpeed(Primitive *W, double gamma, double p_Inf, 
        int N)
{
    // Smax = max_i {abs(u_i^n)+a_i^n}
    double a; 
    double S_max = 0;
    for(int i=0; i<N+2*NGC; i++)
    {
        a = sqrt(gamma*(W[i].p+p_Inf)/W[i].rho); 
        if(fabs(W[i].u)+a > S_max)
                S_max = fabs(W[i].u)+a;
    }
    return S_max;
}

void PrimitiveToConserved(Primitive W, Conserved &U, double gamma,
        double p_Inf)
{
    U.rho = W.rho;
    U.rho_u = W.rho*W.u; 
    U.E = 0.5*W.rho*W.u*W.u+(W.p+gamma*p_Inf)/(gamma-1);
}

void ConservedToPrimitive(Conserved U, Primitive &W, double gamma,
        double p_Inf)
{
    W.rho = U.rho;
    W.u = U.rho_u/U.rho; 
    W.p = (U.E-0.5*U.rho_u*U.rho_u/U.rho)*(gamma-1)-gamma*p_Inf;
}

void boundaryConditions(Primitive *W, Conserved *U, double gamma, 
        double p_Inf, int N)
{
    // Transmissive boundary conditions
    for(int i=0; i<NGC; i++)
    {
        W[i] = W[i+1];
        W[N+NGC-i] = W[N+NGC-i-1];
        PrimitiveToConserved(W[i], U[i], gamma, p_Inf);
        PrimitiveToConserved(W[N+NGC-i], U[N+NGC-i], gamma, 
                p_Inf); 
    }
}

void flux(Primitive W, Conserved U, Conserved &F)
{
    F.rho = U.rho_u;
    F.rho_u = U.rho_u*W.u+W.p;
    F.E = W.u*(U.E+W.p);
}

void laxFriedrich(Primitive W_L, Primitive W_R, Conserved U_L, 
        Conserved U_R, double dt, double dx, Conserved &f)
{
    Conserved F_L, F_R; 
    flux(W_L, U_L, F_L); 
    flux(W_R, U_R, F_R); 
    f = 0.5*(F_L+F_R+dx/dt*(U_L-U_R));
}

void richtmeyer(Primitive W_L, Primitive W_R, Conserved U_L, 
        Conserved U_R, double dt, double dx, double gamma, 
        double p_Inf, Conserved &f)
{
    Conserved F_L, F_R; 
    flux(W_L, U_L, F_L); 
    flux(W_R, U_R, F_R); 
    Conserved U_half = 0.5*((U_L+U_R)
            +dt/dx*(F_L-F_R));
    Primitive W_half;
    ConservedToPrimitive(U_half, W_half, gamma, p_Inf);
    flux(W_half, U_half, f);
}

void force(Primitive W_L, Primitive W_R, Conserved U_L, 
        Conserved U_R, double dt, double dx, double gamma, 
        double p_Inf, Conserved &f)
{
    Conserved f_LF, f_RI;
    laxFriedrich(W_L, W_R, U_L, U_R, dt,dx, f_LF);
    richtmeyer(W_L, W_R, U_L, U_R, dt, dx, gamma, p_Inf, f_RI);
    f = 0.5*(f_LF+f_RI);
}

double minbee(double r)
{
    if(r<=0.0) return 0.0;
    else if(r<=1.0) return r;
    else return 1.0;
}

double vanleer(double r, double a_max)
{
    double phi_g = (1.0-a_max)/(1.0+a_max);
    if(r<=0.0) return 0.0;
    else if(r<=1.0) return 2.0*r/(1.0+r);
    else return phi_g+2*(1.0-phi_g)*r/(1.0+r);
}

double superbee(double r, double a_max)
{
    double phi_g = (1.0-a_max)/(1.0+a_max);
    if(r<=0.0) return 0.0;
    else if(r<=0.5) return 2.0*r;
    else if(r<=1.0) return 1.0; 
    else return (2.0 < phi_g+(1.0-phi_g)*r) ? 
        2.0 : phi_g+(1.0-phi_g)*r;
}

double fluxLimiter(double q_min, double q_0, 
        double q_plus, double q_2p, char *limitFunc)
{
    double Delta_0, Delta_L, Delta_R; 
    double tol = 1e-12;

    if(fabs(q_plus-q_0)<tol) Delta_0 = copysign(tol, q_plus-q_0);
    else Delta_0 = q_plus-q_0;

    if(fabs(q_0-q_min)<tol) Delta_L = copysign(tol, q_0-q_min);
    else Delta_L = q_0-q_min;

    if(fabs(q_2p-q_plus)<tol) Delta_R =copysign(tol, q_2p-q_plus);
    else Delta_R = q_2p-q_plus;
    
    double r_L = Delta_L/Delta_0;
    double r_R = Delta_R/Delta_0;
    double phi_L, phi_R;
    if(!strcmp(limitFunc, "minbee"))
    {
        phi_L = minbee(r_L);
        phi_R = minbee(r_R);
    }
    else if(!strcmp(limitFunc, "vanleer"))
    {
        phi_L = vanleer(r_L, 0.9);
        phi_R = vanleer(r_R, 0.9);
    }
    else
    { 
        phi_L = superbee(r_L, 0.9);
        phi_R = superbee(r_R, 0.9);
    }
    return (phi_L < phi_R) ? phi_L : phi_R;
}

void flic(Primitive W_L, Primitive W_R, Conserved U_L, 
        Conserved U_R, Conserved U_2L, Conserved U_2R, 
        double dt, double dx, char *limitFunc, double gamma, 
        double p_Inf, Conserved &f)
{
    Conserved f_RI, f_FORCE;
    double phi = fluxLimiter(U_2L.rho, U_L.rho, U_R.rho, U_2R.rho,
            limitFunc);
    richtmeyer(W_L, W_R, U_L, U_R, dt, dx, gamma, p_Inf, f_RI);
    force(W_L, W_R, U_L, U_R, dt, dx, gamma, p_Inf, f_FORCE);
    f = f_FORCE+phi*(f_RI-f_FORCE);
}

double sminbee(double r)
{
    if(r <= 0.0) return 0.0;
    else if(r<=1.0) return r;
    else return (1.0<2.0/(1.0+r)) ? 1.0 : 2.0/(1.0+r);
}

double svanleer(double r)
{
    if(r<=0.0) return 0.0;
    else if(r<=1.0) return 2.0*r/(1.0+r);
    else return 2.0/(1.0+r);
}

double ssuperbee(double r)
{
    if(r<=0.0) return 0.0;
    else if(r<=0.5) return 2.0*r;
    else if(r<=1.0) return 1.0; 
    else return (2.0 < r) ? 
        (2.0<2.0/(1.0+r) ? 2.0 : 2.0/(1.0+r)) : 
            (r<2.0/(1.0+r) ? r : 2.0/(1.0+r));
}

double slopeLimiter(double q_min, double q_0, 
        double q_plus, char *limitFunc)
{
    double Delta_L, Delta_R; 
    double tol = 1e-12;

    if(fabs(q_0-q_min)<tol) Delta_L = copysign(tol, q_0-q_min);
    else Delta_L = q_0-q_min;

    if(fabs(q_plus-q_0)<tol) Delta_R = copysign(tol, q_plus-q_0);
    else Delta_R = q_plus-q_0;
    
    double r = Delta_L/Delta_R;
    if(!strcmp(limitFunc, "minbee")) return sminbee(r);
    else if(!strcmp(limitFunc, "vanleer")) return svanleer(r);
    else return ssuperbee(r);
}

void slic2(Conserved U_L, Conserved U_0, Conserved U_R, 
        Conserved U_2R, double dt, double dx, char *limitFunc, 
        double gamma, double p_Inf, Conserved &f)
{
    // Quicker version of SLIC. Omega is always equal to zero. 
    double xi_plus = 
        slopeLimiter(U_0.rho, U_R.rho, U_2R.rho, limitFunc);
    double xi_0 = 
        slopeLimiter(U_L.rho, U_0.rho, U_R.rho, limitFunc);

    Conserved U_L_plus = U_R-0.25*xi_plus*(U_2R-U_0);
    Conserved U_R_plus = U_R+0.25*xi_plus*(U_2R-U_0);
    Conserved U_L_0 = U_0-0.25*xi_0*(U_R-U_L);
    Conserved U_R_0 = U_0+0.25*xi_0*(U_R-U_L);

    Primitive W_L_plus, W_R_plus, W_L_0, W_R_0; 
    ConservedToPrimitive(U_L_plus, W_L_plus, gamma, p_Inf);
    ConservedToPrimitive(U_R_plus, W_R_plus, gamma, p_Inf);
    ConservedToPrimitive(U_L_0, W_L_0, gamma, p_Inf);
    ConservedToPrimitive(U_R_0, W_R_0, gamma, p_Inf);

    Conserved F_L_plus, F_R_plus, F_L_0, F_R_0;
    flux(W_L_plus, U_L_plus, F_L_plus);
    flux(W_R_plus, U_R_plus, F_R_plus);
    flux(W_L_0, U_L_0, F_L_0);
    flux(W_R_0, U_R_0, F_R_0);

    Conserved U_L_bar = U_L_plus+0.5*dt/dx*(F_L_plus-F_R_plus);
    Conserved U_R_bar = U_R_0+0.5*dt/dx*(F_L_0-F_R_0);

    Primitive W_L_bar, W_R_bar; 
    ConservedToPrimitive(U_L_bar, W_L_bar, gamma, p_Inf);
    ConservedToPrimitive(U_R_bar, W_R_bar, gamma, p_Inf);
    force(W_R_bar, W_L_bar, U_R_bar, U_L_bar, dt, dx, gamma, 
            p_Inf, f);
}

void godunov(Primitive W_L, Primitive W_R, double gamma, 
        double p_Inf, Conserved &f)
{
    Primitive W_RP;
    Conserved U_RP;
    double a_L = sqrt(gamma*(W_L.p+p_Inf)/W_L.rho);
    double a_R = sqrt(gamma*(W_R.p+p_Inf)/W_R.rho);
    double p_S, u_S;
    starRegionPressureVelocity(p_S, u_S, 
        W_L.rho, W_L.u, W_L.p, a_L, gamma, p_Inf, 
        W_R.rho, W_R.u, W_R.p, a_R, gamma, p_Inf);
    double e; 
    sample(W_RP.rho, W_RP.u, W_RP.p, e, p_S, u_S, 0, 
        W_L.rho, W_L.u, W_L.p, a_L, gamma, p_Inf, 
        W_R.rho, W_R.u, W_R.p, a_R, gamma, p_Inf);
    PrimitiveToConserved(W_RP, U_RP, gamma, p_Inf);
    flux(W_RP, U_RP, f);
}

void limitedSlopes(Conserved U_L, Conserved U_0, Conserved U_R, 
        char *limitFunc, Conserved &Delta)
{
    double b;
    if(!strcmp(limitFunc, "superbee")) b=2; 
    else b=1;

    Conserved Delta_L = U_0-U_L; 
    Conserved Delta_R = U_R-U_0;

    if(Delta_R.rho > 0)
    {
       double min1 = b*Delta_L.rho < Delta_R.rho ? 
           b*Delta_L.rho : Delta_R.rho;  
       double min2 = Delta_L.rho < b*Delta_R.rho ? 
           Delta_L.rho : b*Delta_R.rho;
       double max = min1 > min2 ? min1 : min2; 
       Delta.rho = max > 0 ? max : 0;
    }
    else
    {
       double max1 = b*Delta_L.rho > Delta_R.rho ? 
           b*Delta_L.rho : Delta_R.rho;  
       double max2 = Delta_L.rho > b*Delta_R.rho ? 
           Delta_L.rho : b*Delta_R.rho;
       double min = max1 < max2 ? max1 : max2; 
       Delta.rho = min < 0 ? min : 0;
    }

    if(Delta_R.rho_u > 0)
    {
       double min1 = b*Delta_L.rho_u < Delta_R.rho_u ? 
           b*Delta_L.rho_u : Delta_R.rho_u;  
       double min2 = Delta_L.rho_u < b*Delta_R.rho_u ? 
           Delta_L.rho_u : b*Delta_R.rho_u;
       double max = min1 > min2 ? min1 : min2; 
       Delta.rho_u = max > 0 ? max : 0;
    }
    else
    {
       double max1 = b*Delta_L.rho_u > Delta_R.rho_u ? 
           b*Delta_L.rho_u : Delta_R.rho_u;  
       double max2 = Delta_L.rho_u > b*Delta_R.rho_u ? 
           Delta_L.rho_u : b*Delta_R.rho_u;
       double min = max1 < max2 ? max1 : max2; 
       Delta.rho_u = min < 0 ? min : 0;
    }

    if(Delta_R.E > 0)
    {
       double min1 = b*Delta_L.E < Delta_R.E ? 
           b*Delta_L.E : Delta_R.E;  
       double min2 = Delta_L.E < b*Delta_R.E ? 
           Delta_L.E : b*Delta_R.E;
       double max = min1 > min2 ? min1 : min2; 
       Delta.E = max > 0 ? max : 0;
    }
    else
    {
       double max1 = b*Delta_L.E > Delta_R.E ? 
           b*Delta_L.E : Delta_R.E;  
       double max2 = Delta_L.E > b*Delta_R.E ? 
           Delta_L.E : b*Delta_R.E;
       double min = max1 < max2 ? max1 : max2; 
       Delta.E = min < 0 ? min : 0;
    }
}

void hllc(Primitive W_L, Primitive W_R, Conserved U_L, 
        Conserved U_R, double gamma, double p_Inf, Conserved &f)
{
    // Approximate HLLC solver
    
    // Sound speeds
    double a_L = sqrt(gamma*(W_L.p+p_Inf)/W_L.rho);
    double a_R = sqrt(gamma*(W_R.p+p_Inf)/W_R.rho);

    // Approximate pressure in star region
    double p_pvrs = 0.5*(W_L.p+W_R.p)
        -0.125*(W_R.u-W_L.u)*(W_L.rho+W_R.rho)*(a_L+a_R);

    // Set pressure to zero if approximation is negative
    double p_S = p_pvrs > 0 ? p_pvrs : 0;

    // Wave speed estimates
    
    // Left
    double q_L = 1;
    if(p_S > W_L.p)
        q_L = sqrt(1+(gamma+1)/gamma
                *((p_S+p_Inf)/(W_L.p+p_Inf)-1));
    double S_L = W_L.u-a_L*q_L;

    // Right
    double q_R = 1;
    if(p_S > W_R.p)
        q_R = sqrt(1+(gamma+1)/gamma
                *((p_S+p_Inf)/(W_R.p+p_Inf)-1));
    double S_R = W_R.u+a_R*q_R;

    // Star
    double S_S = (W_R.p-W_L.p+W_L.rho*W_L.u*(S_L-W_L.u)
            -W_R.rho*W_R.u*(S_R-W_R.u))/(W_L.rho*(S_L-W_L.u)
            -W_R.rho*(S_R-W_R.u));

    // Calculate star states and hllc flux according to speeds
    
    if(S_S>=0)
    { // Left of contact discontinuity
        Conserved F_L; 
        flux(W_L, U_L, F_L);
        if(S_L >=0)
        { // Left state
            f = F_L; 
        }
        else
        { // Left star state
            Conserved U_star_L(1, S_S, U_L.E/U_L.rho+(S_S-W_L.u)
                  *(S_S+W_L.p/(W_L.rho*(S_L-W_L.u))));
            U_star_L = W_L.rho*(S_L-W_L.u)/(S_L-S_S)*U_star_L; 
            f = F_L + S_L*(U_star_L-U_L);
/*          Second version:   
 *            
 *          Conserved D_S(0,1,S_S);
            f = (S_S*(S_L*U_L-F_L)+S_L*(W_L.p+W_L.rho*(S_L-W_L.u)
                        *(S_S-W_L.u))*D_S)/(S_L-S_S);*/
        }
    }
    else
    { // Right of contact discontinuity
        Conserved F_R; 
        flux(W_R, U_R, F_R);
        if(S_R <=0)
        { // Right state
            f = F_R; 
        }
        else
        { // Right star state
            Conserved U_star_R(1, S_S, U_R.E/U_R.rho+(S_S-W_R.u)
                  *(S_S+W_R.p/(W_R.rho*(S_R-W_R.u))));
            U_star_R = W_R.rho*(S_R-W_R.u)/(S_R-S_S)*U_star_R; 
            f = F_R + S_R*(U_star_R-U_R);
/*          Second version: 
 *
            Conserved D_S(0,1,S_S);
            f = (S_S*(S_R*U_R-F_R)+S_R*(W_R.p+W_R.rho*(S_R-W_R.u)
                        *(S_S-W_R.u))*D_S)/(S_R-S_S);*/
        }
    }
}

void muscl2(Conserved U_L, Conserved U_0, Conserved U_R, 
        Conserved U_2R, double dt, double dx, char *limitFunc, 
        double gamma, double p_Inf, Conserved &f)
{
    Conserved Delta_i, Delta_i_plus; 
    limitedSlopes(U_L, U_0, U_R, limitFunc, Delta_i);
    limitedSlopes(U_0, U_R, U_2R, limitFunc, Delta_i_plus);
/*  Slope limiters instead of limited slopes:
 *
    double xi_i = slopeLimiter(U_L.rho, U_0.rho, U_R.rho, 
            limitFunc);
    double xi_i_plus = slopeLimiter(U_0.rho, U_R.rho, U_2R.rho, 
            limitFunc);
    Delta_i = xi_i*0.5*(U_R-U_L);
    Delta_i_plus = xi_i_plus*0.5*(U_2R-U_0);
    */

    Conserved U_L_i, U_R_i, U_L_i_plus, U_R_i_plus;
    U_L_i = U_0-0.5*Delta_i;
    U_R_i = U_0+0.5*Delta_i;
    U_L_i_plus = U_R-0.5*Delta_i_plus;
    U_R_i_plus = U_R+0.5*Delta_i_plus;

    Primitive W_L_i, W_R_i, W_L_i_plus, W_R_i_plus;
    ConservedToPrimitive(U_L_i, W_L_i, gamma, p_Inf);
    ConservedToPrimitive(U_R_i, W_R_i, gamma, p_Inf);
    ConservedToPrimitive(U_L_i_plus, W_L_i_plus, gamma, p_Inf);
    ConservedToPrimitive(U_R_i_plus, W_R_i_plus, gamma, p_Inf);

    Conserved F_L_i, F_R_i, F_L_i_plus, F_R_i_plus;
    flux(W_L_i, U_L_i, F_L_i);
    flux(W_R_i, U_R_i, F_R_i);
    flux(W_L_i_plus, U_L_i_plus, F_L_i_plus);
    flux(W_R_i_plus, U_R_i_plus, F_R_i_plus);

    Conserved U_bar_L, U_bar_R; 
    U_bar_L = U_L_i_plus+0.5*dt/dx*(F_L_i_plus-F_R_i_plus);
    U_bar_R = U_R_i+0.5*dt/dx*(F_L_i-F_R_i);

    Primitive W_bar_L, W_bar_R; 
    ConservedToPrimitive(U_bar_L, W_bar_L, gamma, p_Inf);
    ConservedToPrimitive(U_bar_R, W_bar_R, gamma, p_Inf);

    hllc(W_bar_R, W_bar_L, U_bar_R, U_bar_L, gamma, p_Inf, f);
}

void advance(Primitive *W, Conserved *U, Conserved *U_old, 
        double dt, double dx, int N, double gamma, double p_Inf, 
        char *scheme, char *limitFunc)
{
    Conserved *f = new Conserved[N+1];
    int L, R;

    for(int i=0; i<N+1; i++)
    {
        L = i+NGC-1; R = i+NGC;
        if(!strcmp(scheme, "FORCE"))
            force(W[L], W[R], U_old[L], U_old[R], dt, dx, 
                    gamma, p_Inf, f[i]);
        else if(!strcmp(scheme, "FLIC"))
            flic(W[L], W[R], U_old[L], U_old[R], U_old[L-1], 
                    U_old[R+1], dt, dx, limitFunc, gamma, 
                    p_Inf, f[i]);
        else if(!strcmp(scheme, "SLIC"))
            slic2(U_old[L-1], U_old[L], U_old[R], U_old[R+1],
                    dt, dx, limitFunc, gamma, p_Inf, f[i]);
        else if(!strcmp(scheme, "MUSCL"))
            muscl2(U_old[L-1], U_old[L], U_old[R], 
                    U_old[R+1], dt, dx, limitFunc, gamma, 
                    p_Inf, f[i]);
    }
    for(int i=NGC; i<N+NGC; i++)
    {
        L = i-NGC; R = i-NGC+1;
        U[i] = U_old[i]-dt/dx*(f[R]-f[L]); 
        ConservedToPrimitive(U[i], W[i], gamma, p_Inf); 
    }
    for(int i=NGC; i<N+NGC; i++)
    {
        U_old[i] = U[i];
    }
    delete []f; f = NULL;
}

void initiateTestCase(int testCase, int &nInterfaces, 
        double x[3], double rho[4], double u[4], double p[4], 
        int mat[4], double gamma[2], double p_Inf[2], 
        double &tStop){
    gamma[0]=1.4;
    mat[0]=0;
    mat[1]=1;
    p_Inf[0]=0.0;
    p_Inf[1]=0.0;

    if(testCase<8){
        nInterfaces=1;
        if(testCase==7) gamma[1]=1.67;
        else gamma[1]=1.4;
        if(testCase<6) x[0]=0.5;
        else x[0]=0.25;
    }

    if(testCase==1){
        rho[0] = 1.0;
        rho[1] = 0.125;
        u[0] = 0.0;
        u[1] = 0.0;
        p[0] = 1.0;
        p[1] = 0.1;
        tStop = 0.25; 
    }else if(testCase==2){
        rho[0] = 1.0;
        rho[1] = 1.0;
        u[0] = -2.0;
        u[1] = 2.0;
        p[0] = 0.4;
        p[1] = 0.4;
        tStop = 0.15; 
    }else if(testCase==3){
        rho[0] = 1.0;
        rho[1] = 1.0;
        u[0] = 0.0;
        u[1] = 0.0;
        p[0] = 1000;
        p[1] = 0.01;
        tStop = 0.012; 
    }else if(testCase==4){
        rho[0] = 1.0;
        rho[1] = 1.0;
        u[0] = 0.0;
        u[1] = 0.0;
        p[0] = 0.01;
        p[1] = 100.0;
        tStop = 0.035; 
    }else if(testCase==5){
        rho[0] = 5.99924;
        rho[1] = 5.9942;
        u[0] = 19.5975;
        u[1] = -6.19633;
        p[0] = 460.894;
        p[1] = 46.0950;
        tStop = 0.035; 
    }else if(testCase<8){
        rho[0] = 1.0;
        rho[1] = 0.138;
        u[0] = 0.5;
        u[1] = 0.5;
        p[0] = 1.0;
        p[1] = 1.0;
        tStop = 1.0; 
    }else if(testCase==8){
        nInterfaces = 2;
        x[0] = 0.05;
        x[1] = 0.5;
        rho[0] = 1.3333;
        rho[1] = 1.0;
        rho[2] = 0.1379;
        u[0] = 0.3535*sqrt(100000);
        u[1] = 0;
        u[2] = 0;
        p[0] = 1.5e5;
        p[1] = 1.0e5;
        p[2] = 1.0e5;
        mat[0] = 0;
        mat[1] = 0;
        mat[2] = 1;
        gamma[1] = 1.67;
        tStop = 0.0012;
    }else if(testCase==9){
        nInterfaces = 3;
        x[0] = 0.05;
        x[1] = 0.4;
        x[2] = 0.6;
        rho[0] = 1.3333;
        rho[1] = 1.0;
        rho[2] = 0.1379;
        rho[3] = 1.0;
        u[0] = 0.3535*sqrt(100000);
        u[1] = 0;
        u[2] = 0;
        u[3] = 0;
        p[0] = 1.5e5;
        p[1] = 1.0e5;
        p[2] = 1.0e5;
        p[3] = 1.0e5;
        mat[0] = 0;
        mat[1] = 0;
        mat[2] = 1;
        mat[3] = 0;
        gamma[1] = 1.67;
        tStop = 0.0014;
    }else if(testCase==10){
        nInterfaces = 1; 
        x[0] = 0.7;
        rho[0] = 1000;
        rho[1] = 50;
        u[0] = 0;
        u[1] = 0;
        p[0] = 1.0e9;
        p[1] = 1.0e5;
        mat[0] = 0;
        mat[1] = 1;
        gamma[0] = 4.4;
        gamma[1] = 1.4;
        p_Inf[0] = 6.0e8;
        p_Inf[1] = 0;
        tStop = 237.44e-6; 
    }
}

void setScheme(char *scheme, char *limitFunc, 
        int schemeChoice, int limitChoice)
{
    if(schemeChoice == 0)
        strcpy(scheme, "FORCE");
    else if(schemeChoice == 1)
        strcpy(scheme, "FLIC");
    else if(schemeChoice == 2)
        strcpy(scheme, "SLIC"); 
    else if(schemeChoice == 3)
        strcpy(scheme, "MUSCL");
    if(limitChoice == 0)
        strcpy(limitFunc, "minbee");
    else if(limitChoice == 1)
        strcpy(limitFunc, "vanleer");
    else if(limitChoice == 2)
        strcpy(limitFunc, "superbee");
}

void starRegionPressureVelocity(double &p, double &u,
        double rho_L, double u_L, double p_L, double a_L, 
        double gamma_L, double p_Inf_L, double rho_R, double u_R, 
        double p_R, double a_R, double gamma_R, double p_Inf_R)
{
    // Compute the pressure and velocity in the star region.
    // 
    // First, the pressure in the star region is computed to 
    // an accuracy of tol by Newton-Raphson root finding, 
    // with initial guess p_init = 0.5(p_left+p_right). 
    // The procedure cancels if the amount of iterations 
    // exceeds maxIter. 
    // 
    // Afterwards, the velocity u is computed. 
    int iter = 0, maxIter = 20;
    double p_Old =0.5*(p_L+p_R),
           change = 1, 
           tol = 1e-6,
           f_L, f_L_Diff, f_R, f_R_Diff;

    // Calculation of pressure in star region
    while(change > tol)
    {
        pressureFunctions(f_L, f_L_Diff, p_Old, 
                rho_L, p_L, a_L, gamma_L, p_Inf_L);
        pressureFunctions(f_R, f_R_Diff, p_Old, 
                rho_R, p_R, a_R, gamma_R, p_Inf_R);
        p = p_Old-(f_L+f_R+u_R-u_L)/(f_L_Diff+f_R_Diff);
        if(p < 0) p = tol;
        change = 2*fabs((p-p_Old)/(p+p_Old));
        p_Old = p;
        iter++;
        if(iter > maxIter)
        {
            qDebug() << "Pressure in star region N/A, iter > " 
                << maxIter;
            return;
        }
    }
    u = 0.5*(u_L+u_R+f_R-f_L);
}

void pressureFunctions(double &f, double &f_Diff, double p_Old, 
        double rho, double p, double a, double g, 
        double p_Inf)
{
    // Calculate the pressure functions f_L and f_R as given in 
    // Toro, and their first derivatives wrt pressure. 
    if(p_Old <= p) 
    { 
        // Rarefaction
        f = 2*a/(g-1)
            *(pow((p_Old+p_Inf)/(p+p_Inf), 0.5*(g-1)/g)-1);
        f_Diff = 1/(rho*a)
            *pow((p_Old+p_Inf)/(p+p_Inf), -0.5*(g+1)/g);
    }
    else 
    {
        // Shock
        double A = 2/(rho*(g+1));
        double B = (g-1)/(g+1)*p+2*g*p_Inf/(g+1);
        f = (p_Old-p)*sqrt(A/(B+p_Old));
        f_Diff = sqrt(A/(p_Old+B))*(1-0.5*(p_Old-p)/(p_Old+B));
    }
}

void starRegionRho(double &rho_L_S, double &rho_R_S, double p_S, 
        double rho_L, double p_L, double g_L, double p_Inf_L,
        double rho_R, double p_R, double g_R, double p_Inf_R)
{
    // Find density in star region
    if(p_S <= p_L)   
    { 
        // Left rarefaction
        rho_L_S = rho_L*pow((p_S+p_Inf_L)/(p_L+p_Inf_L), 1.0/g_L);
    }
    else 
    { 
        // Left shock
        rho_L_S = rho_L*((p_S+p_Inf_L)/(p_L+p_Inf_L)
                +(g_L-1)/(g_L+1))
            /((p_S+p_Inf_L)/(p_L+p_Inf_L)*(g_L-1)/(g_L+1)+1);
    }
    if(p_S > p_R)
    { 
        // Right shock
       rho_R_S = rho_R*((p_S+p_Inf_R)/(p_R+p_Inf_R)
               +(g_R-1)/(g_R+1))
            /((p_S+p_Inf_R)/(p_R+p_Inf_R)*(g_R-1)/(g_R+1)+1);
    }
    else
    { 
        // Right rarefaction 
        rho_R_S = rho_R*pow((p_S+p_Inf_R)/(p_R+p_Inf_R), 1.0/g_R);
    }
}

void sample(double &rho, double &u, double &p, double &e,
        double p_S, double u_S, double S,
        double rho_L, double u_L, double p_L, double a_L, 
        double g_L, double p_Inf_L, double rho_R, double u_R, 
        double p_R, double a_R, double g_R, double p_Inf_R)
{
    // Sample the solution of the Riemann problem. See Toro. 
    double S_H_L, a_S_L, S_T_L, S_H_R, a_S_R, S_T_R, a; 
    if(S <= u_S)
    { 
        // Sampling point to left of contact discontinuity
        if(p_S <= p_L)   
        { 
            // Left rarefaction
            // Speed of head:
            S_H_L = u_L - a_L; 
            if(S <= S_H_L) 
            { 
                // Sampling point is left data state
                rho = rho_L;
                u = u_L; 
                p = p_L;
            }
            else
            {   
                // Sound speed behind rarefaction:
                a_S_L = a_L*pow((p_S+p_Inf_L)/(p_L+p_Inf_L), 
                        0.5*(g_L-1)/g_L); 
                // Speed of tail:
                S_T_L = u_S - a_S_L; 
                if(S > S_T_L) 
                { 
                    // Sampling point is left star state
                    rho = rho_L*pow((p_S+p_Inf_L)/(p_L+p_Inf_L), 
                            1.0/g_L);
                    u = u_S;
                    p = p_S;
                }
                else 
                { 
                    // Sampling point is inside left fan
                    u = 2/(g_L+1)*(a_L+0.5*(g_L-1)*u_L+S);
                    a = 2/(g_L+1)*(a_L+0.5*(g_L-1)*(u_L-S));
                    rho = rho_L*pow(a/a_L, 2/(g_L-1));
                    p = (p_L+p_Inf_L)*(pow(a/a_L, 
                                2*g_L/(g_L-1)))-p_Inf_L;
                }
            }
        }
        else 
        { 
            // Left shock
            if(S <= u_L-a_L*sqrt(0.5*(g_L+1)/g_L
                        *(p_S+p_Inf_L)/(p_L+p_Inf_L)
                        +0.5*(g_L-1)/g_L))
            {
                // Left data state
                rho = rho_L;
                u = u_L;
                p = p_L;
            }
            else
            { 
                // Sampling point is left star state
                rho = rho_L*((p_S+p_Inf_L)/(p_L+p_Inf_L)
                        +(g_L-1)/(g_L+1))/((p_S+p_Inf_L)
                        /(p_L+p_Inf_L)*(g_L-1)/(g_L+1)+1);
                u = u_S;
                p = p_S;
            }
        }
        e = (p+g_L*p_Inf_L)/(rho*(g_L-1));
    }
    else
    { 
        // Sampling point is to the right of contact discontinuity
        if(p_S > p_R)
        { 
            // Right shock
            if(S > u_R+a_R*sqrt(0.5*(g_R+1)/g_R
                        *(p_S+p_Inf_R)/(p_R+p_Inf_R)
                        +0.5*(g_R-1)/g_R)) 
            {
                //Right data state
                rho = rho_R;
                u = u_R;
                p = p_R;
            }
            else
            { 
                // Sampled point is right star state
                rho = rho_R*((p_S+p_Inf_R)/(p_R+p_Inf_R)
                        +(g_R-1)/(g_R+1))
                    /((p_S+p_Inf_R)/(p_R+p_Inf_R)
                            *(g_R-1)/(g_R+1)+1);
                u = u_S;
                p = p_S;
            }
        }
        else
        { 
            // Right rarefaction 
            S_H_R = u_R + a_R; 
            // Speed of head
            if(S >= S_H_R)
            {
                // Sampling point is right data state
                rho = rho_R;
                u = u_R;
                p = p_R;
            }
            else
            {  
                // Sound speed behind rarefaction:
                a_S_R = a_R*pow((p_S+p_Inf_R)/(p_R+p_Inf_R), 
                        0.5*(g_R-1)/g_R); 
                // Speed of tail:
                S_T_R = u_S + a_S_R; 
                if(S <= S_T_R)
                { 
                    // Sampling point is right star state
                    rho = rho_R*pow((p_S+p_Inf_R)/(p_R+p_Inf_R), 
                            1.0/g_R);
                    u = u_S;
                    p = p_S;
                }
                else
                { 
                    // Sampling point is inside right fan
                    u = 2/(g_R+1)*(-a_R+0.5*(g_R-1)*u_R+S);
                    a = 2/(g_R+1)*(a_R-0.5*(g_R-1)*(u_R-S));
                    rho = rho_R*pow(a/a_R, 2/(g_R-1));
                    p = (p_R+p_Inf_R)*(pow(a/a_R, 2*g_R/(g_R-1)))
                        -p_Inf_R;
                }
            }
        }
        e = (p+g_R*p_Inf_R)/(rho*(g_R-1));
    }
}


#endif

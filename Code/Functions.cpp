#ifndef __FUNCTIONS_CPP
#define __FUNCTIONS_CPP

#include "Functions.h"

void initialConditions(Primitive *W, Conserved *U, 
        double x_0, int N, double dx,
        double rho_L, double u_L, double p_L, 
        double rho_R, double u_R, double p_R){
   for(int i=0; i<N+2*NGC; i++){
       if((i-NGC+0.5)*dx < x_0){
           W[i].rho = rho_L; 
           W[i].u = u_L; 
           W[i].p = p_L; 
       }else{
           W[i].rho = rho_R; 
           W[i].u = u_R; 
           W[i].p = p_R; 
       }
       PrimitiveToConserved(W[i], U[i]);
   } 
}

double maxWaveSpeed(Primitive *W, double a, int N){
    // Smax = max_i {abs(S_L_i+0.5), abs(S_R_i+0.5)}
    double u_max = 0;
    for(int i=0; i<N+2*NGC; i++){
        if(u_max < fabs(W[i].u)) u_max = fabs(W[i].u);
    }
    return u_max+a; 
}

void PrimitiveToConserved(Primitive W, Conserved &U){
    U.rho = W.rho;
    U.rho_u = W.rho*W.u; 
    U.E = 0.5*W.rho*W.u*W.u + W.p/(g-1);
}

void ConservedToPrimitive(Conserved U, Primitive &W){
    W.rho = U.rho;
    W.u = U.rho_u/U.rho; 
    W.p = (U.E-0.5*U.rho_u*U.rho_u/U.rho)*(g-1);
}

void boundaryConditions(Primitive *W, Conserved *U, int N){
    // Transmissive boundary conditions
    for(int i=0; i<NGC; i++){
        W[i] = W[i+1];
        W[N+NGC-i] = W[N+NGC-i-1];
        PrimitiveToConserved(W[i], U[i]);
        PrimitiveToConserved(W[N+NGC-i], U[N+NGC-i]); 
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
        Conserved &f){
    Conserved F_L, F_R, U_half;
    Primitive W_half;
    flux(W_L, U_L, F_L); 
    flux(W_R, U_R, F_R); 
    U_half = 0.5*((U_L+U_R)
            +dt/dx*(F_L-F_R));
    ConservedToPrimitive(U_half, W_half);
    flux(W_half, U_half, f);
}

void force(Primitive W_L, Primitive W_R, 
        Conserved U_L, Conserved U_R, double dt, double dx,
        Conserved &f){
    Conserved f_LF, f_RI;
    laxFriedrich(W_L, W_R, U_L, U_R, dt,dx, f_LF);
    richtmeyer(W_L, W_R, U_L, U_R, dt, dx, f_RI);
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
    }else{ //superbee
        phi_L = superbee(r_L, 0.9);
        phi_R = superbee(r_R, 0.9);
    }
    return (phi_L < phi_R) ? phi_L : phi_R;
}

void flic(Primitive W_L, Primitive W_R, 
        Conserved U_L, Conserved U_R, 
        Conserved U_2L, Conserved U_2R, 
        double dt, double dx, char *limitFunc, 
        Conserved &f){
    Conserved f_RI, f_FORCE;
    double phi = fluxLimiter(U_2L.rho, U_L.rho, 
            U_R.rho, U_2R.rho, limitFunc);
    richtmeyer(W_L, W_R, U_L, U_R, dt, dx, f_RI);
    force(W_L, W_R, U_L, U_R, dt, dx, f_FORCE);
    f = f_FORCE+phi*(f_RI-f_FORCE);
}

void hllc(Primitive W_L, Primitive W_R, Conserved &f){
    Conserved U_L, U_R, F_L, F_R;
    double S_L, S_R, S_plus, S_star;
    S_L = fabs(W_L.u)+sqrt(g*W_L.p/W_L.rho);
    S_R = fabs(W_R.u)+sqrt(g*W_R.p/W_R.rho);
    S_plus = S_L > S_R ? S_L : S_R;
    S_L = -S_plus;
    S_R = S_plus;
    PrimitiveToConserved(W_L, U_L);
    flux(W_L, U_L, F_L);
    PrimitiveToConserved(W_R, U_R);
    flux(W_R, U_R, F_R);
    S_star = (W_R.p-W_L.p+W_L.rho*W_L.u*(S_L-W_L.u)-W_R.rho*W_R.u*(S_R-W_R.u))/(W_L.rho*(S_L-W_L.u)-W_R.rho*(S_R-W_R.u));
    if(S_star >= 0){
    Conserved U_L_star(1, S_star, 
            U_L.E/U_L.rho+(S_star-W_L.u)
            *(S_star+W_L.p/(W_L.rho*(S_L-W_L.u))));
    U_L_star = U_L_star*W_L.rho*(S_L-W_L.u)/(S_L-S_star);
    f = F_L+S_L*(U_L_star-U_L);
    }else{
    Conserved U_R_star(1, S_star, 
            U_R.E/U_R.rho+(S_star-W_R.u)
            *(S_star+W_R.p/(W_R.rho*(S_R-W_R.u))));
    U_R_star = U_R_star*W_R.rho*(S_R-W_R.u)/(S_R-S_star);
    f = F_R+S_R*(U_R_star-U_R);
    }
}

void godunov(Primitive W_L, Primitive W_R, char *RP, Conserved &f){
    Primitive W_RP;
    Conserved U_RP;
    double p_S, u_S;
    if(!strcmp(RP, "Exact")){
        starRegionPressureVelocity(p_S, u_S, 
            W_L.rho, W_L.u, W_L.p, sqrt(g*W_L.p/W_L.rho), 
            W_R.rho, W_R.u, W_R.p, sqrt(g*W_R.p/W_R.rho));
        sample(W_RP.rho, W_RP.u, W_RP.p, p_S, u_S, 0, 
            W_L.rho, W_L.u, W_L.p, sqrt(g*W_L.p/W_L.rho), 
            W_R.rho, W_R.u, W_R.p, sqrt(g*W_R.p/W_R.rho));
        PrimitiveToConserved(W_RP, U_RP);
        flux(W_RP, U_RP, f);
    }else if(!strcmp(RP, "HLLC")){
        hllc(W_L, W_R, f);
    }
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

void muscl(Conserved *U, 
        double dt, double dx, double omega, 
        int N, char *limitFunc, char *RP, 
        Conserved *f){
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
        ConservedToPrimitive(U_L[i], W_L[i]);
        ConservedToPrimitive(U_R[i], W_R[i]);
        flux(W_L[i], U_L[i], F_L[i]);
            flux(W_R[i], U_R[i], F_R[i]);
        U_L_bar[i] = U_L[i]+0.5*dt/dx*(F_L[i]-F_R[i]);
        U_R_bar[i] = U_R[i]+0.5*dt/dx*(F_L[i]-F_R[i]);
        ConservedToPrimitive(U_L_bar[i], W_L_bar[i]);
        ConservedToPrimitive(U_R_bar[i], W_R_bar[i]);
    }
    for(int i=0; i<N+1; i++){
        L = i+NGC-1; R = i+NGC;
        godunov(W_R_bar[L], W_L_bar[R], RP, f[i]);
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

void slic(Conserved *U, 
        double dt, double dx, double omega, 
        int N, char *limitFunc, Conserved *f){
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
        ConservedToPrimitive(U_L[i], W_L[i]);
        ConservedToPrimitive(U_R[i], W_R[i]);
        flux(W_L[i], U_L[i], F_L[i]);
            flux(W_R[i], U_R[i], F_R[i]);
        U_L_bar[i] = U_L[i]+0.5*dt/dx*(F_L[i]-F_R[i]);
        U_R_bar[i] = U_R[i]+0.5*dt/dx*(F_L[i]-F_R[i]);
        ConservedToPrimitive(U_L_bar[i], W_L_bar[i]);
        ConservedToPrimitive(U_R_bar[i], W_R_bar[i]);
    }
    for(int i=0; i<N+1; i++){
        L = i+NGC-1; R = i+NGC;
        force(W_R_bar[L], W_L_bar[R], U_R_bar[L], U_L_bar[R], 
                dt, dx, f[i]); 
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
        double dt, double dx, int N, 
        char *scheme, char *limitFunc, char *RP){
    Conserved *f = new Conserved[N+1];
    int L, R;
    if(!strcmp(scheme, "MUSCL"))
        muscl(U_old, dt, dx, 0.0, N, limitFunc, RP, f);
    else if(!strcmp(scheme, "SLIC"))
        slic(U_old, dt, dx, 0.0, N, limitFunc, f); 
    else{
        for(int i=0; i<N+1; i++){
            L = i+NGC-1; R = i+NGC;
            if(!strcmp(scheme, "FORCE"))
                force(W[L], W[R], U_old[L], U_old[R], 
                        dt, dx, f[i]);
            else if(!strcmp(scheme, "FLIC"))
                flic(W[L], W[R], U_old[L], U_old[R], U_old[L-1], 
                        U_old[R+1], dt, dx, limitFunc, f[i]);
            else if(!strcmp(scheme, "Godunov"))
                godunov(W[L], W[R], RP, f[i]);
        }
    }
    for(int i=NGC; i<N+NGC; i++){
        L = i-NGC; R = i-NGC+1;
        U[i] = U_old[i]-dt/dx*(f[R]-f[L]); 
        ConservedToPrimitive(U[i], W[i]); 
    }
    for(int i=NGC; i<N+NGC; i++){
        U_old[i] = U[i];
    }
    delete []f; f = NULL;
}

void initiateTestCase(int testCase, double &tStop, double &x_0,
        double &rho_L, double &u_L, double &p_L, 
        double &rho_R, double &u_R, double &p_R){
    if(testCase==1){
        rho_L = 1.0;
        u_L = 0.75;
        p_L = 1.0;
        rho_R = 0.125;
        u_R = 0.0;
        p_R = 0.1;
        x_0 = 0.3;
        tStop = 0.2; 
    }else if(testCase==2){
        rho_L = 1.0;
        u_L = -2.0;
        p_L = 0.4;
        rho_R = 1.0;
        u_R = 2.0;
        p_R = 0.4;
        x_0 = 0.5;
        tStop = 0.15; 
    }else if(testCase==3){
        rho_L = 1.0;
        u_L = 0.0;
        p_L = 1000;
        rho_R = 1.0;
        u_R = 0.0;
        p_R = 0.01;
        x_0 = 0.5;
        tStop = 0.012; 
    }else if(testCase==4){
        rho_L = 5.99924;
        u_L = 19.5975;
        p_L = 460.894;
        rho_R = 5.9942;
        u_R = -6.19633;
        p_R = 46.0950;
        x_0 = 0.4;
        tStop = 0.035; 
    }else if(testCase==5){
        rho_L = 1.0;
        u_L = -19.59745;
        p_L = 1000;
        rho_R = 1.0;
        u_R = -19.59745;
        p_R = 0.01;
        x_0 = 0.8;
        tStop = 0.012; 
    }
}

void starRegionPressureVelocity(double &p, double &u,
        double rho_L, double u_L, double p_L, double a_L,
        double rho_R, double u_R, double p_R, double a_R){
    // Compute the pressure and velocity in the star region.
    // 
    // First, the pressure in the star region is computed to 
    // an accuracy of tol by Newton-Raphson root finding, 
    // with initial guess p_init = 0.5(p_left+p_right). 
    // The procedure cancels if the amount of iterations 
    // exceeds maxIter. 
    // 
    // Afterwards, the velocity u is computed. 
    int iter = 0, maxIter = 100;
    double p_Old =0.5*(p_L+p_R),
           change = 1, 
           tol = 1e-6,
           f_L, f_L_Diff, f_R, f_R_Diff;

/*    printf("---------------------------------------\n");
    printf("Calculation of pressure in star region\n");
    printf("\ti\tp\tchange\n");
    printf("---------------------------------------\n"); */

    while(change > tol){
        pressureFunctions(f_L, f_L_Diff, p_Old, 
                rho_L, u_L, p_L, a_L);
        pressureFunctions(f_R, f_R_Diff, p_Old, 
                rho_R, u_R, p_R, a_R);
        p = p_Old-(f_L+f_R+u_R-u_L)/(f_L_Diff+f_R_Diff);
        if(p < 0) p = tol;
        change = 2*fabs((p-p_Old)/(p+p_Old));
        p_Old = p;
        iter++;
//        printf("\t%d\t%f\t%f\n", iter, p, change);
        if(iter > maxIter){
            printf("---------------------------------------\n");
            printf("Pressure in star region N/A, iter > %d.\n", 
                    maxIter);
            return;
        }
    }
    u = 0.5*(u_L+u_R+f_R-f_L);
}

void pressureFunctions(double &f, double &f_Diff, double p_Old, 
        double rho, double u, double p, double a){
    // Calculate the pressure functions f_L and f_R as given in 
    // Toro, and their first derivatives wrt pressure. 
    if(p_Old <= p){ // Rarefaction
        f = 2*a/(g-1)*(pow(p_Old/p, 0.5*(g-1)/g)-1);
        f_Diff = 1/(rho*a)*pow(p_Old/p, -0.5*(g+1)/g);
    }else{ // Shock
        double A = 2/(rho*(g+1));
        double B = (g-1)/(g+1)*p;
        f = (p_Old-p)*sqrt(A/(B+p_Old));
        f_Diff = sqrt(A/(p_Old+B))*(1-0.5*(p_Old-p)/(p_Old+B));
    }
}


void sample(double &rho, double &u, double &p,
        double p_S, double u_S, double S,
        double rho_L, double u_L, double p_L, double a_L,
        double rho_R, double u_R, double p_R, double a_R){
    // Sample the solution of the Riemann problem. See Toro. 
    double S_H_L, a_S_L, S_T_L, S_H_R, a_S_R, S_T_R, a; 
    if(S <= u_S){ 
        // Sampling point is to the left of contact discontinuity
        if(p_S <= p_L){ 
            // Left rarefaction
            S_H_L = u_L - a_L; 
            // Speed of head
            if(S <= S_H_L){ 
                // Sampling point is left data state
                rho = rho_L;
                u = u_L; 
                p = p_L;
            }else{ 
                a_S_L = a_L*pow(p_S/p_L, 0.5*(g-1)/g); 
                // Sound speed behind rarefaction
                S_T_L = u_S - a_S_L; // Speed of tail 
                if(S > S_T_L){ 
                    // Sampling point is left star state
                    rho = rho_L*pow(p_S/p_L, 1.0/g);
                    u = u_S;
                    p = p_S;
                }else{ 
                    // Sampling point is inside left fan
                    u = 2/(g+1)*(a_L+0.5*(g-1)*u_L+S);
                    a = 2/(g+1)*(a_L+0.5*(g-1)*(u_L-S));
                    rho = rho_L*pow(a/a_L, 2/(g-1));
                    p = p_L*(pow(a/a_L, 2*g/(g-1)));
                }
            }
        }else{ 
            // Left shock
            if(S <= u_L-
                    a_L*sqrt(0.5*(g+1)/g*p_S/p_L+0.5*(g-1)/g)){ 
                // Left data state
                rho = rho_L;
                u = u_L;
                p = p_L;
            }else{ 
                // Sampling point is left star state
                rho = rho_L*(p_S/p_L+(g-1)/(g+1))
                    /(p_S/p_L*(g-1)/(g+1)+1);
                u = u_S;
                p = p_S;
            }
        }
    }else{ 
        // Sampling point is to the right of contact discontinuity
        if(p_S > p_R){ 
            // Right shock
            if(S > u_R
                    + a_R*sqrt(0.5*(g+1)/g*p_S/p_R+0.5*(g-1)/g)){ 
                //Right data state
                rho = rho_R;
                u = u_R;
                p = p_R;
            }else{ 
                // Sampled point is right star state
                rho = rho_R*(p_S/p_R+(g-1)/(g+1))
                    /(p_S/p_R*(g-1)/(g+1)+1);
                u = u_S;
                p = p_S;
            }
        }else{ 
            // Right rarefaction 
            S_H_R = u_R + a_R; 
            // Speed of head
            if(S >= S_H_R){ 
                // Sampling point is right data state
                rho = rho_R;
                u = u_R;
                p = p_R;
            }else{  
                a_S_R = a_R*pow(p_S/p_R, 0.5*(g-1)/g); 
                // Sound speed behind rarefaction
                S_T_R = u_S + a_S_R; // Speed of tail 
                if(S <= S_T_R){ 
                    // Sampling point is right star state
                    rho = rho_R*pow(p_S/p_R, 1.0/g);
                    u = u_S;
                    p = p_S;
                }else{ 
                    // Sampling point is inside right fan
                    u = 2/(g+1)*(-a_R+0.5*(g-1)*u_R+S);
                    a = 2/(g+1)*(a_R-0.5*(g-1)*(u_R-S));
                    rho = rho_R*pow(a/a_R, 2/(g-1));
                    p = p_R*(pow(a/a_R, 2*g/(g-1)));
                }
            }
        }
    }
}

#endif

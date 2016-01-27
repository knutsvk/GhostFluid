#ifndef __TYPES_H
#define __TYPES_H

// Primitive vector W: 
struct Primitive{
    double rho; 
    double u; 
    double p;

    Primitive(){}
    Primitive(double _rho, double _u, double _p);
    virtual ~Primitive(){}

    Primitive& operator=(const Primitive &rhs);
    Primitive operator+(const Primitive &rhs);
    Primitive operator-(const Primitive &rhs);
    Primitive operator*(const double &factor);
    Primitive operator/(const double &factor);
};
Primitive operator*(const double &factor, const Primitive &p);
Primitive operator/(const double &factor, const Primitive &p);

// Conserved vector U:
struct Conserved{
    double rho; 
    double rho_u;
    double E; 

    Conserved(){}
    Conserved(double _rho, double _rho_u, double _E);
    virtual ~Conserved(){}

    Conserved& operator=(const Conserved &rhs);
    Conserved operator+(const Conserved &rhs);
    Conserved operator-(const Conserved &rhs);
    Conserved operator*(const double &factor);
    Conserved operator/(const double &factor);
};
Conserved operator*(const double &factor, const Conserved &c);
Conserved operator/(const double &factor, const Conserved &c);

#endif

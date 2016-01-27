#ifndef __TYPES_CPP
#define __TYPES_CPP

#include "Types.h"

Primitive::Primitive(double _rho, double _u, double _p){
    rho = _rho;
    u = _u; 
    p = _p;
}

Primitive& Primitive::operator=(const Primitive &rhs){
    if(&rhs == this) return *this;
    rho = rhs.rho;
    u = rhs.u;
    p = rhs.p;
    return *this;
}

Primitive Primitive::operator+(const Primitive &rhs){
    Primitive result(*this);
    result.rho += rhs.rho;
    result.u += rhs.u;
    result.p += rhs.p;
    return result;
}

Primitive Primitive::operator-(const Primitive &rhs){
    Primitive result(*this);
    result.rho -= rhs.rho;
    result.u -= rhs.u;
    result.p -= rhs.p;
    return result;
}

Primitive Primitive::operator*(const double &factor){
    Primitive result(*this);
    result.rho *= factor;
    result.u *= factor;
    result.p *= factor;
    return result; 
}

Primitive Primitive::operator/(const double &factor){
    Primitive result(*this);
    result.rho /= factor;
    result.u /= factor;
    result.p /= factor;
    return result; 
}

Primitive operator*(const double &factor, const Primitive &p){
    // To allow multiplication with factor to the right
    Primitive result(p);
    result.rho *= factor;
    result.u *= factor;
    result.p *= factor;
    return result; 
}

Primitive operator/(const double &factor, const Primitive &p){
    // To allow multiplication with factor to the right
    Primitive result(p);
    result.rho /= factor;
    result.u /= factor;
    result.p /= factor;
    return result; 
}

Conserved::Conserved(double _rho, double _rho_u, double _E){
    rho = _rho;
    rho_u = _rho_u; 
    E = _E;
}

Conserved& Conserved::operator=(const Conserved &rhs){
    if(&rhs == this) return *this;
    rho = rhs.rho;
    rho_u = rhs.rho_u;
    E = rhs.E;
    return *this;
}

Conserved Conserved::operator+(const Conserved &rhs){
    Conserved result(*this);
    result.rho += rhs.rho;
    result.rho_u += rhs.rho_u;
    result.E += rhs.E;
    return result;
}

Conserved Conserved::operator-(const Conserved &rhs){
    Conserved result(*this);
    result.rho -= rhs.rho;
    result.rho_u -= rhs.rho_u;
    result.E -= rhs.E;
    return result;
}

Conserved Conserved::operator*(const double &factor){
    Conserved result(*this);
    result.rho *= factor;
    result.rho_u *= factor;
    result.E *= factor;
    return result; 
}

Conserved Conserved::operator/(const double &factor){
    Conserved result(*this);
    result.rho /= factor;
    result.rho_u /= factor;
    result.E /= factor;
    return result; 
}

Conserved operator*(const double &factor, const Conserved &c){
    // To allow multiplication with factor to the right
    Conserved result(c);
    result.rho *= factor;
    result.rho_u *= factor;
    result.E *= factor;
    return result; 
}

Conserved operator/(const double &factor, const Conserved &c){
    // To allow multiplication with factor to the right
    Conserved result(c);
    result.rho /= factor;
    result.rho_u /= factor;
    result.E /= factor;
    return result; 
}

#endif

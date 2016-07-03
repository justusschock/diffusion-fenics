//
// Created by js on 03.07.16.
//

#ifndef DIFFUSION_FENICS_PDEHELPER_H
#define DIFFUSION_FENICS_PDEHELPER_H

#include <dolfin.h>

class Initial : public dolfin::Expression {
    void eval(dolfin::Array<double> &values, const dolfin::Array<double> &x) const{
        values[0] = 1;
    }
};

//sourceTerm (right-hand side)
class Source : public dolfin::Expression {
    void eval(dolfin::Array<double> &values, const dolfin::Array<double> &x) const {
        values[0] = -1;
        if(x[0] > 1 - (0.5+DOLFIN_EPS) and x[0] < 1 - (0.5-DOLFIN_EPS) and x[1] > 1 - (0.5+DOLFIN_EPS) and x[1] < 1 - (0.5-DOLFIN_EPS))
            values[0] = 200;

    }
};

//Normal derivative (used for Neumann boundary condition)
class dUdN : public dolfin::Expression {
    void eval(dolfin::Array<double> &values, const dolfin::Array<double> &x) const {
        values[0] = x[0];
    }
};

//Sub domain for Dirichlet boundary condition
    class DirichletBoundary : public dolfin::SubDomain {
        bool inside(const dolfin::Array<double> &x, bool on_boundary) const {
            return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS or x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS;
        }
    };
#endif //DIFFUSION_FENICS_PDEHELPER_H

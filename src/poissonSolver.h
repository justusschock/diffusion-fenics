//
// Created by js on 13.05.16.
//

#ifndef DIFFUSION_FENICS_POISSONSOLVER_H
#define DIFFUSION_FENICS_POISSONSOLVER_H
#include <dolfin.h>
#include "poissonProblem.h"

namespace Poisson {


//sourceTerm (right-hand side)
    class Source : public dolfin::Expression {
        void eval(dolfin::Array<double> &values, const dolfin::Array<double> &x) const {
            values[0] = x[0]*sin(x[0]);
        }
    };

//Normal derivative (used for Neumann boundary condition)
    class dUdN : public dolfin::Expression {
        void eval(dolfin::Array<double> &values, const dolfin::Array<double> &x) const {
            values[0] = 0;
        }
    };

//Sub domain for Dirichlet boundary condition
    class DirichletBoundary : public dolfin::SubDomain {
        bool inside(const dolfin::Array<double> &x, bool on_boundary) const {
            return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS;
        }
    };

    auto solvePDE(std::shared_ptr<dolfin::Mesh> mesh, dolfin::Constant dirichletBoundary) -> dolfin::Function {
        //create FunctionSpace
        auto V = std::make_shared<poissonProblem::FunctionSpace>(mesh);

        //Define boundary condition
        auto u0 = std::make_shared<dolfin::Constant>(dirichletBoundary);
        auto boundary = std::make_shared<DirichletBoundary>();
        dolfin::DirichletBC bc(V, u0, boundary);

        //Define variational forms
        poissonProblem::BilinearForm a(V, V);
        poissonProblem::LinearForm L(V);


        //Set Boundary Condition for Problem
        auto f = std::make_shared<Source>();
        auto g = std::make_shared<dUdN>();
        L.f = *f;
        L.g = *g;

        dolfin::Function u(V);
        //Compute solution
        dolfin::solve(a == L, u, bc);

        return u;

    }
}

#endif //DIFFUSION_FENICS_POISSONSOLVER_H

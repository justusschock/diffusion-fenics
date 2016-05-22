//
// Created by js on 13.05.16.
//

#ifndef DIFFUSION_FENICS_POISSONSOLVER_H
#define DIFFUSION_FENICS_POISSONSOLVER_H
#include <dolfin.h>


namespace Poisson {

    class Initial : public dolfin::Expression {
        void eval(dolfin::Array<double> &values, const dolfin::Array<double> &x) const{
            values[0] = 1;
        }
    };

//sourceTerm (right-hand side)
    class Source : public dolfin::Expression {
        void eval(dolfin::Array<double> &values, const dolfin::Array<double> &x) const {
            values[0] = -1;
            if(x[0] > 0.5-DOLFIN_EPS and x[0] < 0.5+DOLFIN_EPS and x[1] > 0.5-DOLFIN_EPS and x[1] < 0.5+DOLFIN_EPS)
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

    //TODO: implement solver for multiple dimensions in this function (search solution for usage of different namespaces)
    template <class FunctionSpace, class LinearForm, class BilinearForm>
    auto solvePDE(std::shared_ptr<dolfin::Mesh> mesh, dolfin::Constant dirichletBoundary, dolfin::Expression& initial) -> dolfin::Function {

        //create FunctionSpace

        //default FunctionSpace and variational forms(necessary for setting derived classes for each dimension)
        auto V = std::make_shared<FunctionSpace>(mesh);
        BilinearForm a(V, V);
        LinearForm L(V);


        dolfin::Function u(V);

        //TODO: implement initial values (not working yet)
        u.interpolate(initial);


        //Define boundary condition
        auto u0 = std::make_shared<dolfin::Constant>(dirichletBoundary);
        auto boundary = std::make_shared<DirichletBoundary>();
        dolfin::DirichletBC bc(V, u0, boundary);


        //Set Boundary Condition for Problem
        Source f;
        dUdN g;
        L.g = g;
        L.f = f;

        //Compute solution
        dolfin::solve(a == L, u, bc);

        return u;

    }
}

#endif //DIFFUSION_FENICS_POISSONSOLVER_H

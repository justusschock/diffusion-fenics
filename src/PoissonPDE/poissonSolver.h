//
// Created by js on 13.05.16.
//

#ifndef DIFFUSION_FENICS_POISSONSOLVER_H
#define DIFFUSION_FENICS_POISSONSOLVER_H
#include <dolfin.h>
#include "../pdeTestExamples.h"

//Include generated headers for solving the Poisson-PDE
#include "poisson1D.h"
#include "poisson2D.h"
#include "poisson3D.h"


namespace Poisson {

    //Wrapper class for changing Dimensions of Poissond-Problem
    template <int> class DimensionWrapper;

    //Provide FunctionSpace, Linear and Bilinear Form in 1D
    template <> class DimensionWrapper<1> {
    public:
        auto FunctionSpace(std::shared_ptr<dolfin::Mesh> mesh) -> poisson1D::FunctionSpace
        {return poisson1D::FunctionSpace(mesh);}

        auto LinearForm(std::shared_ptr<poisson1D::FunctionSpace> FunctionSpace) -> poisson1D::LinearForm
        { return poisson1D::LinearForm(FunctionSpace);}

        auto BilinearForm(std::shared_ptr<poisson1D::FunctionSpace> FunctionSpace1,
                          std::shared_ptr<poisson1D::FunctionSpace> FunctionSpace2) -> poisson1D::BilinearForm
        {return poisson1D::BilinearForm(FunctionSpace1, FunctionSpace2);}

    };

    //Provide FunctionSpace, Linear and Bilinear Form in 2D
    template <> class DimensionWrapper<2> {
    public:
        auto FunctionSpace(std::shared_ptr<dolfin::Mesh> mesh) -> poisson2D::FunctionSpace
        {return poisson2D::FunctionSpace(mesh);}

        auto LinearForm(std::shared_ptr<poisson2D::FunctionSpace> FunctionSpace) -> poisson2D::LinearForm
        { return poisson2D::LinearForm(FunctionSpace);}

        auto BilinearForm(std::shared_ptr<poisson2D::FunctionSpace> FunctionSpace1,
                          std::shared_ptr<poisson2D::FunctionSpace> FunctionSpace2) -> poisson2D::BilinearForm
        {return poisson2D::BilinearForm(FunctionSpace1, FunctionSpace2);}

    };


    //Provide FunctionSpace, Linear and Bilinear Form in 3D
    template <> class DimensionWrapper<3> {
    public:
        auto FunctionSpace(std::shared_ptr<dolfin::Mesh> mesh) -> poisson3D::FunctionSpace
        {return poisson3D::FunctionSpace(mesh);}

        auto LinearForm(std::shared_ptr<poisson3D::FunctionSpace> FunctionSpace) -> poisson3D::LinearForm
        { return poisson3D::LinearForm(FunctionSpace);}

        auto BilinearForm(std::shared_ptr<poisson3D::FunctionSpace> FunctionSpace1,
                          std::shared_ptr<poisson3D::FunctionSpace> FunctionSpace2) -> poisson3D::BilinearForm
        {return poisson3D::BilinearForm(FunctionSpace1, FunctionSpace2);}

    };

    //TODO: implement solver for multiple dimensions in this function (search solution for usage of different namespaces)
    template <const int dim>
    auto solvePDE(std::shared_ptr<dolfin::Mesh> mesh, dolfin::Constant& dirichletBoundary, dolfin::Expression& initial,
                  dolfin::Expression& source, dolfin::Expression& neumann) -> dolfin::Function {

        DimensionWrapper<dim> dimensionWrapper;

        //Setup FunctionSpace, Linear and BilinearForm (based on dim)
        auto V = std::make_shared<decltype(dimensionWrapper.FunctionSpace(mesh))>(dimensionWrapper.FunctionSpace(mesh));
        auto a = std::make_shared<decltype(dimensionWrapper.BilinearForm(V,V))>(dimensionWrapper.BilinearForm(V,V));
        auto L = std::make_shared<decltype(dimensionWrapper.LinearForm(V))>(dimensionWrapper.LinearForm(V));

        //setup solution
        dolfin::Function u(V);

        //TODO: implement initial values (not working yet)
        u.interpolate(initial);

        //Define boundary condition
        auto u0 = std::make_shared<dolfin::Constant>(dirichletBoundary);
        auto boundary = std::make_shared<DirichletBoundary>();
        dolfin::DirichletBC bc(V, u0, boundary);

        auto g = std::make_shared<dolfin::Expression>(neumann);
        auto f = std::make_shared<dolfin::Expression>(source);
        //Set Boundary Condition for Problem
        L->g = g;
        L->f = f;

        dolfin::Equation equation(a, L);

        //Compute solution
        dolfin::solve(equation, u, bc);


        return u;

    }
}

#endif //DIFFUSION_FENICS_POISSONSOLVER_H

//
// Created by js on 13.05.16.
//

#ifndef DIFFUSION_FENICS_POISSONSOLVER_H
#define DIFFUSION_FENICS_POISSONSOLVER_H
#include <dolfin.h>
#include "../pdeHelper.h"

//Include generated headers for solving the Poisson-PDE
#include "poissonProblem1D.h"
#include "poissonProblem2D.h"
#include "poissonProblem3D.h"


namespace Poisson {

    //Wrapper class for changing Dimensions of Poissond-Problem
    template <int> class DimensionWrapper;

    //Provide FunctionSpace, Linear and Bilinear Form in 1D
    template <> class DimensionWrapper<1> {
    public:
        auto FunctionSpace(std::shared_ptr<dolfin::Mesh> mesh) -> poissonProblem1D::FunctionSpace
        {return poissonProblem1D::FunctionSpace(mesh);}

        auto LinearForm(std::shared_ptr<poissonProblem1D::FunctionSpace> FunctionSpace) -> poissonProblem1D::LinearForm
        { return poissonProblem1D::LinearForm(FunctionSpace);}

        auto BilinearForm(std::shared_ptr<poissonProblem1D::FunctionSpace> FunctionSpace1,
                          std::shared_ptr<poissonProblem1D::FunctionSpace> FunctionSpace2) -> poissonProblem1D::BilinearForm
        {return poissonProblem1D::BilinearForm(FunctionSpace1, FunctionSpace2);}

    };

    //Provide FunctionSpace, Linear and Bilinear Form in 2D
    template <> class DimensionWrapper<2> {
    public:
        auto FunctionSpace(std::shared_ptr<dolfin::Mesh> mesh) -> poissonProblem2D::FunctionSpace
        {return poissonProblem2D::FunctionSpace(mesh);}

        auto LinearForm(std::shared_ptr<poissonProblem2D::FunctionSpace> FunctionSpace) -> poissonProblem2D::LinearForm
        { return poissonProblem2D::LinearForm(FunctionSpace);}

        auto BilinearForm(std::shared_ptr<poissonProblem2D::FunctionSpace> FunctionSpace1,
                          std::shared_ptr<poissonProblem2D::FunctionSpace> FunctionSpace2) -> poissonProblem2D::BilinearForm
        {return poissonProblem2D::BilinearForm(FunctionSpace1, FunctionSpace2);}

    };


    //Provide FunctionSpace, Linear and Bilinear Form in 3D
    template <> class DimensionWrapper<3> {
    public:
        auto FunctionSpace(std::shared_ptr<dolfin::Mesh> mesh) -> poissonProblem3D::FunctionSpace
        {return poissonProblem3D::FunctionSpace(mesh);}

        auto LinearForm(std::shared_ptr<poissonProblem3D::FunctionSpace> FunctionSpace) -> poissonProblem3D::LinearForm
        { return poissonProblem3D::LinearForm(FunctionSpace);}

        auto BilinearForm(std::shared_ptr<poissonProblem3D::FunctionSpace> FunctionSpace1,
                          std::shared_ptr<poissonProblem3D::FunctionSpace> FunctionSpace2) -> poissonProblem3D::BilinearForm
        {return poissonProblem3D::BilinearForm(FunctionSpace1, FunctionSpace2);}

    };

    //TODO: implement solver for multiple dimensions in this function (search solution for usage of different namespaces)
    template <const int dim>
    auto solvePDE(std::shared_ptr<dolfin::Mesh> mesh, dolfin::Constant& dirichletBoundary, dolfin::Expression& initial,
                  dolfin::Expression& source, dolfin::Expression& neumann) -> dolfin::Function {

        DimensionWrapper<dim> dimensionWrapper;

        //Setup FunctionSpace, Linear and BilinearForm (based on dim)
        auto V = std::make_shared<decltype(dimensionWrapper.FunctionSpace(mesh))>(dimensionWrapper.FunctionSpace(mesh));
        auto a = dimensionWrapper.BilinearForm(V,V);
        auto L = dimensionWrapper.LinearForm(V);

        //setup solution
        dolfin::Function u(V);

        //TODO: implement initial values (not working yet)
        u.interpolate(initial);

        //Define boundary condition
        auto u0 = std::make_shared<dolfin::Constant>(dirichletBoundary);
        auto boundary = std::make_shared<DirichletBoundary>();
        dolfin::DirichletBC bc(V, u0, boundary);


        //Set Boundary Condition for Problem
        L.g = neumann;
        L.f = source;

        //Compute solution
        dolfin::solve(a == L, u, bc);


        return u;

    }
}

#endif //DIFFUSION_FENICS_POISSONSOLVER_H

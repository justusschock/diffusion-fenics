//
// Created by js on 03.07.16.
//

#ifndef DIFFUSION_FENICS_CONVECTIONDIFFUSIONSOLVER_H
#define DIFFUSION_FENICS_CONVECTIONDIFFUSIONSOLVER_H

#include <dolfin.h>

//Include generated headers for solving the Convection-Diffusion-PDE
//#include "convectionDiffusion1D.h"
//#include "convectionDiffusion2D.h"
#include "convectionDiffusion3D.h"
//#include "velocity1D.h"
//#include "velocity2D.h"
#include "velocity3D.h"
#include "pdeSetupClasses.h"

namespace ConvectionDiffusion{

    //Wrapper class for changing Dimensions of Convection-Diffusion-Problem
    template <int> class DimensionWrapper;

    /*
    //Provide FunctionSpace, Linear and Bilinear Form in 1D
    template <> class DimensionWrapper<1> {
    public:
        auto FunctionSpace(std::shared_ptr<dolfin::Mesh> mesh) -> convectionDiffusion1D::FunctionSpace
        {return convectionDiffusion1D::FunctionSpace(mesh);}

        auto LinearForm(std::shared_ptr<convectionDiffusion1D::FunctionSpace> FunctionSpace) -> convectionDiffusion1D::LinearForm
        { return convectionDiffusion1D::LinearForm(FunctionSpace);}

        auto BilinearForm(std::shared_ptr<convectionDiffusion1D::FunctionSpace> FunctionSpace1,
                          std::shared_ptr<convectionDiffusion1D::FunctionSpace> FunctionSpace2) -> convectionDiffusion1D::BilinearForm
        {return convectionDiffusion1D::BilinearForm(FunctionSpace1, FunctionSpace2);}

        auto VelocityFunctionSpace(std::shared_ptr<dolfin::Mesh> mesh) -> velocity1D::FunctionSpace
        {return velocity1D::FunctionSpace(mesh);}
    };

    //Provide FunctionSpace, Linear and Bilinear Form in 2D
    template <> class DimensionWrapper<2> {
    public:
        auto FunctionSpace(std::shared_ptr<dolfin::Mesh> mesh) -> convectionDiffusion2D::FunctionSpace
        {return convectionDiffusion2D::FunctionSpace(mesh);}

        auto LinearForm(std::shared_ptr<convectionDiffusion2D::FunctionSpace> FunctionSpace) -> convectionDiffusion2D::LinearForm
        { return convectionDiffusion2D::LinearForm(FunctionSpace);}

        auto BilinearForm(std::shared_ptr<convectionDiffusion2D::FunctionSpace> FunctionSpace1,
                          std::shared_ptr<convectionDiffusion2D::FunctionSpace> FunctionSpace2) -> convectionDiffusion2D::BilinearForm
        {return convectionDiffusion2D::BilinearForm(FunctionSpace2, FunctionSpace2);}

        auto VelocityFunctionSpace(std::shared_ptr<dolfin::Mesh> mesh) -> velocity2D::FunctionSpace
        {return velocity2D::FunctionSpace(mesh);}
    };
     */

    //Provide FunctionSpace, Linear and Bilinear Form in 2D
    template <> class DimensionWrapper<3> {
    public:
        auto FunctionSpace(std::shared_ptr<dolfin::Mesh> mesh) -> convectionDiffusion3D::FunctionSpace
        {return convectionDiffusion3D::FunctionSpace(mesh);}

        auto LinearForm(std::shared_ptr<convectionDiffusion3D::FunctionSpace> FunctionSpace) -> convectionDiffusion3D::LinearForm
        { return convectionDiffusion3D::LinearForm(FunctionSpace);}

        //::FIXME:: FunctionSpace1 wurde nicht verwendet
        auto BilinearForm(std::shared_ptr<convectionDiffusion3D::FunctionSpace> FunctionSpace1,
                          std::shared_ptr<convectionDiffusion3D::FunctionSpace> FunctionSpace2) -> convectionDiffusion3D::BilinearForm
        {return convectionDiffusion3D::BilinearForm(FunctionSpace1, FunctionSpace2);}

        auto VelocityFunctionSpace(std::shared_ptr<dolfin::Mesh> mesh) -> velocity3D::FunctionSpace
        {return velocity3D::FunctionSpace(mesh);}
    };

    template <const int dim, class SetupCase>
    void solvePDE(SetupCase& setup, dolfin::Constant _k = dolfin::Constant(0.01), const double T = 2.0, double t = 0.00)
    {
        DimensionWrapper<dim> dimensionwrapper;

        // Create velocity FunctionSpace and velocity function
        auto V_u = std::make_shared<decltype(dimensionwrapper.VelocityFunctionSpace(setup.getMesh()))>
                (dimensionwrapper.VelocityFunctionSpace(setup.getMesh()));
        //    auto velocity = std::make_shared<dolfin::Function>(V_u);
        //   *velocity = velocityFunction;

        // Create function space and function (to store solution)
        auto V = std::make_shared<decltype(dimensionwrapper.FunctionSpace(setup.getMesh()))>
                (dimensionwrapper.FunctionSpace(setup.getMesh()));

        setup.setU(new dolfin::Function(V));

        auto k =std::make_shared<dolfin::Constant>(_k);
        // Set up forms
        auto a = dimensionwrapper.BilinearForm(V, V);
        a.b = setup.getVelocity();
        a.c = setup.getDiffusivity();
        a.k = k;

        //Set velocityfunction, initial values and source term
        auto L = dimensionwrapper.LinearForm(V);

        auto ds = setup.getFacetFunction();
        auto dx = setup.getSubDomainFunction();
        L.u0 = setup.getInitial();
        L.b = setup.getVelocity();
        L.f = setup.getSource();
        L.g = setup.getNeumann();
        L.k = k;
        L.c = setup.getDiffusivity();

        a.ds = ds;
        L.ds = ds;
        a.dx = dx;
        L.dx = dx;

        // Set up boundary condition
        //dolfin::DirichletBC bc(*V, *setup.getDirichletValue(), ds,6); // works only without bc

        // Linear system
        std::shared_ptr<dolfin::Matrix> A(new dolfin::Matrix);
        dolfin::Vector b;

        // Assemble matrix
        assemble(*A, a);
        //bc.apply(*A);

        // LU solver
        dolfin::LUSolver lu(A);
        lu.parameters["reuse_factorization"] = true;

        double dt = _k;
        dolfin::File file ("../output/convection_diffusion.pvd","compressed");

        // Time-stepping
        dolfin::Progress p("Time-stepping");
        while (t < T)
        {
            std::cout<<"Simulating time: " << t << " of " << T << std::endl;
            // Assemble vector and apply boundary conditions
            assemble(b, L);
            //bc.apply(b);

            // Solve the linear system (re-use the already factorized matrix A)
            lu.solve(*(setup.getU()->vector()), b);

            file << std::pair<dolfin::Function*,double>(&(*setup.getU()),t);
            L.u0=setup.getU();

            // Move to next interval
            p = t/T;
            t += dt;
        }

    }

}
#endif //DIFFUSION_FENICS_CONVECTIONDIFFUSIONSOLVER_H



#ifndef DIFFUSION_FENICS_CURRENTSOLVER_H
#define DIFFUSION_FENICS_CURRENTSOLVER_H

#include <algorithm>
#include <functional>
#include <iterator>
#include <random>

#include "current3DExperimental.h"
#include "current3DLiquid.h"
#include "current3DSolid.h"

namespace Current {

    // Wrapper class for changing Dimensions
    template <int>
    class DimensionWrapper;

    /*  //Provide FunctionSpace, Linear and Bilinear Form in 1D
      template <> class DimensionWrapper<1> {
      public:
          auto FunctionSpace(std::shared_ptr<dolfin::Mesh> mesh) ->
      current1D::FunctionSpace
          {return current1D::FunctionSpace(mesh);}

          auto LinearForm(std::shared_ptr<current1D::FunctionSpace>
      FunctionSpace)
      -> current1D::LinearForm
          { return current1D::LinearForm(FunctionSpace);}

          auto BilinearForm(std::shared_ptr<current1D::FunctionSpace>
      FunctionSpace1,
                            std::shared_ptr<current1D::FunctionSpace>
      FunctionSpace2) -> current1D::BilinearForm
          {return current1D::BilinearForm(FunctionSpace1, FunctionSpace2);}

      };

      //Provide FunctionSpace, Linear and Bilinear Form in 2D
      template <> class DimensionWrapper<2> {
      public:
          auto FunctionSpace(std::shared_ptr<dolfin::Mesh> mesh) ->
      current2D::FunctionSpace
          {return current2D::FunctionSpace(mesh);}

          auto LinearForm(std::shared_ptr<current2D::FunctionSpace>
      FunctionSpace)
      -> current2D::LinearForm
          { return current2D::LinearForm(FunctionSpace);}

          auto BilinearForm(std::shared_ptr<current2D::FunctionSpace>
      FunctionSpace1,
                            std::shared_ptr<current2D::FunctionSpace>
      FunctionSpace2) -> current2D::BilinearForm
          {return current2D::BilinearForm(FunctionSpace1, FunctionSpace2);}

      };

    */
    /*
// Provide FunctionSpace, Linear and Bilinear Form in 3D
    template <> class DimensionWrapper<3> {
    public:
        auto FunctionSpaceLiquid(std::shared_ptr<dolfin::Mesh> mesh)
        -> current3DLiquid::FunctionSpace {
            return current3DLiquid::FunctionSpace(mesh);
        }

        auto FunctionSpaceSolid(std::shared_ptr<dolfin::Mesh> mesh)
                ->current3DSolid::FunctionSpace{
            return current3DSolid::FunctionSpace(mesh);
        }

        auto LinearFormSolid(std::shared_ptr<current3DSolid::FunctionSpace>
FunctionSpace)
        -> current3DSolid::LinearForm {
            return current3DSolid::LinearForm(FunctionSpace);
        }

        auto LinearFormLiquid(std::shared_ptr<current3DLiquid::FunctionSpace>
FunctionSpace)
                -> current3DLiquid::LinearForm {
            return current3DLiquid::LinearForm(FunctionSpace);
        }

        auto BilinearFormSolid(std::shared_ptr<current3DSolid::FunctionSpace>
FunctionSpace1,
                          std::shared_ptr<current3DSolid::FunctionSpace>
FunctionSpace2)
        -> current3DSolid::BilinearForm {
            return current3DSolid::BilinearForm(FunctionSpace1, FunctionSpace2);
        }

        auto BilinearFormLiquid(std::shared_ptr<current3DLiquid::FunctionSpace>
FunctionSpace1,
                          std::shared_ptr<current3DLiquid::FunctionSpace>
FunctionSpace2)
        -> current3DLiquid::BilinearForm {
            return current3DLiquid::BilinearForm(FunctionSpace1,
FunctionSpace2);
        }
    };

    template <const int dim, class SetupCase>
    auto solvePDE(SetupCase& setup) -> void {

        DimensionWrapper<dim> dimensionWrapper;

        // Setup FunctionSpace, Linear and BilinearForm (based on dim)
        auto V_l =
std::make_shared<decltype(dimensionWrapper.FunctionSpaceLiquid(setup.getMesh()))>(
                dimensionWrapper.FunctionSpaceLiquid(setup.getMesh()));
        auto V_s =
std::make_shared<decltype(dimensionWrapper.FunctionSpaceSolid(setup.getMesh()))>(
                dimensionWrapper.FunctionSpaceSolid((setup.getMesh())));

        auto a = dimensionWrapper.BilinearFormSolid(V_s, V_s);
        auto L = dimensionWrapper.LinearFormSolid(V_s);
        auto b = dimensionWrapper.BilinearFormLiquid(V_l, V_l);
        auto K = dimensionWrapper.LinearFormLiquid(V_l);

        setup.setUL(std::make_shared<dolfin::Function>(V_l));
        setup.setUS(std::make_shared<dolfin::Function>(V_s));

        //auto ul = dolfin::Function(V_l);
        //auto us = dolfin::Function(V_s);

        auto dx = setup.getSubdomainFunction();
        auto ds = setup.getFacetFunction();

        a.dx = dx;
        L.dx = dx;
        b.dx = dx;
        K.dx = dx;

        a.ds = ds;
        L.ds = ds;
        b.ds = ds;
        K.ds = ds;

        a.sigmaSolid = *setup.getSigmaS();
        b.sigmaLiquid = *setup.getSigmaL();
        L.fSolid = *setup.getSourceS();
        K.fLiquid = *setup.getSourceL();

        // Compute solutions
        dolfin::solve(a == L, *setup.getUS());
        dolfin::solve(b == K, *setup.getUL());

        dolfin::File file("../output/current_liquid.pvd", "compressed");
        file << *setup.getUL();

        dolfin::File file_s("../output/current_solve.pvd", "compressed");
        file_s << *setup.getUS();


    }
    */

    template <const int dim, class SetupCase>
    auto solvePDE(SetupCase& setup) -> void
    {
        /*        // First create an instance of an engine.
                std::random_device rnd_device;
                // Specify the engine and distribution.
                std::mt19937 mersenne_engine(rnd_device());
                std::uniform_real_distribution<double> dist(-0.5, 0.5);
                auto gen = std::bind(dist, mersenne_engine);
        */
        auto dx = setup.getSubdomainFunction();
        auto ds = setup.getFacetFunction();
        // Setup FunctionSpace, Linear and BilinearForm (based on dim)
        auto mesh = setup.getMesh();
        auto V = std::make_shared<current3DExperimental::FunctionSpace>(mesh);

        // auto V = std::make_shared<
        //     current3DExperimental::MultiMeshForm_F::CoefficientSpace_u>(
        //     mesh_multi);

        auto u = std::make_shared<dolfin::Function>(V);
        //         auto u_s = std::make_shared<dolfin::Funiction>(V_solid);

        /*        std::vector<double> tmp(u->vector()->size());
                std::generate(tmp.begin(), tmp.end(), gen);
                u->vector()->set_local(tmp);
                u->vector()->apply("");
        */
                u->interpolate(dolfin::Constant(0.9, 0.00, -0.9));

    //    auto umin = std::make_shared<dolfin::Function>(V);
    //    umin->interpolate(dolfin::Constant(0.0, -1.0, -2.0));

      //  auto umax = std::make_shared<dolfin::Function>(V);
      //  umax->interpolate(dolfin::Constant(2.0, 1.0, 0.0));

        auto F = std::make_shared<current3DExperimental::Form_F>(V);
        auto J = std::make_shared<current3DExperimental::Form_J>(V, V);

        //        auto u = std::make_shared<dolfin::MultiMeshFunction>(V);
        J->dx = dx;
        F->dx = dx;

 //       F->ds = ds;
 //       F->dS = ds;

   //     J->ds = ds;
   //     J->dS = ds;

        auto bc1 = std::make_shared<const dolfin::DirichletBC>(
            V->sub(0),
            std::make_shared<dolfin::Constant>(1.),
            ds,
            137,
            "geometric");
        auto bc2 = std::make_shared<const dolfin::DirichletBC>(
            V->sub(1),
            std::make_shared<dolfin::Constant>(0.),
            ds,
            139,
            "geometric");
        auto bc3 = std::make_shared<const dolfin::DirichletBC>(
            V->sub(2),
            std::make_shared<dolfin::Constant>(-1.),
            ds,
            136,
            "geometric");
        std::vector<std::shared_ptr<const dolfin::DirichletBC>> bcs{bc1,bc3};

        J->u = u;
        F->u = u;

        J->sigma = setup.getSigma();
        F->sigma = setup.getSigma();

        J->U_eq = setup.getU_eq();
        J->RT = setup.getRT();
        J->alpha = setup.getAlpha();
        J->i0 = setup.getI0();

        F->U_eq = setup.getU_eq();
        F->RT = setup.getRT();
        F->alpha = setup.getAlpha();
        F->i0 = setup.getI0();
        // F.n = setup.getNeumann();
        // F.k = std::make_shared<dolfin::Constant>(1e-7);

        // Compute solutions
        auto problem =
            std::make_shared<dolfin::NonlinearVariationalProblem>(F, u, bcs, J);
        dolfin::NonlinearVariationalSolver solver(problem);

        solver.parameters["nonlinear_solver"] = "snes";
        solver.parameters("snes_solver")["linear_solver"] = "bicgstab";
        solver.parameters("snes_solver")["maximum_iterations"] = 20;
        solver.parameters("snes_solver")["relative_tolerance"] = 1e-8;
        solver.parameters("snes_solver")["absolute_tolerance"] = 1e-8;
        solver.parameters("snes_solver")["preconditioner"] = "icc";
        solver.parameters("snes_solver")["report"] = true;
        solver.parameters("snes_solver")["line_search"] = "basic";
        solver.parameters("snes_solver")["error_on_nonconvergence"] = false;
        solver.parameters("snes_solver")("krylov_solver")["relative_tolerance"] = 1e-8;
        solver.parameters("snes_solver")("krylov_solver")["absolute_tolerance"] = 1e-8;
        solver.parameters("snes_solver")("krylov_solver")["error_on_nonconvergence"] = false;
        solver.parameters("snes_solver")("krylov_solver")["maximum_iterations"] = 100;
        //solver.parameters("snes_solver")["linear_solver"]["relative_tolerance"] = 1e-8;
        //dolfin::info(solver.parameters,true);
        solver.solve();
        dolfin::File file("../output/current.pvd");
        file << (*u);
        dolfin::File fileXML("../output/currrent.xml");
        fileXML << (*u);

        dolfin::File file_pos("../output/current_pos.pvd", "compressed");
        file_pos << (*u)[0];

        dolfin::File file_electrolyte("../output/current_electrolyte.pvd",
                                      "compressed");
        file_electrolyte << (*u)[1];

        dolfin::File file_neg("../output/current_neg.pvd", "compressed");
        file_neg << (*u)[2];

        dolfin::File file_domains("../output/current_domains.pvd",
                                  "compressed");
        file_domains << *dx;
        /*
        dolfin::Function uSolid = u[1];
        dolfin::Function uLiquid = u[0];

        fileSolid<<uSolid;
        fileLiquid << uLiquid;
        */
    }
}

#endif  // DIFFUSION_FENICS_CURRENTSOLVER_H

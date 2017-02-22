

#ifndef DIFFUSION_FENICS_CURRENTSOLVER_H
#define DIFFUSION_FENICS_CURRENTSOLVER_H

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
        auto dx = setup.getSubdomainFunction();
        auto ds = setup.getFacetFunction();
        // Setup FunctionSpace, Linear and BilinearForm (based on dim)
        auto mesh = setup.getMesh();
        /* dolfin::SubDomain solid_domain;
         dolfin::SubDomain electrolyte_domain;
         solid_domain.mark(*dx, 1);
         solid_domain.mark(*dx, 2);
         solid_domain.mark(*dx, 3);
         solid_domain.mark(*dx, 4);

         electrolyte_domain.mark(*dx, 2);
         electrolyte_domain.mark(*dx, 3);
         electrolyte_domain.mark(*dx, 5);

         auto mesh_solid =
             std::make_shared<dolfin::SubMesh>(*mesh, solid_domain);
         auto mesh_electroylte =
             std::make_shared<dolfin::SubMesh>(*mesh, electrolyte_domain);

         dolfin::MultiMesh mesh_multi;
         mesh_multi.add(mesh_solid);
         mesh_multi.add(mesh_electroylte);
         mesh_multi.build();*/
        auto V =
            std::make_shared<current3DExperimental::CoefficientSpace_u>(mesh);

        // auto a = dimensionWrapper.BilinearFormExperimental(V, V);
        //    auto L = dimensionWrapper.LinearFormExperimental(V);
        //    auto J = dimensionWrapper.JacobianFormExperimental(V, V);

        //        auto V = std::make_shared<
        //            current3DExperimental::MultiMeshForm_F::CoefficientSpace_u>(
        //            mesh_multi);

        auto u = std::make_shared<dolfin::Function>(V);
        //         auto u_s = std::make_shared<dolfin::Funiction>(V_solid);

        u->interpolate(dolfin::Constant(2.1e-2,1e-2));
        auto F = current3DExperimental::Form_F(V);
        auto J = current3DExperimental::Form_J(V, V);

        //        auto u = std::make_shared<dolfin::MultiMeshFunction>(V);
        J.dx = dx;
        F.dx = dx;

        J.ds = ds;
        F.ds = ds;

        J.dS = ds;
        F.dS = ds;
        dolfin::DirichletBC bc1(V->sub(0),
                                std::make_shared<dolfin::Constant>(2e-2),
                                ds,
                                2);
        std::vector<const dolfin::DirichletBC*> bcs{&bc1};

        J.u = u;
        J.sigma = setup.getSigma();
        J.U_eq = setup.getU_eq();
        J.RT = setup.getRT();
        J.alpha = setup.getAlpha();
        J.i0 = setup.getI0();

        F.u = u;
        F.sigma = setup.getSigma();
        F.U_eq = setup.getU_eq();
        F.RT = setup.getRT();
        F.alpha = setup.getAlpha();
        F.i0 = setup.getI0();
        F.n=setup.getNeumann();
        F.k=std::make_shared<dolfin::Constant>(1e-7);

        dolfin::Parameters params("nonlinear_variational_solver");
        dolfin::Parameters newton_params("newton_solver");
        newton_params.add("relative_tolerance", 1e-6);
        newton_params.add("linear_solver", "mumps");
        params.add(newton_params);

        /*dolfin::NewtonSolver newtonSolver;
        newtonSolver.parameters = params;

        // Linear system
        std::shared_ptr<dolfin::Matrix> A(new dolfin::Matrix);
        dolfin::Vector b;

        // Assemble matrix
        assemble(*A, a);
        assemble(b, L);
        bc0.apply(*A);
        bc0.apply(b);

        dolfin::NonlinearProblem;
        newtonSolver.solve(A,setup.getU()->vector(),b);
        */

        // Compute solutions
        dolfin::solve(F == 0, *u, bcs, J,params);

        // auto u = dolfin::Function(V);
        // dolfin::solve(a==L,u,bc);

        dolfin::File file("../output/current.pvd");
        file << (*u);
        dolfin::File fileXML("../output/currrent.xml");
        fileXML << (*u);
        dolfin::File file_electrolyte("../output/current_electrolyte.pvd",
                                      "compressed");
        dolfin::File file_solid("../output/current_solid.pvd", "compressed");
        file_solid << (*u)[0];
        file_electrolyte << (*u)[1];

        /*
        dolfin::Function uSolid = u[1];
        dolfin::Function uLiquid = u[0];

        fileSolid<<uSolid;
        fileLiquid << uLiquid;
        */
    }
}

#endif  // DIFFUSION_FENICS_CURRENTSOLVER_H

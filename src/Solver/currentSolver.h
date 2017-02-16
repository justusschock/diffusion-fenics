

#ifndef DIFFUSION_FENICS_CURRENTSOLVER_H
#define DIFFUSION_FENICS_CURRENTSOLVER_H

#include "current3DSolid.h"
#include "current3DLiquid.h"
#include "current3DExperimental.h"

namespace Current {

// Wrapper class for changing Dimensions
    template <int> class DimensionWrapper;

/*  //Provide FunctionSpace, Linear and Bilinear Form in 1D
  template <> class DimensionWrapper<1> {
  public:
      auto FunctionSpace(std::shared_ptr<dolfin::Mesh> mesh) ->
  current1D::FunctionSpace
      {return current1D::FunctionSpace(mesh);}

      auto LinearForm(std::shared_ptr<current1D::FunctionSpace> FunctionSpace)
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

      auto LinearForm(std::shared_ptr<current2D::FunctionSpace> FunctionSpace)
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

        auto LinearFormSolid(std::shared_ptr<current3DSolid::FunctionSpace> FunctionSpace)
        -> current3DSolid::LinearForm {
            return current3DSolid::LinearForm(FunctionSpace);
        }

        auto LinearFormLiquid(std::shared_ptr<current3DLiquid::FunctionSpace> FunctionSpace)
                -> current3DLiquid::LinearForm {
            return current3DLiquid::LinearForm(FunctionSpace);
        }

        auto BilinearFormSolid(std::shared_ptr<current3DSolid::FunctionSpace> FunctionSpace1,
                          std::shared_ptr<current3DSolid::FunctionSpace> FunctionSpace2)
        -> current3DSolid::BilinearForm {
            return current3DSolid::BilinearForm(FunctionSpace1, FunctionSpace2);
        }

        auto BilinearFormLiquid(std::shared_ptr<current3DLiquid::FunctionSpace> FunctionSpace1,
                          std::shared_ptr<current3DLiquid::FunctionSpace> FunctionSpace2)
        -> current3DLiquid::BilinearForm {
            return current3DLiquid::BilinearForm(FunctionSpace1, FunctionSpace2);
        }
    };

    template <const int dim, class SetupCase>
    auto solvePDE(SetupCase& setup) -> void {

        DimensionWrapper<dim> dimensionWrapper;

        // Setup FunctionSpace, Linear and BilinearForm (based on dim)
        auto V_l = std::make_shared<decltype(dimensionWrapper.FunctionSpaceLiquid(setup.getMesh()))>(
                dimensionWrapper.FunctionSpaceLiquid(setup.getMesh()));
        auto V_s = std::make_shared<decltype(dimensionWrapper.FunctionSpaceSolid(setup.getMesh()))>(
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

    template <> class DimensionWrapper<3> {
    public:
        auto FunctionSpaceExperimental(std::shared_ptr<dolfin::Mesh> mesh)
        -> current3DExperimental::FunctionSpace {
            return current3DExperimental::FunctionSpace(mesh);
        }

        auto LinearFormExperimental(std::shared_ptr<current3DExperimental::FunctionSpace> FunctionSpace)
        -> current3DExperimental::LinearForm {
            return current3DExperimental::LinearForm(FunctionSpace);
        }

        auto BilinearFormExperimental(std::shared_ptr<current3DExperimental::FunctionSpace> FunctionSpace1,
                                      std::shared_ptr<current3DExperimental::FunctionSpace> FunctionSpace2)
        -> current3DExperimental::BilinearForm {
            return current3DExperimental::BilinearForm(FunctionSpace1, FunctionSpace2);
        }

        auto JacobianFormExperimental(std::shared_ptr<current3DExperimental::FunctionSpace> FunctionSpace1,
                                      std::shared_ptr<current3DExperimental::FunctionSpace> FunctionSpace2)
        -> current3DExperimental::JacobianForm {
            return current3DExperimental::JacobianForm(FunctionSpace1,FunctionSpace2);
        }

    };

    template <const int dim, class SetupCase>
    auto solvePDE(SetupCase& setup) -> void {

        DimensionWrapper<dim> dimensionWrapper;

        // Setup FunctionSpace, Linear and BilinearForm (based on dim)
        auto V = std::make_shared<decltype(dimensionWrapper.FunctionSpaceExperimental(setup.getMesh()))>(
                dimensionWrapper.FunctionSpaceExperimental((setup.getMesh())));

        //auto a = dimensionWrapper.BilinearFormExperimental(V, V);
        auto L = dimensionWrapper.LinearFormExperimental(V);
        auto J = dimensionWrapper.JacobianFormExperimental(V,V);

        setup.setU(std::make_shared<dolfin::Function>(V));

        //auto ul = dolfin::Function(V_l);
        //auto us = dolfin::Function(V_s);



        auto dx = setup.getSubdomainFunction();
        auto ds = setup.getFacetFunction();

        //a.dx = dx;
        L.dx = dx;
        J.dx = dx;

        //a.ds = ds;
        L.ds = ds;
        J.ds = ds;

        dolfin::DirichletBC bc0(V, std::make_shared<dolfin::Constant>(0.0, 0.0), ds, 11);
        dolfin::DirichletBC bc1(V, std::make_shared<dolfin::Constant>(6.0, 6.0), ds, 12);
        std::vector<const dolfin::DirichletBC*> bcs{&bc0, &bc1};

        //a.sigma = setup.getSigma();
        L.f = setup.getSource();
        L.g = setup.getNeumann();
        J.sigma = setup.getSigma();

        dolfin::Parameters params("nonlinear_variational_solver");
        dolfin::Parameters newton_params("newton_solver");
        newton_params.add("relative_tolerance", 1e-6);
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
        dolfin::solve(L==0, *setup.getU(), bc, J, params);

        //auto u = dolfin::Function(V);
        //dolfin::solve(a==L,u,bc);


        dolfin::File file("../../output/current.pvd");
        file << (*(setup.getU().get()));
        dolfin::File fileXML("../../output/currrent.xml");
        fileXML << (*(setup.getU().get()));
        std::cout<<(*setup.getU()).value_rank()<<std::endl;
        dolfin::File fileSolid("../../output/currentSolid.pvd", "compressed");
        dolfin::File fileLiquid("../../output/currentLiquid.pvd", "compressed");
        fileSolid << (*(setup.getU().get()))[1];
        fileLiquid << (*(setup.getU().get()))[0];


        /*
        dolfin::Function uSolid = u[1];
        dolfin::Function uLiquid = u[0];

        fileSolid<<uSolid;
        fileLiquid << uLiquid;
        */

    }
}

#endif //DIFFUSION_FENICS_CURRENTSOLVER_H



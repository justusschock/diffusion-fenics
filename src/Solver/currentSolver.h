

#ifndef DIFFUSION_FENICS_CURRENTSOLVER_H
#define DIFFUSION_FENICS_CURRENTSOLVER_H

#include "current3DSolid.h"
#include "current3DLiquid.h"

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

        setup.setUL(new dolfin::Function(V_s));
        setup.setUS(new dolfin::Function(V_l));

        //auto ul = dolfin::Function(V_l);
        //auto us = dolfin::Function(V_s);

        auto dx = setup.getSubdomainFunction();
        auto ds = setup.getFacetFunction();

        a.dx = *dx;
        L.dx = *dx;
        b.dx = *dx;
        K.dx = *dx;

        a.ds = *ds;
        L.ds = *ds;
        b.ds = *ds;
        K.ds = *ds;

        a.sigmaSolid = *setup.getSigmaS();
        b.sigmaLiquid = *setup.getSigmaL();
        L.fSolid = *setup.getSourceS();
        K.fLiquid = *setup.getSourceL();

        // Compute solutions
        dolfin::solve(a == L, *setup.getUS());
        dolfin::solve(b == K, *setup.getUL());

        dolfin::File file("../output/current.pvd", "compressed");
        file << *setup.getUS();
        file << *setup.getUL();


    }
}

#endif //DIFFUSION_FENICS_CURRENTSOLVER_H

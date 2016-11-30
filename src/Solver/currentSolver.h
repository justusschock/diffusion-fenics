

#ifndef DIFFUSION_FENICS_CURRENTSOLVER_H
#define DIFFUSION_FENICS_CURRENTSOLVER_H

#include "UFL/current3D.h"

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
        auto FunctionSpace(std::shared_ptr<dolfin::Mesh> mesh)
        -> current3D::FunctionSpace {
            return current3D::FunctionSpace(mesh);
        }

        auto LinearForm(std::shared_ptr<current3D::FunctionSpace> FunctionSpace)
        -> current3D::LinearForm {
            return current3D::LinearForm(FunctionSpace);
        }

        auto BilinearForm(std::shared_ptr<current3D::FunctionSpace> FunctionSpace1,
                          std::shared_ptr<current3D::FunctionSpace> FunctionSpace2)
        -> current3D::BilinearForm {
            return current3D::BilinearForm(FunctionSpace1, FunctionSpace2);
        }
    };

    template <const int dim, class SetupCase>
    auto solvePDE(SetupCase& setup) -> void {

        DimensionWrapper<dim> dimensionWrapper;

        // Setup FunctionSpace, Linear and BilinearForm (based on dim)
        auto V = std::make_shared<decltype(dimensionWrapper.FunctionSpace(setup.getMesh()))>(
                dimensionWrapper.FunctionSpace(setup.getMesh()));
        auto a = dimensionWrapper.BilinearForm(V, V);
        auto L = dimensionWrapper.LinearForm(V);
        auto b = dimensionWrapper.BilinearForm(V, V);
        auto K = dimensionWrapper.LinearForm(V);

        setup.setU1(new dolfin::Function(V));
        setup.setU2(new dolfin::Function(V));

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

        a.sigmaSolid = *setup.getSigma1();
        b.sigmaLiquid = *setup.getSigma2();
        L.fSolid = *setup.getSource1();
        K.fLiquid = *setup.getSource2();


        // Compute solutions
        dolfin::solve(a == L, *setup.getU1());
        dolfin::solve(b == K, *setup.getU2());

        dolfin::File file("../output/current.pvd", "compressed");
        file << *setup.getU1();
        file << *setup.getU2();


    }
}

#endif //DIFFUSION_FENICS_CURRENTSOLVER_H

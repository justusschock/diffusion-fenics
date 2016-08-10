//
// Created by js on 05.08.16.
//

#ifndef DIFFUSION_FENICS_CURRENTSOLVER_H
#define DIFFUSION_FENICS_CURRENTSOLVER_H
#include <dolfin.h>

//Include generated headers for solving div(j) = 0 <=> div(-sigma*grad(phi) = 0
#include "current1D.h"
#include "current2D.h"
#include "current3D.h"

namespace Current {

	//Wrapper class for changing Dimensions
	template <int> class DimensionWrapper;

	//Provide FunctionSpace, Linear and Bilinear Form in 1D
	template <> class DimensionWrapper<1> {
	public:
		auto FunctionSpace(std::shared_ptr<dolfin::Mesh> mesh) -> current1D::FunctionSpace
		{return current1D::FunctionSpace(mesh);}

		auto LinearForm(std::shared_ptr<current1D::FunctionSpace> FunctionSpace) -> current1D::LinearForm
		{return current1D::LinearForm(FunctionSpace);}

		auto BilinearForm(std::shared_ptr<current1D::FunctionSpace> FunctionSpace1,
				  std::shared_ptr<current1D::FunctionSpace> FunctionSpace2) -> current1D::BilinearForm
		{return current1D::BilinearForm(FunctionSpace1, FunctionSpace2);}
	};

	//Provide FunctionSpace, Linear and Bilinear Form in 2D
	template <> class DimensionWrapper<2> {
	public:
		auto FunctionSpace(std::shared_ptr<dolfin::Mesh> mesh) -> current2D::FunctionSpace
		{return current2D::FunctionSpace(mesh);}

		auto LinearForm(std::shared_ptr<current2D::FunctionSpace> FunctionSpace) -> current2D::LinearForm
		{return current2D::LinearForm(FunctionSpace);}

		auto BilinearForm(std::shared_ptr<current2D::FunctionSpace> FunctionSpace1,
				  std::shared_ptr<current2D::FunctionSpace> FunctionSpace2) -> current2D::BilinearForm
		{return current2D::BilinearForm(FunctionSpace1, FunctionSpace2);}
	};

//Provide FunctionSpace, Linear and Bilinear Form in 3D
	template <> class DimensionWrapper<3> {
	public:
		auto FunctionSpace(std::shared_ptr<dolfin::Mesh> mesh) -> current3D::FunctionSpace
		{return current3D::FunctionSpace(mesh);}

		auto LinearForm(std::shared_ptr<current3D::FunctionSpace> FunctionSpace) -> current3D::LinearForm
		{return current3D::LinearForm(FunctionSpace);}

		auto BilinearForm(std::shared_ptr<current3D::FunctionSpace> FunctionSpace1,
				  std::shared_ptr<current3D::FunctionSpace> FunctionSpace2) -> current3D::BilinearForm
		{return current3D::BilinearForm(FunctionSpace1, FunctionSpace2);}
	};

	template<const int dim>
	auto solvePDE(std::shared_ptr<dolfin::Mesh> mesh, std::shared_ptr<dolfin::Constant> dirichletValue, std::shared_ptr<dolfin::expression> source,
		      std::shared_ptr<dolfin::Expression> neumann, std::shared_ptr<dolfin::SubDomain> dirichletBoundary) -> dolfin::Function {

		DimensionWrapper<dim> dimensionWrapper;

		//Setup FunctionSpace, Linear and Bilinear Form (based on dim)
		auto V = std::make_shared<decltype(dimensionWrapper.FunctionSpace(mesh))>(dimensionWrapper.FunctionSpace(mesh));
		auto a = dimensionWrapper.BilinearForm(V,V);
		auto L = dimensionWrapper.LinearForm(V);

		//setup solution
		dolfin::Function u(V);

		//Define boundary condition
		dolfin::DirichletBC bc(V;dirichletValue, dirichletBoundary);

		//Set Boundary Condition for Problem
		L.g = *neumann;
		L.f = *source;

		//Compute solution
		dolfin::solve(a==L, u, bc);

		dolfin::File file ("../output/current.pvd","compressed");
		file << u;

		return u;

	}

}
#endif //DIFFUSION_FENICS_CURRENTSOLVER_H

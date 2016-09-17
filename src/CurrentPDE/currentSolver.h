//
// Created by js on 05.08.16.
//

#ifndef DIFFUSION_FENICS_CURRENTSOLVER_H
#define DIFFUSION_FENICS_CURRENTSOLVER_H
#include <dolfin.h>
#include <list>

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
		auto FunctionSpaceCG(std::shared_ptr<dolfin::Mesh> mesh) -> current1D::FunctionSpace
		{return current1D::FunctionSpaceCG(mesh);}

		auto LinearFormCG(std::shared_ptr<current1D::FunctionSpace> FunctionSpace) -> current1D::LinearForm
		{return current1D::LinearForm(FunctionSpace);}

		auto BilinearFormCG(std::shared_ptr<current1D::FunctionSpace> FunctionSpace1,
				  std::shared_ptr<current1D::FunctionSpace> FunctionSpace2) -> current1D::BilinearForm
		{return current1D::BilinearForm(FunctionSpace1, FunctionSpace2);}
	};

	//Provide FunctionSpace, Linear and Bilinear Form in 2D
	template <> class DimensionWrapper<2> {
	public:
		auto FunctionSpaceCG(std::shared_ptr<dolfin::Mesh> mesh) -> current2D::FunctionSpace
		{return current2D::FunctionSpace(mesh);}

		auto LinearFormCG(std::shared_ptr<current2D::FunctionSpace> FunctionSpace) -> current2D::LinearForm
		{return current2D::LinearForm(FunctionSpace);}

		auto BilinearFormCG(std::shared_ptr<current2D::FunctionSpace> FunctionSpace1,
				  std::shared_ptr<current2D::FunctionSpace> FunctionSpace2) -> current2D::BilinearForm
		{return current2D::BilinearForm(FunctionSpace1, FunctionSpace2);}
	};

//Provide FunctionSpace, Linear and Bilinear Form in 3D
	template <> class DimensionWrapper<3> {
	public:
		auto FunctionSpaceCG(std::shared_ptr<dolfin::Mesh> mesh) -> current3D::FunctionSpace
		{return current3D::FunctionSpace(mesh);}

		auto LinearFormCG(std::shared_ptr<current3D::FunctionSpace> FunctionSpace) -> current3D::LinearForm
		{return current3D::LinearForm(FunctionSpace);}

		auto BilinearFormCG(std::shared_ptr<current3D::FunctionSpace> FunctionSpace1,
				  std::shared_ptr<current3D::FunctionSpace> FunctionSpace2) -> current3D::BilinearForm
		{return current3D::BilinearForm(FunctionSpace1, FunctionSpace2);}
	};

	template<const int dim>
	auto solvePDE(std::shared_ptr<dolfin::Mesh> mesh, std::shared_ptr<dolfin::Constant> dirichletValue, std::shared_ptr<dolfin::expression> source,
		      std::shared_ptr<dolfin::Expression> neumann, std::shared_ptr<dolfin::SubDomain> dirichletBoundary, 
		      std::shared_ptr<dolfin::SubDomain> elektrolyt, std::shared_ptr<dolfin::SubDomain> conductors, 
		      std::shared_ptr<dolfin::SubDomain> firstElectrode, std::shared_ptr<dolfin::SubDomain> secondElectrode, 
		      std::shared_ptr<dolfin::Subdomain) -> dolfin::Function 
	{
		
		DimensionWrapper<dim> dimensionwrapper;

		//Create functionSpace and function for solution (in CG)
		auto W = std::make_shared<decltype(dimensionwrapper.FunctionSpaceCG(mesh))>(dimensionwrapper.FunctionSpace(mesh));
		auto w = dolfin::Function(V);

		//Set up CG forms
		auto aCG = dimensionWrapper.BilinearFormCG(W,W);
		auto LCG = dimensionWrapper.LinearFormCG(W);
		LCG.f = *source;
		LCG.g = *neumann;


		

		dolfin::DirichletBC bc(V, dirichletValue, dirichletBoundary);
		dolfin::CellFunction<std::size_t> sub_domains(mesh, mesh->topology().dim -1);
		
		sub_domains = 4;
		
		//Mark SubDomains
		elektrolyt->mark(sub_domains,0);
		conductors->mark(sub_domains,1);
		firstElectrode->mark(sub_domains,2);
		secondElectrode->mark(sub_domains,3);
		
		aCG.dx = sub_domains;
		LCG.dx = sub_domains;

		//Linear System
		std::shared_ptr<dolfin::Matrix> A (new dolfin::Matrix);
		dolfin::Vector b;

		//Assemble Matrix
		assemble(*A,aCG);
		bc.apply(*A);
		A.ident_zeros();
		

		//Assemble Vector and apply boundary conditions
		assemble(b, LCG);
		bc.apply(b);
		b.ident_zeros();

		dolfin::solve(A, w.vector(), b);
		
		//Create FunctionSoace and function for solution (in DG)
		auto V = //not implemented yet
		auto u = dolfin::Function(V);

		//Create DG Forms
		auto aDG = //not implemented yet
		auto LDG = //not implemented yet

		dolfin::solve(aDG == LDG, u);



		dolfin::File file("../output/current.pvd","compressed");
		file << u;




}
#endif //DIFFUSION_FENICS_CURRENTSOLVER_H

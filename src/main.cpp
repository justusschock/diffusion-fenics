#include <iostream>
#include <dolfin.h>

//includes for PDEs:
#include "PoissonPDE/poissonSolver.h"
#include "ConvectionDiffusionPDE/convectionDiffusionSolver.h"
#include "pdeTestExamples.h"


int main(int argc, char* argv[]) {
    try {

        enum{
            poisson,
            diffusion
        };

	enum{	
		general,
		inOut,
		constSides,
		randomSource
	};
        int equation = diffusion;
	int testCase = randomSource;

        dolfin::init(argc, argv);

        //dimension
        const int dim = 2;


        //initial mesh and solution (used for overwriting it for each dimension)
        std::shared_ptr<dolfin::Mesh> mesh;

        std::shared_ptr<dolfin::Expression> initial;
        std::shared_ptr<dolfin::Constant> dirichlet;
        std::shared_ptr<dolfin::Expression> neumann;
        std::shared_ptr<dolfin::Expression> source;       
	std::shared_ptr<dolfin::SubDomain> dirichletBoundary;
	std::shared_ptr<dolfin::Expression> velocity;
	std::shared_ptr<dolfin::Expression> diffusivity;
/*
	if(testCase == inOut && equation == diffusion){ 
		ConvectionDiffusion::CaseInOut caseInOut(dim);
		initial = caseInOut.getInitial();
		dirichletBoundary = caseInOut.getDirichletBoundary();
		neumann = caseInOut.getNeumann();
		source = caseInOut.getSource();
		dirichlet.reset (new dolfin::Constant(0.0));
		velocity = caseInOut.getVelocity();
		diffusivity = caseInOut.getDiffusivity();
		mesh.reset(new dolfin::UnitSquareMesh(50,50));// = caseInOut.getMesh();
		
	}
*/
	if(testCase == inOut && equation == diffusion) {
		initial.reset(new ConvectionDiffusion::InOut::Initial);
		dirichlet.reset(new dolfin::Constant(0.0));
		neumann.reset(new ConvectionDiffusion::InOut::Neumann);
		source.reset(new ConvectionDiffusion::InOut::Source);
		dirichletBoundary.reset(new ConvectionDiffusion::InOut::DirichletBoundary);
		velocity.reset(new ConvectionDiffusion::InOut::Velocity(dim));
		diffusivity.reset (new ConvectionDiffusion::InOut::Diffusivity);
		mesh.reset(new dolfin::UnitSquareMesh(50,50));
	}
	else if(testCase == constSides && equation == diffusion) {
		initial.reset(new ConvectionDiffusion::ConstSides::Initial);
		dirichlet.reset(new dolfin::Constant(0.9));
		neumann.reset(new ConvectionDiffusion::ConstSides::Neumann);
		source.reset(new ConvectionDiffusion::ConstSides::Source);
		dirichletBoundary.reset(new ConvectionDiffusion::ConstSides::DirichletBoundary);
		velocity.reset(new ConvectionDiffusion::ConstSides::Velocity(dim));
		diffusivity.reset(new ConvectionDiffusion::ConstSides::Diffusivity);
		mesh.reset(new dolfin::UnitSquareMesh(50,50));
	}
	else if(testCase == randomSource && equation == diffusion) {
		initial.reset(new ConvectionDiffusion::RandomSource::Initial);
		dirichlet.reset(new dolfin::Constant(0.0));
		neumann.reset(new ConvectionDiffusion::RandomSource::Neumann);
		source.reset(new ConvectionDiffusion::RandomSource::Source);
		dirichletBoundary.reset(new ConvectionDiffusion::RandomSource::DirichletBoundary);
		velocity.reset(new ConvectionDiffusion::RandomSource::Velocity(dim));
		diffusivity.reset(new ConvectionDiffusion::RandomSource::Diffusivity);
		mesh.reset(new dolfin::UnitSquareMesh(50,50));
	}
	else if(testCase == general){

	//set meshes and solutions for different dimensions (1D-3D)
        	if(dim == 1)
            		mesh.reset(new dolfin::UnitIntervalMesh(32));
       		else if(dim == 2)
            		mesh.reset(new dolfin::UnitSquareMesh(32, 32));
        	else if (dim == 3)
            		mesh.reset(new dolfin::UnitCubeMesh(32,32,32));
        	else
            		throw std::string("Wrong dimension Parameter");

		initial.reset(new TestInitial);
		dirichlet.reset(new dolfin::Constant(0.0));
		neumann.reset(new TestdUdN);
		source.reset(new TestSource);
		dirichletBoundary.reset(new TestDirichletBoundary);
		velocity.reset(new TestVelocity(dim));
		diffusivity.reset(new TestDiffusivity);
	}


        if(equation == poisson) {
            dolfin::Function u = Poisson::solvePDE<dim>(mesh, dirichlet, initial, source, neumann, dirichletBoundary);
            //dolfin::plot(u);
            //dolfin::interactive();
        }
        else if(equation == diffusion){
            	dolfin::Constant dt(0.5);
		double T = 50.0;
		dolfin::Function u = ConvectionDiffusion::solvePDE<dim>(mesh,dirichlet,initial,velocity,source,neumann,dirichletBoundary,diffusivity,dt,T);
            	//dolfin::plot(u);
            	//dolfin::interactive();
        }
        else
            throw std::string("Wrong equation selected");



        return 0;

    }
    catch (std::exception &e) {
        std::cout << e.what() << std::endl;
        return 1;
    }
    catch (std::string &e) {
        std::cout << e << std::endl;
        return 1;
    }
    catch (...) {
        std::cout << "An unknown error has occured" << std::endl;
        return 1;
    }
}

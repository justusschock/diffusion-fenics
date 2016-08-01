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
            		randomSource,
            		differentDiffusivities,
            		transferFunction
        	};
        	int equation = diffusion;
        	int testCase = differentDiffusivities;

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

		dolfin::Constant dt(1e-1);
        	double T = 5.0;

		if(testCase == inOut && equation == diffusion){ 
			ConvectionDiffusion::SetupCase<ConvectionDiffusion::InOut> setup(dim);
			dolfin::Function u = ConvectionDiffusion::solvePDE<dim>(setup.getMesh(),dirichlet,setup.getInitial(),setup.getVelocity(),setup.getSource(),
				    						setup.getNeumann(), setup.getDirichletBoundary(), setup.getDiffusivity(),
										dt, T);
		}
       
        	else if(testCase == constSides && equation == diffusion) {
        	   	ConvectionDiffusion::SetupCase<ConvectionDiffusion::ConstSides> setup(dim);
			dolfin::Function u = ConvectionDiffusion::solvePDE<dim>(setup.getMesh(),dirichlet,setup.getInitial(),setup.getVelocity(),setup.getSource(),
				    						setup.getNeumann(), setup.getDirichletBoundary(), setup.getDiffusivity(),
										dt, T);
        	}
        	else if(testCase == randomSource && equation == diffusion) {
        	    	ConvectionDiffusion::SetupCase<ConvectionDiffusion::RandomSource> setup(dim);
			dolfin::Function u = ConvectionDiffusion::solvePDE<dim>(setup.getMesh(),dirichlet,setup.getInitial(),setup.getVelocity(),setup.getSource(),
				    						setup.getNeumann(), setup.getDirichletBoundary(), setup.getDiffusivity(),
										dt, T);
        	}
        	else if(testCase == differentDiffusivities && equation == diffusion) {
        	    	ConvectionDiffusion::SetupCase<ConvectionDiffusion::DifferentDiffusivities> setup(dim);
			dolfin::Function u = ConvectionDiffusion::solvePDE<dim>(setup.getMesh(),dirichlet,setup.getInitial(),setup.getVelocity(),setup.getSource(),
				    						setup.getNeumann(), setup.getDirichletBoundary(), setup.getDiffusivity(),
										dt, T);
        	}
        	else if(testCase == transferFunction && equation == diffusion) {
        	    	ConvectionDiffusion::SetupCase<ConvectionDiffusion::TransferFunction> setup(dim);
			dolfin::Function u = ConvectionDiffusion::solvePDE<dim>(setup.getMesh(),dirichlet,setup.getInitial(),setup.getVelocity(),setup.getSource(),
				    						setup.getNeumann(), setup.getDirichletBoundary(), setup.getDiffusivity(),
										dt, T);
        	}
        	else if(testCase == general){
        	    	ConvectionDiffusion::SetupCase<ConvectionDiffusion::General> setup(dim);
			dolfin::Function u = ConvectionDiffusion::solvePDE<dim>(setup.getMesh(),dirichlet,setup.getInitial(),setup.getVelocity(),setup.getSource(),
				    						setup.getNeumann(), setup.getDirichletBoundary(), setup.getDiffusivity(),
										dt, T);
        	}
        	else if(equation == poisson) {
		    	Poisson::SetupCase<Poisson::General> setup(dim);
        	    	dolfin::Function u = Poisson::solvePDE<dim>(mesh, dirichlet, initial, source, neumann, dirichletBoundary);
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

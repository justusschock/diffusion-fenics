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

        int equation = diffusion;

        dolfin::init(argc, argv);

        //dimension
        const int dim = 2;


        //initial mesh and solution (used for overwriting it for each dimension)
        auto mesh = std::make_shared<dolfin::Mesh>();


        //set meshes and solutions for different dimensions (1D-3D)
        if(dim == 1) {
            mesh.reset(new dolfin::UnitIntervalMesh(32));
        }
        else if(dim == 2) {
            mesh.reset(new dolfin::UnitSquareMesh(32, 32));
        }
        else if (dim == 3) {
            mesh.reset(new dolfin::UnitCubeMesh(32,32,32));
        }
        else
            throw std::string("Wrong dimension Parameter");
        Initial initial;
        dolfin::Constant dirichlet(0.0);
        dUdN neumann;
        Source source;
        auto velocity=std::make_shared<Velocity>(dim);


        if(equation == poisson) {
            dolfin::Function u = Poisson::solvePDE<dim>(mesh, dirichlet, initial, source, neumann);
            dolfin::plot(u);
            dolfin::interactive();
        }
        else if(equation == diffusion){
            dolfin::Function u = convectionDiffusion::solvePDE<dim>(mesh,dirichlet,initial,velocity,source, neumann);
            dolfin::plot(u);
            dolfin::interactive();
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

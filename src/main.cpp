#include <iostream>

//includes for Poisson-PDE:
#include "poissonSolver.h"


int main(int argc, char* argv[]) {
    try {

        dolfin::init(argc, argv);

        //dimension and initial values
        const int dim = 3;
        Poisson::Initial initial;
        dolfin::Constant dirichlet(0.0);

        dolfin::FunctionSpace fSpace();
        //initial mesh and solution (used for overwriting it in each dimension
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

        auto u = std::make_shared<dolfin::Function>(Poisson::solvePDE<dim>(mesh,dirichlet,initial));
        dolfin::plot(*u);
        dolfin::interactive();


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

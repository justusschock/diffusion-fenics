#include <iostream>

//includes for Poisson-PDE:
#include "poissonSolver.h"
#include "poissonProblem1D.h"
#include "poissonProblem2D.h"
#include "poissonProblem3D.h"


int main(int argc, char* argv[]) {
    try {

        dolfin::init(argc, argv);

        //dimension and initial values
        const int dim = 2;
        Poisson::Initial initial;


        //set meshes and solutions for different dimensions (1D-3D)
        if(dim == 1) {
            using namespace poissonProblem1D;
            auto mesh = std::make_shared<dolfin::UnitSquareMesh>(32, 32);
            auto u = std::make_shared<dolfin::Function>
                    (Poisson::solvePDE<FunctionSpace, LinearForm, BilinearForm>(mesh, dolfin::Constant(0), initial));
            dolfin::plot(*u);
            dolfin::interactive();
        }
        else if(dim == 2) {
            using namespace poissonProblem2D;
            auto mesh = std::make_shared<dolfin::UnitSquareMesh>(32, 32);
            auto u = std::make_shared<dolfin::Function>
                    (Poisson::solvePDE<FunctionSpace, LinearForm, BilinearForm>(mesh, dolfin::Constant(0), initial));
            dolfin::plot(*u);
            dolfin::interactive();
        }
        else if (dim == 3) {
            using namespace poissonProblem3D;
            auto mesh = std::make_shared<dolfin::UnitSquareMesh>(32, 32);
            auto u = std::make_shared<dolfin::Function>
                    (Poisson::solvePDE<FunctionSpace, LinearForm, BilinearForm>(mesh, dolfin::Constant(0), initial));
            dolfin::plot(*u);
            dolfin::interactive();
        }
        else
            throw std::string("Wrong dimension Parameter");

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
    }
}

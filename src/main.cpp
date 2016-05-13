#include <iostream>
#include "poissonSolver.h"


int main() {
    try {

        //setup mesh and solution
        auto mesh = std::make_shared<dolfin::UnitSquareMesh>(32, 32);

        dolfin::Function u = Poisson::solvePDE(mesh, dolfin::Constant(0.0));

        dolfin::File file("poisson.pvd");
        file << u;

        dolfin::plot(u);
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
    }
}

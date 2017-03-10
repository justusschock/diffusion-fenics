#include <dolfin.h>
#include <iostream>
#include </usr/include/dolfin/la/solve.h>

// includes for PDEs:
//#include "Solver/poissonSolver.h"
//#include "Solver/convectionDiffusionSolver.h"
#include "pdeSetupClasses.h"
#include "Solver/currentSolver.h"

int main(int argc, char *argv[]) {
    try {

        dolfin::init(argc, argv);
        dolfin::list_linear_solver_methods();
        dolfin::list_krylov_solver_preconditioners();

        // dimension
        const int dim = 3;


        std::string meshFile = "../test_mesh/test_mesh.xml";
        std::string subdomainFile = "../test_mesh/test_mesh_physical_region.xml";
        std::string facetFile = "../test_mesh/test_mesh_facet_region.xml";

        std::cout<<"Simulation running"<<std::endl;
        Current::SetupCase<Current::General> setup(dim,meshFile,subdomainFile,facetFile);

        Current::solvePDE<dim, Current::SetupCase<Current::General>>(setup);
        // dolfin::plot(u);
        // dolfin::interactive();
    } catch (std::exception &e) {
        std::cout << e.what() << std::endl;
        return 1;
    } catch (std::string &e) {
        std::cout << e << std::endl;
        return 1;
    } catch (...) {
        std::cout << "An unknown error has occured" << std::endl;
        return 1;
    }
}

#include <dolfin.h>
#include <iostream>

// includes for PDEs:
#include "PoissonPDE_Mesh/poissonSolver.h"
#include "ConvectionDiffusionPDE_Mesh/convectionDiffusionSolver.h"
#include "pdeTestExamples.h"

int main(int argc, char *argv[]) {
    try {

        dolfin::init(argc, argv);

        // dimension
        const int dim = 3;

        // initial mesh and solution (used for overwriting it for each dimension)
        std::shared_ptr<dolfin::Mesh> mesh =
                std::make_shared<dolfin::Mesh>("../test_mesh/test_mesh.xml");

        auto subdomain_function = std::make_shared<dolfin::MeshFunction<size_t>>(
                mesh, "../test_mesh/test_mesh_physical_region.xml");
        auto facet_function = std::make_shared<dolfin::MeshFunction<size_t>>(
                mesh, "../test_mesh/test_mesh_facet_region.xml");

        dolfin::File file2("../output/subdomains.pvd", "compressed");
        file2 << *subdomain_function;
        dolfin::File file3("../output/facets.pvd", "compressed");
        file3 << *facet_function;
        dolfin::Constant dt(1e-1);
        double T = 5.0;

        // dirichlet.reset(new dolfin::Constant(0.9));

        /*Poisson::SetupCase<Poisson::General> setup(dim);
        std::shared_ptr<dolfin::Expression> neumann = setup.getNeumann();
        std::shared_ptr<dolfin::Expression> source = setup.getSource();

         */

        std::cout<<"Simulation running"<<std::endl;
        ConvectionDiffusion::SetupCase<ConvectionDiffusion::General> setup(dim);
        std::shared_ptr<dolfin::Expression> neumann = setup.getNeumann();
        std::shared_ptr<dolfin::Expression> source = setup.getSource();
        std::shared_ptr<dolfin::Expression> initial= setup.getInitial();
        std::shared_ptr<dolfin::Expression> velocity = setup.getVelocity();
        std::shared_ptr<dolfin::Expression> diffusivity = setup.getDiffusivity();

        dolfin::Function u =
                ConvectionDiffusion::solvePDE<dim>(mesh, source, neumann, subdomain_function,facet_function,initial,velocity,diffusivity);
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

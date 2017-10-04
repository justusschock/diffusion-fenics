

#ifndef DIFFUSION_FENICS_CURRENTSOLVER_H
#define DIFFUSION_FENICS_CURRENTSOLVER_H

#include <algorithm>
#include <functional>
#include <iterator>
#include <random>

#include "current3DFixPotential.h"
#include "current3DFixCurrent.h"

namespace Current {


    void setSolverParameters(dolfin::Variable& solver, int maxIt = 20, double relTol = 1e-8, double absTol = 1e-8,
                             double relTolKrylov = 1e-8, double absTolKrylov = 1e-8){
        solver.parameters["nonlinear_solver"] = "snes";
        solver.parameters("snes_solver")["linear_solver"] = "bicgstab";
        solver.parameters("snes_solver")["maximum_iterations"] = maxIt;
        solver.parameters("snes_solver")["relative_tolerance"] = relTol;
        solver.parameters("snes_solver")["absolute_tolerance"] = absTol;
        solver.parameters("snes_solver")["preconditioner"] = "icc";
        solver.parameters("snes_solver")["report"] = true;
        solver.parameters("snes_solver")["line_search"] = "basic";
        solver.parameters("snes_solver")["error_on_nonconvergence"] = false;
        solver.parameters("snes_solver")("krylov_solver")["relative_tolerance"] = relTolKrylov;
        solver.parameters("snes_solver")("krylov_solver")["absolute_tolerance"] = absTolKrylov;
        solver.parameters("snes_solver")("krylov_solver")["error_on_nonconvergence"] = false;
        solver.parameters("snes_solver")("krylov_solver")["maximum_iterations"] = 1000;
        solver.parameters("snes_solver")("krylov_solver")["report"] = false;
    }

    template <const int dim, class SetupCaseFixPotential, class SetupCaseFixCurrent>
    auto solvePDE(SetupCaseFixPotential& setupFixPotential, SetupCaseFixCurrent& setupFixCurrent, dolfin::Constant _k = dolfin::Constant(0.01),
                  const double T = 2.0, double t = 0.00, std::vector<double> caseToggleTimes = std::vector<double>()) -> void
    {

        auto dxFixPotential = setupFixPotential.getSubdomainFunction();
        auto dsFixPotential = setupFixPotential.getFacetFunction();

        auto dxFixCurrent = setupFixCurrent.getSubdomainFunction();
        auto dsFixCurrent = setupFixCurrent.getFacetFunction();

        // Setup FunctionSpace, Linear and BilinearForm (based on dim)
        auto meshFixPotential = setupFixPotential.getMesh();
        auto V = std::make_shared<current3DFixPotential::FunctionSpace>(meshFixPotential);

        auto u = std::make_shared<dolfin::Function>(V);

        u->interpolate(dolfin::Constant(0.9, 0.00, -0.9));

        //    auto umin = std::make_shared<dolfin::Function>(V);
        //    umin->interpolate(dolfin::Constant(0.0, -1.0, -2.0));

        //  auto umax = std::make_shared<dolfin::Function>(V);
        //  umax->interpolate(dolfin::Constant(2.0, 1.0, 0.0));

        auto FFixPotential = std::make_shared<current3DFixPotential::Form_F>(V);
        auto JFixPotential = std::make_shared<current3DFixPotential::Form_J>(V, V);

        auto FFixCurrent = std::make_shared<current3DFixCurrent::Form_F>(V);
        auto JFixCurrent = std::make_shared<current3DFixCurrent::Form_J>(V, V);

        FFixPotential->dx = dxFixPotential;
        JFixPotential->dx = dxFixPotential;

        FFixCurrent->dx = dxFixCurrent;
        JFixCurrent->dx = dxFixCurrent;

        FFixCurrent->ds = dsFixCurrent;
        JFixCurrent->ds = dsFixCurrent;

        //       F->ds = ds;
        //       F->dS = ds;

        //     J->ds = ds;
        //     J->dS = ds;

        auto bc1 = std::make_shared<const dolfin::DirichletBC>(
                V->sub(0),
                std::make_shared<dolfin::Constant>(1.),
                dsFixPotential,
                137,
                "geometric");
        auto bc2 = std::make_shared<const dolfin::DirichletBC>(
                V->sub(1),
                std::make_shared<dolfin::Constant>(0.),
                dsFixPotential,
                139,
                "geometric");
        auto bc3 = std::make_shared<const dolfin::DirichletBC>(
                V->sub(2),
                std::make_shared<dolfin::Constant>(-1.),
                dsFixPotential,
                136,
                "geometric");
        std::vector<std::shared_ptr<const dolfin::DirichletBC>> dirichletBCsFixPotential{bc1,bc3};
        std::vector<std::shared_ptr<const dolfin::DirichletBC>> dirichletBCsFixCurrent{bc3};

        FFixPotential->u = u;
        JFixPotential->u = u;

        FFixCurrent->u = u;
        JFixCurrent->u = u;

        FFixPotential->sigma = setupFixPotential.getSigma();
        JFixPotential->sigma = setupFixPotential.getSigma();

        FFixCurrent->sigma = setupFixCurrent.getSigma();
        JFixCurrent->sigma = setupFixCurrent.getSigma();


        FFixPotential->U_eq = setupFixPotential.getU_eq();
        FFixPotential->RT = setupFixPotential.getRT();
        FFixPotential->alpha = setupFixPotential.getAlpha();
        FFixPotential->i0 = setupFixPotential.getI0();

        JFixPotential->U_eq = setupFixPotential.getU_eq();
        JFixPotential->RT = setupFixPotential.getRT();
        JFixPotential->alpha = setupFixPotential.getAlpha();
        JFixPotential->i0 = setupFixPotential.getI0();

        FFixCurrent->U_eq = setupFixCurrent.getU_eq();
        FFixCurrent->RT = setupFixCurrent.getRT();
        FFixCurrent->alpha = setupFixCurrent.getAlpha();
        FFixCurrent->i0 = setupFixCurrent.getI0();
        FFixCurrent->neumann = setupFixCurrent.getNeumann();

        JFixCurrent->U_eq = setupFixCurrent.getU_eq();
        JFixCurrent->RT = setupFixCurrent.getRT();
        JFixCurrent->alpha = setupFixCurrent.getAlpha();
        JFixCurrent->i0 = setupFixCurrent.getI0();

        // Compute solutions
        auto problemFixPotential =
                std::make_shared<dolfin::NonlinearVariationalProblem>(FFixPotential, u, dirichletBCsFixPotential, JFixPotential);
        auto problemFixCurrent =
                std::make_shared<dolfin::NonlinearVariationalProblem>(FFixCurrent, u, dirichletBCsFixCurrent, JFixCurrent);
        dolfin::NonlinearVariationalSolver solverFixPotential(problemFixPotential);
        dolfin::NonlinearVariationalSolver solverFixCurrent(problemFixCurrent);

        setSolverParameters(solverFixCurrent);
        setSolverParameters(solverFixPotential);

        //If not dirichlet Case use Neumann Case
        bool dirichletCase = true;

        std::sort(caseToggleTimes.begin(),caseToggleTimes.end());
        if(caseToggleTimes.empty() or *caseToggleTimes.end() < T) caseToggleTimes.push_back(T);

        auto dt = _k;
        dolfin::Progress p("Time-stepping");
        for (auto i : caseToggleTimes){
            dt = _k;
            while(t <= i){
                if(t+dt >i) dt = i-t;

                std::cout << "Simulating time: " << t << " of " << T << std::endl;
                if(dirichletCase){
                    solverFixPotential.solve();
                }
                else {
                    solverFixCurrent.solve();
                }

                p = t / T;
                //dt = rungeKutta(u, dt, meshFixPotential);
                t += dt;
                setupFixPotential.setTime(t);
                setupFixCurrent.setTime(t);
            }
            dirichletCase = not dirichletCase;
        }

        dolfin::File file("../output/current.pvd");
        file << (*u);
        dolfin::File fileXML("../output/currrent.xml");
        fileXML << (*u);

        dolfin::File file_pos("../output/current_pos.pvd", "compressed");
        file_pos << (*u)[0];

        dolfin::File file_electrolyte("../output/current_electrolyte.pvd",
                                      "compressed");
        file_electrolyte << (*u)[1];

        dolfin::File file_neg("../output/current_neg.pvd", "compressed");
        file_neg << (*u)[2];

        dolfin::File file_domains("../output/current_domains.pvd",
                                  "compressed");
        file_domains << *dxFixPotential;
        /*
        dolfin::Function uSolid = u[1];
        dolfin::Function uLiquid = u[0];

        fileSolid<<uSolid;
        fileLiquid << uLiquid;
        */
    }
}

#endif  // DIFFUSION_FENICS_CURRENTSOLVER_H

//
// Created by js on 22.04.17.
//

#ifndef PROJECT_RUNGEKUTTA_H
#define PROJECT_RUNGEKUTTA_H

#include <dolfin.h>
#include <cmath>
#include "rungeKuttaDiffusion.h"

#define EPS_RUNGE_KUTTA 1e-10

// using anonymous namespace to limit scope to current file
namespace {
    // Define Arrays and constant factors
    const unsigned int N = 6;

    std::array<std::shared_ptr<dolfin::FunctionAXPY>, N> k;
    std::array<std::array<double, 5>, N> factors = {
        std::array<double, 5>{0.0, 0.0, 0.0, 0.0, 0.0},
        std::array<double, 5>{0.25, 0.0, 0.0, 0.0, 0.0},
        std::array<double, 5>{3.0 / 32.0, 9.0 / 32.0, 0.0, 0.0, 0.0},
        std::array<double, 5>{
            1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0, 0.0, 0.0},
        std::array<double, 5>{
            439.0 / 216.0, 8.0, 3680.0 / 513.0, 845.0 / 4104.0, 0.0},
        std::array<double, 5>{
            -8.0 / 27.0, 2.0, -3544.0 / 2565.0, 1859.0 / 4104, -11.0 / 40.0}};
    std::array<double, N> eval_factors = {
        EPS_RUNGE_KUTTA, 1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 0.5};
    std::array<double, N> y_factors = {
        25.0 / 216.0, 0.0, 1408.0 / 2565.0, 2197.0 / 4104.0, -1.0 / 5.0, 0.0};
    std::array<double, N> z_factors = {16.0 / 135.0,
                                       0.0,
                                       6656.0 / 12825.0,
                                       28561.0 / 56430.0,
                                       -9.0 / 50.0,
                                       2.0 / 55.0};
}

auto evaluate_func(std::shared_ptr<dolfin::Mesh> mesh,
                   std::shared_ptr<dolfin::MeshFunction<size_t>> facets,
                   std::shared_ptr<dolfin::MeshFunction<size_t>> subdomains,
                   std::shared_ptr<dolfin::Function> y,
                   std::shared_ptr<dolfin::Expression> diffusivity,
                   std::shared_ptr<dolfin::Expression> velocity,
                   std::shared_ptr<dolfin::Expression> source,
                    std::shared_ptr<dolfin::Expression> neumann)
    -> std::shared_ptr<dolfin::Function>
{
    std::shared_ptr<dolfin::Function> u =
        std::make_shared<dolfin::Function>(y->function_space());

    // Setup Linear and Bilinear Form
    auto V = std::make_shared<dolfin::FunctionSpace>(
        rungeKuttaDiffusion::FunctionSpace(mesh));
    auto L = rungeKuttaDiffusion::LinearForm(V);
    auto a = rungeKuttaDiffusion::BilinearForm(V, V);

    // Assign Coefficients
    L.y = y;
    L.D = diffusivity;
    L.b = velocity;
    L.S = source;
    L.g = neumann;

    a.ds = facets;
    L.ds = facets;
    a.dx = subdomains;
    L.dx = subdomains;

    // Create Solver
    std::shared_ptr<dolfin::Matrix> A(new dolfin::Matrix);
    dolfin::Vector b;

    dolfin::assemble(*A, a);
    dolfin::assemble(b, L);

    dolfin::LUSolver lu(A);
    lu.parameters["reuse_factorization"] = true;
    lu.solve(*(u->vector()), b);

    return u;
}

double calc_l2_norm(dolfin::Function &func)
{
    std::vector<double> values;
    func.vector()->get_local(values);
    std::size_t counter = 0;
    for (unsigned int i = 0; i < values.size(); i++) {
        if (values[i] != values[i]) {
            values[i] = 0.0;
        }
        else {
            counter++;
        }
    }
    double norm = 0;
    for (auto i : values) {
        norm += i * i;
    }
    norm /= counter;
    return std::sqrt(norm);
}
void calc_l1_norm(dolfin::Function &func, double &norm, unsigned int &counter)
{
    std::vector<double> values;
    func.vector()->get_local(values);
    counter = 0;
    for (unsigned int i = 0; i < values.size(); i++) {
        if (values[i] != values[i]) {
            values[i] = 0.0;
        }
        else {
            counter++;
        }
    }
    norm = 0;
    for (auto i : values) {
        norm += std::abs(i);
    }
}

// FIXME Add velocity and diffusivity term
double rungeKuttaFifthOrder(
    std::shared_ptr<dolfin::Mesh> mesh,
    std::shared_ptr<dolfin::MeshFunction<size_t>> facets,
    std::shared_ptr<dolfin::MeshFunction<size_t>> subdomains,
    std::shared_ptr<dolfin::Function> y,
    std::shared_ptr<dolfin::Expression> source,
    std::shared_ptr<dolfin::Expression> diffusivity,
    std::shared_ptr<dolfin::Expression> velocity,
    std::shared_ptr<dolfin::Expression> neumann,
    std::shared_ptr<dolfin::Function> u,
    double max_stepsize,
    double eps_rel,
    double eps_abs)
{
    while (true) {
        auto source_fkt =
            std::make_shared<dolfin::Function>(y->function_space());
        source_fkt->interpolate(*source);
        auto z_next = dolfin::FunctionAXPY(y, 1.0);
        auto y_next = dolfin::FunctionAXPY(y, 1.0);
        for (unsigned int i = 0; i < N; i++) {
            auto ev_func_axpy = dolfin::FunctionAXPY(y, 1.0);
            if (i > 0) {
                for (unsigned int j = 0; j < i - 1; j++) {
                    if (j > 0) {
                        ev_func_axpy =
                            ev_func_axpy + (*k.at(j) * factors.at(i).at(j));
                    }
                }
            }
            std::shared_ptr<dolfin::Function> rhs =
                std::make_shared<dolfin::Function>(y->function_space());
            *rhs = ev_func_axpy;

            //        auto dt = max_stepsize * eval_factors.at(i);

            auto tmp = std::make_shared<dolfin::FunctionAXPY>(
                evaluate_func(mesh,
                              facets,
                              subdomains,
                              rhs,
                              diffusivity,
                              velocity,
                              source, neumann),
                1.0);
            *tmp = (*tmp + source_fkt) * max_stepsize;
            k.at(i) = tmp;

            y_next = y_next + (*k.at(i) * y_factors.at(i));
            z_next = z_next + (*k.at(i) * z_factors.at(i));
        }

        std::vector<double> values_z;
        auto z_func = dolfin::Function(y->function_space());
        auto diff = dolfin::Function(y->function_space());
        z_func = z_next;
        diff = z_next - y_next;
        const double norm_diff = calc_l2_norm(diff);
        const double norm_y = calc_l2_norm(*y);

        if (norm_diff >= (eps_rel * norm_y + eps_abs)) {
            // manually calculate norm after set NaNs to zero
            double s = (eps_rel * norm_y + eps_abs) / (2.0 * norm_diff);
            s = std::pow(s, 0.25);
            max_stepsize *= s;
        }
        else {
            u = std::make_shared<dolfin::Function>(z_func);
            break;
        }
    }
    return max_stepsize;
}

#endif  // PROJECT_RUNGEKUTTA_H

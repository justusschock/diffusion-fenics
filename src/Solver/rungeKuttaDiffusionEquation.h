//
// Created by js on 22.04.17.
//

#ifndef PROJECT_RUNGEKUTTA_H
#define PROJECT_RUNGEKUTTA_H

#include <dolfin.h>
#include <cmath>
#include "convectionDiffusion3D.h"



auto evaluate_func(std::shared_ptr<dolfin::Function> u0, double start_time, double eval_time,
                   std::shared_ptr<convectionDiffusion3D::Form_L> L, std::shared_ptr<convectionDiffusion3D::Form_a> a)
-> std::shared_ptr<dolfin::Function> {
    std::shared_ptr<dolfin::Function> u = std::make_shared<dolfin::Function>(u0->function_space());
    L->u0 = u0;
    double dt = start_time-eval_time;
    auto k = std::make_shared<dolfin::Constant>(dt);
    L->k = k;
    a->k = k;
    std::shared_ptr<dolfin::Matrix> A(new dolfin::Matrix);
    dolfin::Vector b;

    dolfin::assemble(*A, *a);

    dolfin::LUSolver lu(A);
    lu.parameters["reuse_factorization"] = true;
    lu.solve(*u->vector(), b);

    return u;
}


double rungeKuttaFifthOrder(std::shared_ptr<dolfin::Function> y, double max_stepsize,
                            double tol, double start_time, std::shared_ptr<convectionDiffusion3D::Form_L> L,
                            std::shared_ptr<convectionDiffusion3D::Form_a> a) {
    const unsigned int N = 6;

    std::array<std::shared_ptr<dolfin::FunctionAXPY>, N> k;
    std::array<std::array<double, 5>, N> factors = {std::array<double, 5>{0.0, 0.0, 0.0, 0.0, 0.0},
                                                    std::array<double, 5>{0.25, 0.0, 0.0, 0.0, 0.0},
                                                    std::array<double, 5>{3.0 / 32.0, 9.0 / 32.0, 0.0, 0.0, 0.0},
                                                    std::array<double, 5>{1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0, 0.0, 0.0},
                                                    std::array<double, 5>{439.0 / 216.0, 8.0, 3680.0 / 513.0, 845.0 / 4104.0, 0.0},
                                                    std::array<double, 5>{-8.0 / 27.0, 2.0, -3544.0 / 2565.0, 1859.0 / 4104,
                                                                          -11.0 / 40.0}};
    std::array<double, N> eval_factors = {0, 1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 0.5};
    std::array<double, N> y_factors = {25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0, 0.0};
    std::array<double, N> z_factors = {16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0};
    auto z_next = dolfin::FunctionAXPY(y, 1.0);
    auto y_next = dolfin::FunctionAXPY(y, 1.0);
    for (unsigned int i = 0; i < N; i++) {
        auto ev_func_axpy = dolfin::FunctionAXPY(y, 1.0);
        if (i > 0) {
            for (unsigned int j = 0; j < i - 1; j++) {

                if (j > 0)
                    ev_func_axpy = ev_func_axpy + (*k.at(j) * factors.at(i).at(j));
            }
        }
        std::shared_ptr<dolfin::Function> ev_func = std::make_shared<dolfin::Function>(y->function_space());
        *ev_func = ev_func_axpy;

        auto eval_offset = max_stepsize * eval_factors.at(i);
        auto eval_time = start_time + eval_offset;

        k.at(i) = std::make_shared<dolfin::FunctionAXPY>(evaluate_func(ev_func, start_time, eval_time, L, a), 1.0);

        y_next = y_next + (*k.at(i)* y_factors.at(i));
        z_next = z_next + (*k.at(i)* z_factors.at(i));

    }

    auto diff = z_next-y_next;
    auto diff_func = dolfin::Function(y->function_space());
    diff_func = diff;
    double norm_diff = dolfin::norm(*diff_func.vector(), std::string("l1"));
    double s = tol/(2.0*norm_diff);
    s = pow(s, 0.25);

    return s;
}
#endif //PROJECT_RUNGEKUTTA_H

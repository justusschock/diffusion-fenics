//
// Created by js on 22.04.17.
//

#ifndef PROJECT_RUNGEKUTTA_H
#define PROJECT_RUNGEKUTTA_H

#include <dolfin.h>
#include "UFL/gradientCurrentEquation.h"
#include "UFL/l2ErrorCurrentEquation.h"
#include <cmath>

/*
class CustomFunction: public dolfin::Function{

    explicit CustomFunction(std::shared_ptr<const dolfin::FunctionSpace> V)
            : dolfin::Function(V) { };

    CustomFunction(std::shared_ptr<const dolfin::FunctionSpace> V,
                   std::shared_ptr<dolfin::GenericVector> x) : dolfin::Function(V, x) { };

    CustomFunction(std::shared_ptr<const dolfin::FunctionSpace> V,
                   std::string filename) : dolfin::Function(V, filename) { };

    CustomFunction(const dolfin::Function& v) : dolfin::Function(v) { };

    CustomFunction(const dolfin::Function& v, std::size_t i) : dolfin::Function(v, i) { };

    CustomFunction operator*(dolfin::Function func, dolfin::Mesh mesh){
        dolfin::Function f(this->function_space())

        for (dolfin::MeshEntityIterator e(mesh, 0); !e.end(); ++e)
        {
            auto tmp = func.vector();
            for (dolfin::MeshEntityIterator e(mesh, 0); !e.end(); ++e)
            {
                for (dolfin::MeshEntityIterator e(mesh, 0); !e.end(); ++e)
                {
                    f[e->index()] = (*this)[e->index()]*func[e->index()]
                }
            }
        }
    }
};
 */

auto calc_grad(std::shared_ptr<dolfin::Function> f, std::shared_ptr<dolfin::Mesh> mesh, double grad_faktor) -> std::shared_ptr<dolfin::Function>
{
    using gradient=gradientCurrentEquation;
    auto V = std::make_shared<gradient::FunctionSpace>(mesh);
    auto L = std::make_shared<gradient::LinearForm>(V);
    auto a = std::make_shared<gradient::BilinearForm>(V, V);
    auto grad_func = std::make_shared<dolfin::Function>(V);

    L->u = f;
    L->factor = grad_faktor;
    dolfin::solve(a == L, grad_func);

    return grad_func;
}

double rungeKuttaFifthOrder(std::shared_ptr<dolfin::Function> y, double max_stepsize, std::shared_ptr<dolfin::Mesh> mesh, double tol=1e-3) {
    unsigned int N = 6;

    std::vector<std::shared_ptr<dolfin::FunctionAXPY>> k(N);
    std::vector<std::vector<double>> factors = {std::vector<double>{0.0},
                                                std::vector<double>{0.25},
                                                std::vector<double>{3.0 / 32.0, 9.0 / 32.0},
                                                std::vector<double>{1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0},
                                                std::vector<double>{439.0 / 216.0, 8.0, 3680.0 / 513.0, 845.0 / 4104.0},
                                                std::vector<double>{-8.0 / 27.0, 2.0, -3544.0 / 2565.0, 1859.0 / 4104,
                                                                    -11.0 / 40.0}};
    std::vector<double> grad_factors = {0, 1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 0.5};
    std::vector<double> y_factors = {25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0, 0.0};
    std::vector<double> z_factors = {16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0};
    auto z_next = dolfin::FunctionAXPY(y, 1.0);
    auto y_next = dolfin::FunctionAXPY(y, 1.0);
    for (unsigned int i = 0; i < N; i++) {
        auto grad_fkt_axpy = dolfin::FunctionAXPY(y, 1.0);
        for (unsigned int j = 0; j < i - 1; j++) {
            grad_fkt_axpy = grad_fkt_axpy + (*k.at(j) * factors.at(i).at(j));
        }
        std::shared_ptr<dolfin::Function> grad_fkt = std::make_shared<dolfin::Function>(grad_fkt_axpy);

        auto grad_factor = max_stepsize * grad_factors.at(i);
        k.at(i) = calc_grad(grad_fkt, mesh, grad_factor);

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

//
// Created by js on 22.04.17.
//

#ifndef PROJECT_RUNGEKUTTA_H
#define PROJECT_RUNGEKUTTA_H

#include <dolfin.h>
#include "UFL/gradientCurrentEquation.h"
#include "UFL/l2ErrorCurrentEquation.h"
#include <cmath>

double rungeKutta(std::shared_ptr<dolfin::Function> y, double max_stepsize, std::shared_ptr<dolfin::Mesh> mesh, double tol=1e-3){
    double b[] = {16/135, 0, 6656/12825, 28561/56430, (-9)/50, 2/55};
    double b_[] = {25/216, 0, 1408/2565, 2197/4104, (-1)/5, 0};

    // TODO Find Way to multiply Functions
    auto V = std::make_shared<gradientCurrentEquation::FunctionSpace>(mesh);

    auto gradY1 = std::make_shared<dolfin::Function>(V);
    auto L    = std::make_shared<gradientCurrentEquation::LinearForm>(V);
    auto a    = std::make_shared<gradientCurrentEquation::BilinearForm>(V,V);
    L->u = y;

    dolfin::solve(a == L, gradY1);

    dolfin::Function k1((*gradY1)*max_stepsize);

    auto gradY2 = std::make_shared<dolfin::Function>(V);
    L->u = std::make_shared<dolfin::Function>((*y+(k1*(1/3.0))+((*gradY1)*k1)*(max_stepsize/18.0));
    dolfin::solve(a==L, gradY2);

    dolfin::Function k2((*gradY2)*max_stepsize);

    auto gradY3 = std::make_shared<dolfin::Function>(V);
    L->u = std::make_shared<dolfin::Function>(*y-(k1*152/125.0)+(k2*252/125.0)-(((*gradY1)*k1)*(max_stepsize*44/125.0)));
    dolfin::solve(a==L, gradY3);
    dolfin::Function k3((*gradY3)*max_stepsize);

    auto gradY4 = std::make_shared<dolfin::Function>(V);
    L->u = std::make_shared<dolfin::Function>(*y+(k1*19/2.0)-(k2*72/7.0)+(k3*25/14.0)+(((*gradY1)*k1)*(max_stepsize*15/2.0)));
    dolfin::solve(a==L, gradY4);
    dolfin::Function k4((*gradY4)*max_stepsize);

    dolfin::Function y_next(*y + (k1*5/48.0) + (k2*27/56.0) + (k3*125/336.0) + (k4/24.0));

    l2ErrorCurrentEquation::Functional L2error_form (mesh);
    L2error_form.Function1 = y;
    L2error_form.Function2= y_next;
    double L2Error = sqrt(dolfin::assemble(L2error_form));

    return max_stepsize*std::pow(tol/L2Error, 0.25);

}

#endif //PROJECT_RUNGEKUTTA_H

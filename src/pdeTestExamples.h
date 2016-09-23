//
// Created by js on 03.07.16.
//

#ifndef DIFFUSION_FENICS_PDETESTEXAMPLES_H
#define DIFFUSION_FENICS_PDETESTEXAMPLES_H

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <dolfin.h>

namespace Poisson {
template <class Case> class SetupCase {
public:
  std::shared_ptr<dolfin::Mesh> getMesh() { return this->mesh; }
  /*    	std::shared_ptr<dolfin::Expression> getInitial()
              {
                      return this->initial;
              }
              std::shared_ptr<dolfin::SubDomain> getDirichletBoundary()
              {
                      return this->dirichletBoundary;
              }*/
  std::shared_ptr<dolfin::Expression> getSource() { return this->source; }
  std::shared_ptr<dolfin::Expression> getNeumann() { return this->neumann; }

  SetupCase(std::size_t dim)
      : source(std::make_shared<typename Case::Source>()),
        neumann(std::make_shared<typename Case::Neumann>()) {
    if (dim == 1) {
      mesh = std::make_shared<dolfin::UnitIntervalMesh>(50);
    } else if (dim == 2) {
      mesh = std::make_shared<dolfin::UnitSquareMesh>(50, 50);
    } else if (dim == 3) {
      mesh = std::make_shared<dolfin::UnitCubeMesh>(50, 50, 50);
    } else {
      mesh = std::make_shared<dolfin::UnitSquareMesh>(50, 50);
    }
  }

protected:
  std::shared_ptr<dolfin::Mesh> mesh;
  //      	std::shared_ptr<dolfin::Expression> initial;
  //        	std::shared_ptr<dolfin::SubDomain> dirichletBoundary;
  std::shared_ptr<dolfin::Expression> source;
  std::shared_ptr<dolfin::Expression> neumann;
};

class General {
public:
  class Initial : public dolfin::Expression {
    void eval(dolfin::Array<double> &values,
              const dolfin::Array<double> &x) const {
      values[0] = 0;
    }
  };

  // sourceTerm (right-hand side)
  class Source : public dolfin::Expression {
    void eval(dolfin::Array<double> &values,
              const dolfin::Array<double> &x) const {
      values[0] = 0;
      if (x[0] > 1 - (0.5 + DOLFIN_EPS) and x[0] < 1 - (0.5 - DOLFIN_EPS) and
          x[1] > 1 - (0.5 + DOLFIN_EPS) and x[1] < 1 - (0.5 - DOLFIN_EPS))
        values[0] = 20000;
    }
  };

  // Normal derivative (used for Neumann boundary condition)
  class Neumann : public dolfin::Expression {
    void eval(dolfin::Array<double> &values,
              const dolfin::Array<double> &x) const {
      values[0] = 0.0;
      ;
    }
  };
};
}

#endif // DIFFUSION_FENICS_PDETESTEXAMPLES_H

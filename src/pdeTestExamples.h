//
// Created by js on 03.07.16.
//

#ifndef DIFFUSION_FENICS_PDETESTEXAMPLES_H
#define DIFFUSION_FENICS_PDETESTEXAMPLES_H

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <dolfin.h>
#include <algorithm>

namespace Poisson {
    template <class Case> class SetupCase {
    public:
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
        }

    protected:
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
                values[0] = 1.0;
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

namespace ConvectionDiffusion {

    template<class Case>
    class SetupCase {
    public:
        std::shared_ptr <dolfin::Expression> getInitial() {
            return this->initial;
        }

        std::shared_ptr <dolfin::Expression> getSource() { return this->source; }

        std::shared_ptr <dolfin::Expression> getNeumann() {
            return this->neumann;
        }

        std::shared_ptr <dolfin::Expression> getVelocity() {
            return this->velocity;
        }

        std::shared_ptr <dolfin::Expression> getDiffusivity() {
            return this->diffusivity;
        }

        SetupCase(std::size_t dim)
                : initial(std::make_shared<typename Case::Initial>()),
                  source(std::make_shared<typename Case::Source>()),
                  neumann(std::make_shared<typename Case::Neumann>()),
                  velocity(std::make_shared<typename Case::Velocity>(dim)),
                  diffusivity(std::make_shared<typename Case::Diffusivity>()) {
        }

    protected:
        std::shared_ptr <dolfin::Expression> initial;
        std::shared_ptr <dolfin::Expression> source;
        std::shared_ptr <dolfin::Expression> neumann;
        std::shared_ptr <dolfin::Expression> velocity;
        std::shared_ptr <dolfin::Expression> diffusivity;
    };

    class General {

    public:
        class Initial : public dolfin::Expression {
            void eval(dolfin::Array<double> &values, const dolfin::Array<double> &x) const {
                values[0] = 0;
            }
        };

        // sourceTerm (right-hand side)
        class Source : public dolfin::Expression {
            void eval(dolfin::Array<double> &values, const dolfin::Array<double> &x) const {
                values[0] = 10;
            }
        };

        // Normal derivative (used for Neumann boundary condition)
        class Neumann : public dolfin::Expression {
            void eval(dolfin::Array<double> &values, const dolfin::Array<double> &x) const {
                values[0] = 0;
            }
        };


        class Velocity : public dolfin::Expression {
            void eval(dolfin::Array<double> &values, const dolfin::Array<double> &x) const {
                for (int i = 0; i < values.size(); i++) {
                    values[i] = 0;
                }
            }

        public:
            Velocity(std::size_t dim) : dolfin::Expression(dim) {};
        };

        class Diffusivity : public dolfin::Expression {
            void eval(dolfin::Array<double> &values, const dolfin::Array<double> &x) const {
                values[0] = 5;
            }
        };
    };
}

namespace Current{
    template <class Case> class SetupCase {
    public:

        std::shared_ptr<dolfin::Expression> getSource() { return this->source; }
        std::shared_ptr<dolfin::Expression> getNeumann() { return this->neumann; }
        std::shared_ptr<dolfin::Constant> getSigma1() { return this->sigma1; }
        std::shared_ptr<dolfin::Constant> getSigma2() { return this->sigma2; }

        SetupCase(std::size_t dim)
                : source(std::make_shared<typename Case::Source>()),
                  neumann(std::make_shared<typename Case::Neumann>()),
                  sigma1(std::make_shared<dolfin::Constant>(Case::sigma1)),
                  sigma2(std::make_shared<dolfin::Constan>(Case::sigma2))
        {
        }

    protected:
        std::shared_ptr<dolfin::Expression> source;
        std::shared_ptr<dolfin::Expression> neumann;
        std::shared_ptr<dolfin::Constant> sigma1;
        std::shared_ptr<dolfin::Constant> sigma2;
    };

    class General {
    public:

        // Normal derivative (used for Neumann boundary condition)
        class Neumann : public dolfin::Expression {
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &x) const {
                values[0] = 0.0;
            }
        };

        dolfin::Constant sigma1 = 84.935;
        dolfin::Constant sigma2 = 5.3e5;

        // sourceTerm (right-hand side)
        class Source : public dolfin::Expression {
        public:
            void setFunctions(dolfin::Function& u1, dolfin::Function& u2) {
                this->u1 = u1;
                this->u2 = u2;
            }

        private:
            dolfin::Function u1;
            dolfin::Function u2; // würde ich als const Referenz/ Pointer/ Pointer realisieren
            double alpha = 0.3;
            double n = 2.0;
            double j = 0.005;
            double R = 8.3144598; //allgemeine Gaskonstante
            double T = 293.15; //20°C in Kelvin
            double F = 96485.309; //Farady-Konstante

            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &x) const {
                dolfin::Array<double> values1;
                dolfin::Array<double> values2;
                u1.eval(values,x);
                u2.eval(values,x);
                dolfin::Array<double> diff;
#ifndef NDEBUG
                assert(values1.size()== values2.size());
                assert(values1.size()== 1);
#endif
                for(int i = 0; i < values1.size(); i++){
                    diff[i] = values2[i]-values1[i]-1.7;
                    values[i] = j*(std::exp(alpha*n*F*(diff[i])/(R*T)) - std::exp(-(1.0-alpha)*n*F*(diff[i])/(R*T)));
                }
            }
        };
    };
}
#endif // DIFFUSION_FENICS_PDETESTEXAMPLES_H

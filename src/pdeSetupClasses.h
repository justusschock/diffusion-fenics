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
        SetupCase(std::size_t dim, std::string meshFile, std::string subdomainFile, std::string facetFile)
                : mesh(std::make_shared<dolfin::Mesh>(meshFile)),
                  subdomain_function(std::make_shared<dolfin::MeshFunction<size_t>>(mesh, subdomainFile)),
                  facet_function(std::make_shared<dolfin::MeshFunction<size_t>>(mesh, facetFile)),
                  source1(std::make_shared<typename Case::Source1>(u1,u2)),
                  source2(std::make_shared<typename Case::Source2>(u1, u2)),
                  sigma1(std::make_shared<typename Case::Sigma1>()),
                  sigma2(std::make_shared<typename Case::Sigma2>())

        {

        }
        std::shared_ptr<typename Case::Source1> getSource1() { return this->source1; }
        std::shared_ptr<typename Case::Source2> getSource2() { return this->source2; }
        std::shared_ptr<typename Case::Sigma1> getSigma1() { return this->sigma1; }
        std::shared_ptr<typename Case::Sigma2> getSigma2() { return this->sigma2; }
        std::shared_ptr<dolfin::Mesh> getMesh() { return this->mesh; }
        std::shared_ptr<dolfin::MeshFunction<size_t>> getSubdomainFunction() { return this->subdomain_function; }
        std::shared_ptr<dolfin::MeshFunction<size_t>> getFacetFunction() { return this->facet_function; }
        std::shared_ptr<dolfin::Function> getU1() { return this->u1; }
        std::shared_ptr<dolfin::Function> getU2() { return this->u2; }

        void setU1(dolfin::Function* u) { u1.reset(u); }
        void setU2(dolfin::Function* u) { u2.reset(u); }

    protected:
        std::shared_ptr<typename Case::Source1> source1;
        std::shared_ptr<typename Case::Source2> source2;
        std::shared_ptr<typename Case::Sigma1> sigma1;
        std::shared_ptr<typename Case::Sigma2> sigma2;
        std::shared_ptr<dolfin::Mesh> mesh;
        std::shared_ptr<dolfin::MeshFunction<size_t>> subdomain_function;
        std::shared_ptr<dolfin::MeshFunction<size_t>> facet_function;

        std::shared_ptr<dolfin::Function> u1;
        std::shared_ptr<dolfin::Function> u2;
    };

    class General {
    public:

        static constexpr double R = 8.3144598; //allgemeine Gaskonstante
        static constexpr double T = 293.15; //20°C in Kelvin
        static constexpr double F = 96485.309; //Farady-Konstante

        // sourceTerms (right-hand side)
        class Source1: public dolfin::Expression {
        public:
            Source1(std::shared_ptr<dolfin::Function> u1, std::shared_ptr<dolfin::Function> u2)
                    :dolfin::Expression(),u1(u1), u2(u2){ }

            virtual ~Source1(){
                u1.reset();
                u2.reset();
            }

        private:
            double alpha = 0.3;
            double n = 2.0;
            double j = 0.005;
            std::shared_ptr<dolfin::Function> u1;
            std::shared_ptr<dolfin::Function> u2;

            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &x) const {
                dolfin::Array<double> values1(values.size());
                dolfin::Array<double> values2(values.size());
                u1->eval(values,x);
                u2->eval(values,x);
                dolfin::Array<double> diff(values.size());

                for(int i = 0; i < values1.size(); i++){
                    diff[i] = values2[i]-values1[i]-1.7;
                    values[i] = j*(std::exp(alpha*n*F*(diff[i])/(R*T)) - std::exp(-(1.0-alpha)*n*F*(diff[i])/(R*T)));
                }
            }
        };

        class Source2 : public dolfin::Expression {

        public:
            Source2(std::shared_ptr<dolfin::Function> u1, std::shared_ptr<dolfin::Function> u2)
                    :dolfin::Expression(),u1(u1), u2(u2){ }

            virtual ~Source2(){
                u1.reset();
                u2.reset();
            }

        private:
            double alpha = 0.3;
            double n = 2.0;
            double j = 0.005;
            std::shared_ptr<dolfin::Function> u1;
            std::shared_ptr<dolfin::Function> u2;

            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &x) const {
                dolfin::Array<double> values1(values.size());
                dolfin::Array<double> values2(values.size());
                u1->eval(values,x);
                u2->eval(values,x);
                dolfin::Array<double> diff(values.size());

                for(int i = 0; i < values1.size(); i++){
                    diff[i] = values1[i]-values2[i]-0.4;
                    values[i] = j*(std::exp(alpha*n*F*(diff[i])/(R*T)) - std::exp(-(1.0-alpha)*n*F*(diff[i])/(R*T)));
                }
            }
        };

        class Sigma1 : public dolfin::Expression {

        private:
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &x) const {
                values[0] = 84.935;
            }
        };

        class Sigma2 : public dolfin::Expression {

        private:
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &x) const {
                values[0] = 5.3e5;
            }
        };

    };
}
#endif // DIFFUSION_FENICS_PDETESTEXAMPLES_H

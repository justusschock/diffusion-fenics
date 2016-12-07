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
        SetupCase(std::size_t dim, std::string meshFile, std::string subdomainFile, std::string facetFile)
        : mesh(std::make_shared<dolfin::Mesh>(meshFile)),
          subdomain_function(std::make_shared<dolfin::MeshFunction<size_t>>(mesh, subdomainFile)),
          facet_function(std::make_shared<dolfin::MeshFunction<size_t>>(mesh, facetFile)),
          source(std::make_shared<typename Case::Source>()),
          neumann(std::make_shared<typename Case::Neumann>()),
          initial(std::make_shared<typename Case::Initial>())
        {

        }

        virtual ~SetupCase(){
            source.reset();
            neumann.reset();
            initial.reset();
            mesh.reset();
            facet_function.reset();
            subdomain_function.reset();
            u.reset();
        }

        std::shared_ptr<dolfin::Mesh> getMesh() { return mesh;}
        std::shared_ptr<typename Case::Initial> getInitial() { return initial; }
        std::shared_ptr<typename Case::Source> getSource() { return source; }
        std::shared_ptr<typename Case::Neumann> getNeumann() { return neumann; }
        std::shared_ptr<dolfin::MeshFunction<size_t>> getSubDomainFunction() { return subdomain_function; }
        std::shared_ptr<dolfin::MeshFunction<size_t>> getFacetFunction() { return facet_function; }
        std::shared_ptr<dolfin::Function> getU(){ return u; }

        void setU(dolfin::Function* u_) { u.reset(u_); }

    protected:
        //      	std::shared_ptr<dolfin::Expression> initial;
        //        	std::shared_ptr<dolfin::SubDomain> dirichletBoundary;
        std::shared_ptr<typename Case::Source> source;
        std::shared_ptr<typename Case::Neumann> neumann;
        std::shared_ptr<typename Case::Initial> initial;
        std::shared_ptr<dolfin::Mesh> mesh;
        std::shared_ptr<dolfin::MeshFunction<size_t>> facet_function;
        std::shared_ptr<dolfin::MeshFunction<size_t>> subdomain_function;

        std::shared_ptr<dolfin::Function>u;
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

    template <class Case> class SetupCase {
    public:

        SetupCase(std::size_t dim, std::string meshFile, std::string subdomainFile, std::string facetFile,
                  dolfin::Constant dirichletValue = dolfin::Constant(0.0))
                : mesh(std::make_shared<dolfin::Mesh>(meshFile)),
                  subdomain_function(std::make_shared<dolfin::MeshFunction<size_t>>(mesh, subdomainFile)),
                  facet_function(std::make_shared<dolfin::MeshFunction<size_t>>(mesh, facetFile)),
                  initial(std::make_shared<typename Case::Initial>()),
                  velocity(std::make_shared<typename Case::Velocity>(dim)),
                  source(std::make_shared<typename Case::Source>()),
                  neumann(std::make_shared<typename Case::Neumann>()),
                  diffusivity(std::make_shared<typename Case::Diffusivity>()),
                  dirichletValue(std::make_shared<dolfin::Constant>(dirichletValue)){ }

        virtual ~SetupCase()
        {
            mesh.reset();
            dirichletValue.reset();
            initial.reset();
            velocity.reset();
            source.reset();
            neumann.reset();
            diffusivity.reset();
            subdomain_function.reset();
            facet_function.reset();
            u.reset();
        }

        std::shared_ptr<dolfin::Mesh> getMesh() { return mesh;}
        std::shared_ptr<dolfin::Constant> getDirichletValue() { return dirichletValue; }
        std::shared_ptr<typename Case::Initial> getInitial() { return initial; }
        std::shared_ptr<typename Case::Velocity> getVelocity() { return velocity; }
        std::shared_ptr<typename Case::Source> getSource() { return source; }
        std::shared_ptr<typename Case::Neumann> getNeumann() { return neumann; }
        std::shared_ptr<typename Case::Diffusivity> getDiffusivity() { return diffusivity; }
        std::shared_ptr<dolfin::MeshFunction<size_t>> getSubDomainFunction() { return subdomain_function; }
        std::shared_ptr<dolfin::MeshFunction<size_t>> getFacetFunction() { return facet_function; }
        std::shared_ptr<dolfin::Function> getU(){ return u; }

        void setU(dolfin::Function* u_) { u.reset(u_); }

    protected:
        std::shared_ptr<dolfin::Mesh> mesh;
        std::shared_ptr<dolfin::Constant> dirichletValue;
        std::shared_ptr<typename Case::Initial> initial;
        std::shared_ptr<typename Case::Velocity> velocity;
        std::shared_ptr<typename Case::Source> source;
        std::shared_ptr<typename Case::Neumann> neumann;
        std::shared_ptr<typename Case::Diffusivity> diffusivity;
        std::shared_ptr<dolfin::MeshFunction<size_t>> subdomain_function;
        std::shared_ptr<dolfin::MeshFunction<size_t>> facet_function;

        std::shared_ptr<dolfin::Function> u;
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
                  source_l(std::make_shared<typename Case::Source_L>(u_L,u_S)),
                  source_s(std::make_shared<typename Case::Source_S>(u_L, u_S)),
                  sigma_l(std::make_shared<typename Case::Sigma_L>()),
                  sigma_s(std::make_shared<typename Case::Sigma_S>())

        { }

        virtual ~SetupCase(){
            source_l.reset();
            source_s.reset();
            sigma_l.reset();
            sigma_s.reset();
            mesh.reset();
            subdomain_function.reset();
            facet_function.reset();
            u_L.reset();
            u_S.reset();
        }

        std::shared_ptr<typename Case::Source_L> getSourceL() { return this->source_l; }
        std::shared_ptr<typename Case::Source_S> getSourceS() { return this->source_s; }
        std::shared_ptr<typename Case::Sigma_L> getSigmaL() { return this->sigma_l; }
        std::shared_ptr<typename Case::Sigma_S> getSigmaS() { return this->sigma_s; }
        std::shared_ptr<dolfin::Mesh> getMesh() { return this->mesh; }
        std::shared_ptr<dolfin::MeshFunction<size_t>> getSubdomainFunction() { return this->subdomain_function; }
        std::shared_ptr<dolfin::MeshFunction<size_t>> getFacetFunction() { return this->facet_function; }
        std::shared_ptr<dolfin::Function> getUL() { return this->u_L; }
        std::shared_ptr<dolfin::Function> getUS() { return this->u_S; }

        void setUS(dolfin::Function* u)
        {
            u_S.reset(u);
            source_s->setUS(u_S);
        }
        void setUL(dolfin::Function* u)
        {
            u_L.reset(u);
            source_l->setUL(u_L);
        }

    protected:
        std::shared_ptr<typename Case::Source_L> source_l;
        std::shared_ptr<typename Case::Source_S> source_s;
        std::shared_ptr<typename Case::Sigma_L> sigma_l;
        std::shared_ptr<typename Case::Sigma_S> sigma_s;
        std::shared_ptr<dolfin::Mesh> mesh;
        std::shared_ptr<dolfin::MeshFunction<size_t>> subdomain_function;
        std::shared_ptr<dolfin::MeshFunction<size_t>> facet_function;

        std::shared_ptr<dolfin::Function> u_L;
        std::shared_ptr<dolfin::Function> u_S;
    };

    class General {
    public:

        static constexpr double R = 8.3144598; //allgemeine Gaskonstante
        static constexpr double T = 293.15; //20°C in Kelvin
        static constexpr double F = 96485.309; //Farady-Konstante

        // sourceTerms (right-hand side)
        class Source_L: public dolfin::Expression {
        public:
            Source_L(std::shared_ptr<dolfin::Function> u_L, std::shared_ptr<dolfin::Function> u_S)
                    :dolfin::Expression(),u_L(u_L), u_S(u_S){ }

            virtual ~Source_L(){
                u_L.reset();
                u_S.reset();
            }

            void setUL(std::shared_ptr<dolfin::Function> U_L)
            {
                u_L.reset();
                u_L = U_L;
            }

            void setUS(std::shared_ptr<dolfin::Function> U_S)
            {
                u_S.reset();
                u_S = U_S;
            }

        private:
            double alpha = 0.3;
            double n = 2.0;
            double j = 0.005;
            std::shared_ptr<dolfin::Function> u_L;
            std::shared_ptr<dolfin::Function> u_S;

            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &x) const {
                dolfin::Array<double> values_L(values.size());
                dolfin::Array<double> values_S(values.size());
                u_L->eval(values_L,x);
                u_S->eval(values_S,x);
                dolfin::Array<double> diff(values.size());

                for(int i = 0; i < values.size(); i++){
                    diff[i] = values_L[i]-values_S[i]-1.7;
                    values[i] = j*(std::exp(alpha*n*F*(diff[i])/(R*T)) - std::exp(-(1.0-alpha)*n*F*(diff[i])/(R*T)));
                }
            }
        };

        class Source_S : public dolfin::Expression {

        public:
            Source_S(std::shared_ptr<dolfin::Function> u_L, std::shared_ptr<dolfin::Function> u_S)
                    :dolfin::Expression(),u_L(u_L), u_S(u_S){ }

            virtual ~Source_S(){
                u_L.reset();
                u_S.reset();
            }

            void setUL(std::shared_ptr<dolfin::Function> U_L)
            {
                u_L.reset();
                u_L = U_L;
            }

            void setUS(std::shared_ptr<dolfin::Function> U_S)
            {
                u_S.reset();
                u_S = U_S;
            }

        private:
            double alpha = 0.3;
            double n = 2.0;
            double j = 0.005;
            std::shared_ptr<dolfin::Function> u_L;
            std::shared_ptr<dolfin::Function> u_S;

            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &x) const {
                dolfin::Array<double> values_S(values.size());
                dolfin::Array<double> values_L(values.size());
                u_L->eval(values_L,x);
                u_S->eval(values_S,x);
                dolfin::Array<double> diff(values.size());

                for(int i = 0; i < values.size(); i++){
                    diff[i] = values_S[i]-values_L[i]-0.4;
                    values[i] = j*(std::exp(alpha*n*F*(diff[i])/(R*T)) - std::exp(-(1.0-alpha)*n*F*(diff[i])/(R*T)));
                }
            }
        };

        class Sigma_L : public dolfin::Expression {

        private:
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &x) const {
                values[0] = 84.935;
            }
        };

        class Sigma_S : public dolfin::Expression {

        private:
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &x) const {
                values[0] = 5.3e5;
            }
        };

    };
}
#endif // DIFFUSION_FENICS_PDETESTEXAMPLES_H

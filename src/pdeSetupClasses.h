//
// Created by js on 03.07.16.
//

#ifndef DIFFUSION_FENICS_PDETESTEXAMPLES_H
#define DIFFUSION_FENICS_PDETESTEXAMPLES_H

#include <dolfin.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>

namespace Poisson {
    template <class Case>
    class SetupCase {
       public:
        SetupCase(std::size_t dim,
                  std::string meshFile,
                  std::string subdomainFile,
                  std::string facetFile)
            : mesh(std::make_shared<dolfin::Mesh>(meshFile)),
              subdomain_function(std::make_shared<dolfin::MeshFunction<size_t>>(
                  mesh, subdomainFile)),
              facet_function(std::make_shared<dolfin::MeshFunction<size_t>>(
                  mesh, facetFile)),
              source(std::make_shared<typename Case::Source>()),
              neumann(std::make_shared<typename Case::Neumann>()),
              initial(std::make_shared<typename Case::Initial>())
        {
        }

        virtual ~SetupCase()
        {
            source.reset();
            neumann.reset();
            initial.reset();
            mesh.reset();
            facet_function.reset();
            subdomain_function.reset();
            u.reset();
        }

        std::shared_ptr<dolfin::Mesh> getMesh() { return mesh; }
        std::shared_ptr<typename Case::Initial> getInitial() { return initial; }
        std::shared_ptr<typename Case::Source> getSource() { return source; }
        std::shared_ptr<typename Case::Neumann> getNeumann() { return neumann; }
        std::shared_ptr<dolfin::MeshFunction<size_t>> getSubDomainFunction()
        {
            return subdomain_function;
        }
        std::shared_ptr<dolfin::MeshFunction<size_t>> getFacetFunction()
        {
            return facet_function;
        }
        std::shared_ptr<dolfin::Function> getU() { return u; }
        void setU(dolfin::Function *u_) { u.reset(u_); }
       protected:
        //      	std::shared_ptr<dolfin::Expression> initial;
        //        	std::shared_ptr<dolfin::SubDomain> dirichletBoundary;
        std::shared_ptr<typename Case::Source> source;
        std::shared_ptr<typename Case::Neumann> neumann;
        std::shared_ptr<typename Case::Initial> initial;
        std::shared_ptr<dolfin::Mesh> mesh;
        std::shared_ptr<dolfin::MeshFunction<size_t>> facet_function;
        std::shared_ptr<dolfin::MeshFunction<size_t>> subdomain_function;

        std::shared_ptr<dolfin::Function> u;
    };

    class General {
       public:
        class Initial : public dolfin::Expression {
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = 0;
            }
        };

        // sourceTerm (right-hand side)
        class Source : public dolfin::Expression {
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = 1.0;
            }
        };

        // Normal derivative (used for Neumann boundary condition)
        class Neumann : public dolfin::Expression {
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = 0.0;
                ;
            }
        };
    };
}

namespace ConvectionDiffusion {

    template <class Case>
    class SetupCase {
       public:
        SetupCase(std::size_t dim,
                  std::string meshFile,
                  std::string subdomainFile,
                  std::string facetFile,
                  dolfin::Constant dirichletValue = dolfin::Constant(0.0))
            : mesh(std::make_shared<dolfin::Mesh>(meshFile)),
              dirichletValue(
                  std::make_shared<dolfin::Constant>(dirichletValue)),
              initial(std::make_shared<typename Case::Initial>()),
              velocity(std::make_shared<typename Case::Velocity>(dim)),
              source(std::make_shared<typename Case::Source>()),
              neumann(std::make_shared<typename Case::Neumann>()),
              diffusivity(std::make_shared<typename Case::Diffusivity>()),
              subdomain_function(std::make_shared<dolfin::MeshFunction<size_t>>(
                  mesh, subdomainFile)),
              facet_function(std::make_shared<dolfin::MeshFunction<size_t>>(
                  mesh, facetFile))
        {
        }

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

        std::shared_ptr<dolfin::Mesh> getMesh() { return mesh; }
        std::shared_ptr<dolfin::Constant> getDirichletValue()
        {
            return dirichletValue;
        }
        std::shared_ptr<typename Case::Initial> getInitial() { return initial; }
        std::shared_ptr<typename Case::Velocity> getVelocity()
        {
            return velocity;
        }
        std::shared_ptr<typename Case::Source> getSource() { return source; }
        std::shared_ptr<typename Case::Neumann> getNeumann() { return neumann; }
        std::shared_ptr<typename Case::Diffusivity> getDiffusivity()
        {
            return diffusivity;
        }
        std::shared_ptr<dolfin::MeshFunction<size_t>> getSubDomainFunction()
        {
            return subdomain_function;
        }
        std::shared_ptr<dolfin::MeshFunction<size_t>> getFacetFunction()
        {
            return facet_function;
        }
        std::shared_ptr<dolfin::Function> getU() { return u; }
        void setU(dolfin::Function *u_) { u.reset(u_); }
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
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = 0;
            }
        };

        // sourceTerm (right-hand side)
        class Source : public dolfin::Expression {
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = 10;
            }
        };

        // Normal derivative (used for Neumann boundary condition)
        class Neumann : public dolfin::Expression {
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = 0;
            }
        };

        class Velocity : public dolfin::Expression {
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                for (std::size_t i = 0; i < values.size(); i++) {
                    values[i] = 0;
                }
            }

           public:
            Velocity(std::size_t dim) : dolfin::Expression(dim){};
        };

        class Diffusivity : public dolfin::Expression {
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = 5;
            }
        };
    };
}

namespace Current {
    template <class Case>
    class SetupCase {
       public:
        SetupCase(std::size_t,
                  std::string meshFile,
                  std::string subdomainFile,
                  std::string facetFile)
            : mesh(std::make_shared<dolfin::Mesh>(meshFile)),
              subdomain_function(std::make_shared<dolfin::MeshFunction<size_t>>(
                  mesh, subdomainFile)),
              facet_function(std::make_shared<dolfin::MeshFunction<size_t>>(
                  mesh, facetFile)),
              source_(std::make_shared<typename Case::Source>(u_)),
              sigma_(std::make_shared<typename Case::Sigma>())
        {
        }

        virtual ~SetupCase()
        {
            source_.reset();
            sigma_.reset();
            mesh.reset();
            subdomain_function.reset();
            facet_function.reset();
            u_.reset();
        }

        std::shared_ptr<typename Case::Source> getSource()
        {
            return this->source_;
        }

        std::shared_ptr<typename Case::Sigma> getSigma()
        {
            return this->sigma_;
        }
        std::shared_ptr<dolfin::Mesh> getMesh() { return this->mesh; }
        std::shared_ptr<dolfin::MeshFunction<size_t>> getSubdomainFunction()
        {
            return this->subdomain_function;
        }
        std::shared_ptr<dolfin::MeshFunction<size_t>> getFacetFunction()
        {
            return this->facet_function;
        }
        std::shared_ptr<dolfin::Function> getU() { return this->u_; }

        void setU(std::shared_ptr<dolfin::Function> u)
        {
            u_.reset();
            u_ = u;
            source_->setU(u_);
        }

       protected:
        std::shared_ptr<dolfin::Mesh> mesh;
        std::shared_ptr<dolfin::MeshFunction<size_t>> subdomain_function;
        std::shared_ptr<dolfin::MeshFunction<size_t>> facet_function;
        std::shared_ptr<dolfin::Function> u_;
        std::shared_ptr<typename Case::Source> source_;
        std::shared_ptr<typename Case::Sigma> sigma_;
    };

    class General {
       public:
        static constexpr double R = 8.3144598;  // allgemeine Gaskonstante
        static constexpr double T = 293.15;     // 20°C in Kelvin
        static constexpr double F = 96485.309;  // Farady-Konstante

        // sourceTerms (right-hand side)
        /*class Source_L : public dolfin::Expression {
           public:
            Source_L(std::shared_ptr<dolfin::Function> u_L,
                     std::shared_ptr<dolfin::Function> u_S)
                : dolfin::Expression(), u_L(u_L), u_S(u_S)
            {
            }

            virtual ~Source_L()
            {
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
                      const dolfin::Array<double> &x) const
            {
                dolfin::Array<double> values_L(values.size());
                dolfin::Array<double> values_S(values.size());
                u_L->eval(values_L, x);
                u_S->eval(values_S, x);
                dolfin::Array<double> diff(values.size());

                for (std::size_t i = 0; i < values.size(); i++) {
                    diff[i] = values_S[i] - values_L[i] - 1.7;
                    values[i] =
                        -j * (std::exp(alpha * n * F * (diff[i]) / (R * T)) -
                              std::exp(-(1.0 - alpha) * n * F * (diff[i]) /
                                       (R * T)));
                }
            }
        };

        class Source_S : public dolfin::Expression {
           public:
            Source_S(std::shared_ptr<dolfin::Function> u_L,
                     std::shared_ptr<dolfin::Function> u_S)
                : dolfin::Expression(), u_L(u_L), u_S(u_S)
            {
            }

            virtual ~Source_S()
            {
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
                      const dolfin::Array<double> &x) const
            {
                dolfin::Array<double> values_S(values.size());
                dolfin::Array<double> values_L(values.size());
                u_L->eval(values_L, x);
                u_S->eval(values_S, x);
                dolfin::Array<double> diff(values.size());

                for (std::size_t i = 0; i < values.size(); i++) {
                    diff[i] = values_S[i] - values_L[i] - 1.7;
                    values[i] =
                        j * (std::exp(alpha * n * F * (diff[i]) / (R * T)) -
                             std::exp(-(1.0 - alpha) * n * F * (diff[i]) /
                                      (R * T)));
                }
            }
        };*/

        class Source : public dolfin::Expression {
        public:
            Source(std::shared_ptr<dolfin::Function> u)
                    : dolfin::Expression(), u_(u)
            {
            }

            virtual ~Source()
            {
                u_.reset();
            }

            void setU(std::shared_ptr<dolfin::Function> U_)
            {
                u_.reset();
                u_ = U_;
            }


        private:
            double alpha = 0.3;
            double n = 2.0;
            double j = 0.005;
            std::shared_ptr<dolfin::Function> u_;

            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &x) const
            {
                dolfin::Array<double> values_(values.size());
                u_->eval(values_, x);
                dolfin::Array<double> diff(values.size());


                    diff[0] = values_[0] - values_[1] - 1.7;
                    values[0] = -j * (std::exp(alpha * n * F * (diff[0]) / (R * T)) -
                                 std::exp(-(1.0 - alpha) * n * F * (diff[0]) /
                                          (R * T)));
                    values[1] = -values[0];

            }
        };

        /*class Sigma_L : public dolfin::Expression {
           private:
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = 84.935;
            }
        };

        class Sigma_S : public dolfin::Expression {
           private:
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = 5.3e5;
            }
        }; */
        class Sigma : public dolfin::Expression {
        private:
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = 84.935; //Liquid
                values[1] = 5.3e5; //Solid
            }
        };
    };
}
#endif  // DIFFUSION_FENICS_PDETESTEXAMPLES_H

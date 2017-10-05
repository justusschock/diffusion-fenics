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
                values[0] = 1;
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
              sigma_(std::make_shared<typename Case::Sigma>()),
              RT_(std::make_shared<typename Case::RT>()),
              U_eq_(std::make_shared<typename Case::U_eq>()),
              Alpha_(std::make_shared<typename Case::Alpha>()),
              I0_(std::make_shared<typename Case::I0>()),
              neumann_(std::make_shared<typename Case::Neumann>()),
              t(0.0)
        {
        }

        virtual ~SetupCase() {}
        std::shared_ptr<typename Case::Sigma> getSigma()
        {
            return this->sigma_;
        }
        std::shared_ptr<typename Case::RT> getRT() { return this->RT_; }
        std::shared_ptr<typename Case::U_eq> getU_eq() { return this->U_eq_; }
        std::shared_ptr<typename Case::I0> getI0() { return this->I0_; }
        std::shared_ptr<typename Case::Alpha> getAlpha()
        {
            return this->Alpha_;
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
        std::shared_ptr<typename Case::Neumann> getNeumann()
        {
            return this->neumann_;
        }

        void setU(std::shared_ptr<dolfin::Function> u)
        {
            u_.reset();
            u_ = u;
        }

        void setTime(double time)
        {
            t = time;
            Alpha_->setTime(time);
            I0_->setTime(time);
            RT_->setTime(time);
        }

       protected:
        std::shared_ptr<dolfin::Mesh> mesh;
        std::shared_ptr<dolfin::MeshFunction<size_t>> subdomain_function;
        std::shared_ptr<dolfin::MeshFunction<size_t>> facet_function;
        std::shared_ptr<dolfin::Function> u_;
        std::shared_ptr<typename Case::Sigma> sigma_;
        std::shared_ptr<typename Case::RT> RT_;
        std::shared_ptr<typename Case::U_eq> U_eq_;
        std::shared_ptr<typename Case::Alpha> Alpha_;
        std::shared_ptr<typename Case::I0> I0_;
        std::shared_ptr<typename Case::Neumann> neumann_;
        double t;
    };

    class FixPotentialDifference {
       public:
        static constexpr double R = 8.3144598;  // allgemeine Gaskonstante
        static constexpr double T = 293.15;     // 20°C in Kelvin
        static constexpr double F = 96485.309;  // Farady-Konstante

        class Sigma : public dolfin::Expression {
           public:
            Sigma() : dolfin::Expression(3) {}
           private:
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = 5.3e5;   // Solid Pos
                values[1] = 84.935;  // Liquid
                values[2] = 1e6;     // Solid Neg
            }
        };
        class U_eq : public dolfin::Expression {
           public:
            U_eq() : dolfin::Expression() {}
           private:
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = 1.;
            }
        };
        class Alpha : public dolfin::Expression {
           public:
            Alpha() : dolfin::Expression() {}
            void setTime(double t) { time = t; }
           private:
            double time;

            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = 0.5;
            }
        };
        class I0 : public dolfin::Expression {
           public:
            I0() : dolfin::Expression() {}
            void setTime(double t) { time = t; }
           private:
            double time;

            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = 1e-2;
            }
        };
        class RT : public dolfin::Expression {
           public:
            RT() : dolfin::Expression() {}
            void setTime(double t) { time = t; }
           private:
            double time;

            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = R * T / (2.0 * F);
            }
        };

        class Neumann : public dolfin::Expression {
           public:
            Neumann() : dolfin::Expression() {}
           private:
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = 0;
            }
        };
    };

    class FixCurrent {
       public:
        static constexpr double R = 8.3144598;  // allgemeine Gaskonstante
        static constexpr double T = 293.15;     // 20°C in Kelvin
        static constexpr double F = 96485.309;  // Farady-Konstante

        class Sigma : public dolfin::Expression {
           public:
            Sigma() : dolfin::Expression(3) {}
           private:
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = 5.3e5;   // Solid Pos
                values[1] = 84.935;  // Liquid
                values[2] = 1e6;     // Solid Neg
            }
        };
        class U_eq : public dolfin::Expression {
           public:
            U_eq() : dolfin::Expression() {}
           private:
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = 1.;
            }
        };
        class Alpha : public dolfin::Expression {
           public:
            Alpha() : dolfin::Expression() {}
            void setTime(double t) { time = t; }
           private:
            double time;

            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = 0.5;
            }
        };
        class I0 : public dolfin::Expression {
           public:
            I0() : dolfin::Expression() {}
            void setTime(double t) { time = t; }
           private:
            double time;

            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = 1e-2;
            }
        };
        class RT : public dolfin::Expression {
           public:
            RT() : dolfin::Expression() {}
            void setTime(double t) { time = t; }
           private:
            double time;

            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = R * T / (2.0 * F);
            }
        };

        class Neumann : public dolfin::Expression {
           public:
            Neumann() : dolfin::Expression() {}
           private:
            void eval(dolfin::Array<double> &values,
                      const dolfin::Array<double> &) const
            {
                values[0] = 1;
            }
        };
    };
}
#endif  // DIFFUSION_FENICS_PDETESTEXAMPLES_H

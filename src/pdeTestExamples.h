//
// Created by js on 03.07.16.
//

#ifndef DIFFUSION_FENICS_PDETESTEXAMPLES_H
#define DIFFUSION_FENICS_PDETESTEXAMPLES_H

#include <dolfin.h>
#include <cstdlib>
#include <ctime>
#include <cmath>

namespace Poisson {
	template <class Case>
	class SetupCase {
		public:
        	std::shared_ptr<dolfin::Mesh> getMesh() { return this->mesh; }
        	std::shared_ptr<dolfin::Expression> getInitial()
        	{
            		return this->initial;
        	}
        	std::shared_ptr<dolfin::SubDomain> getDirichletBoundary()
        	{
            		return this->dirichletBoundary;
        	}
        	std::shared_ptr<dolfin::Expression> getSource() { return this->source; }
        	std::shared_ptr<dolfin::Expression> getNeumann()
        	{
            		return this->neumann;
        	}
        

		SetupCase(std::size_t dim)
            	: initial(std::make_shared<typename Case::Initial>()),
              	dirichletBoundary(
                  std::make_shared<typename Case::DirichletBoundary>()),
              	source(std::make_shared<typename Case::Source>()),
              	neumann(std::make_shared<typename Case::Neumann>())
        	{
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
        	std::shared_ptr<dolfin::Expression> initial;
        	std::shared_ptr<dolfin::SubDomain> dirichletBoundary;
        	std::shared_ptr<dolfin::Expression> source;
        	std::shared_ptr<dolfin::Expression> neumann;
	};

	class General {
		public:

		class Initial : public dolfin::Expression {
    			void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    			{
        			values[0] = 0;
    			}
		};

		// sourceTerm (right-hand side)
		class Source : public dolfin::Expression {
    			void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    			{
        			values[0] = 0;
        			if (x[0] > 1 - (0.5 + DOLFIN_EPS) and x[0] < 1 - (0.5 - DOLFIN_EPS) and
            			x[1] > 1 - (0.5 + DOLFIN_EPS) and x[1] < 1 - (0.5 - DOLFIN_EPS))
            			values[0] = 20000;
    		}
		};

		// Normal derivative (used for Neumann boundary condition)
		class Neumann : public dolfin::Expression {
    			void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    			{
        			values[0] = x[0];
    			}
		};

		// Sub domain for Dirichlet boundary condition
		class DirichletBoundary : public dolfin::SubDomain {
    			bool inside(const dolfin::Array<double>& x, bool on_boundary) const
    			{
        			return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS or
               			x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS;
    			}
		};

	};

}

namespace ConvectionDiffusion {

	template <class Case>
    	class SetupCase {
       		public:
        	std::shared_ptr<dolfin::Mesh> getMesh() { return this->mesh; }
        	std::shared_ptr<dolfin::Expression> getInitial()
        	{
            		return this->initial;
        	}
        	std::shared_ptr<dolfin::SubDomain> getDirichletBoundary()
        	{
            		return this->dirichletBoundary;
        	}
        	std::shared_ptr<dolfin::Expression> getSource() { return this->source; }
        	std::shared_ptr<dolfin::Expression> getNeumann()
        	{
            		return this->neumann;
        	}
        	std::shared_ptr<dolfin::Expression> getVelocity()
        	{
            		return this->velocity;
        	}
        	std::shared_ptr<dolfin::Expression> getDiffusivity()
        	{
            		return this->diffusivity;
        	}

		SetupCase(std::size_t dim)
            	: initial(std::make_shared<typename Case::Initial>()),
              	dirichletBoundary(
                  std::make_shared<typename Case::DirichletBoundary>()),
              	source(std::make_shared<typename Case::Source>()),
              	neumann(std::make_shared<typename Case::Neumann>()),
              	velocity(std::make_shared<typename Case::Velocity>(dim)),
              	diffusivity(std::make_shared<typename Case::Diffusivity>())
        	{
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
        	std::shared_ptr<dolfin::Expression> initial;
        	std::shared_ptr<dolfin::SubDomain> dirichletBoundary;
        	std::shared_ptr<dolfin::Expression> source;
        	std::shared_ptr<dolfin::Expression> neumann;
        	std::shared_ptr<dolfin::Expression> velocity;
        	std::shared_ptr<dolfin::Expression> diffusivity;
	};

    	class General {

		public:
		class Initial : public dolfin::Expression {
    			void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    			{
        			values[0] = 0;
    			}
		};

		// sourceTerm (right-hand side)
		class Source : public dolfin::Expression {
    			void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    			{
        			values[0] = 0;
        			if (x[0] > 1 - (0.5 + DOLFIN_EPS) and x[0] < 1 - (0.5 - DOLFIN_EPS) and
            			x[1] > 1 - (0.5 + DOLFIN_EPS) and x[1] < 1 - (0.5 - DOLFIN_EPS))
            				values[0] = 20000;
    			}
		};

		// Normal derivative (used for Neumann boundary condition)
		class Neumann : public dolfin::Expression {
    			void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    			{
        			values[0] = x[0];
 		   	}
		};

		// Sub domain for Dirichlet boundary condition
		class DirichletBoundary : public dolfin::SubDomain {
    			bool inside(const dolfin::Array<double>& x, bool on_boundary) const
    			{
        			return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS or
               			x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS;
    			}
		};

		class Velocity : public dolfin::Expression {
    			void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    			{
        			for (int i = 0; i < values.size(); i++) {
            				values[i] = 0;
        			}
    			}

   			public:
    			Velocity(std::size_t dim) : dolfin::Expression(dim){};
		};

		class Diffusivity : public dolfin::Expression {
    			void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    			{
        			values[0] = 0.05;
    			}
		};
	};

    	class InOut{

       		public:
        	class Initial : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		values[0] = 1;
            		}
        	};

       		class DirichletBoundary : public dolfin::SubDomain {
            		bool inside(const dolfin::Array<double>& x, bool on_boundary) const
            		{
                		// return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS or x[1] <
                		// DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS;
                		return false;
            		}

           		public:
            		// DirichletBoundary(DirichletBoundary& copy) :
            		// dolfin::SubDomain(copy.map_tolerance){}
        	};

        	class Source : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		values[0] = 0;
                		if (x[0] < DOLFIN_EPS)
                    			values[0] = 3;
                		else if (x[0] > 1.0 - DOLFIN_EPS)
                    			values[0] = -2;
            		}
        	};

        	class Neumann : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		values[0] = 0;
            		}
        	};

        	class Velocity : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		for (int i = 0; i < values.size(); i++) {
                    			values[i] = 0;
               		 	}
            		}

           		public:
            		Velocity(std::size_t dim) : dolfin::Expression(dim) {}
        	};

        	class Diffusivity : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		values[0] = 0.5;
            		}
        	};
    	};

    	class ConstSides{

		public:
		class Initial : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		values[0] = 1.0;
            		}
        	};

        	class DirichletBoundary : public dolfin::SubDomain {
            		bool inside(const dolfin::Array<double>& x, bool on_boundary) const
            		{
                		return x[0] < 0.05 or x[0] > 1.0 - 0.05;
            		}
        	};

        	class Source : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		values[0] = 0;
                		if (x[0] > 0.45 && x[0] < 0.55 && x[1] > 0.45 && x[1] < 0.55)
                    			values[0] = 1;
            		}
        	};

        	class Neumann : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		values[0] = 0;
            		}
        	};

        	class Velocity : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		for (int i = 0; i < values.size(); i++) {
                    			values[0] = 0;
                		}
            		}

          		public:
            		Velocity(std::size_t dim) : dolfin::Expression(dim) {}
        	};

        	class Diffusivity : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		values[0] = 0.5;
            		}
        	};
	};

    	class RandomSource {

		public:
        	class Initial : public dolfin::Expression {
			void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		values[0] = -0.5 + x[0];
            		}
        	};

        	class DirichletBoundary : public dolfin::SubDomain {
            		bool inside(const dolfin::Array<double>& x, bool on_boundary) const
            		{
                		return false;
            		}
       	 	};

        	class Source : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		values[0] = 0;
                		if (x[0] > 1.0 - (0.5 + DOLFIN_EPS) &&
                    		x[0] < 1.0 - (0.5 - DOLFIN_EPS) &&
                    		x[1] > 1.0 - (0.5 + DOLFIN_EPS) &&
                    		x[1] < 1.0 - (0.5 - DOLFIN_EPS))
                    			values[0] = 10;
            		}
        	};

        	class Neumann : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		values[0] = 0;
            		}
        	};

        	class Velocity : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		for (int i = 0; i < values.size(); i++) {
                    			values[0] = 0;
                		}
            		}

           		public:
            		Velocity(std::size_t dim) : dolfin::Expression(dim) {}
        	};

        	class Diffusivity : public dolfin::Expression {
			void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
          		{
               			values[0] = 1;
          		}
        	};
	};

    	class DifferentDiffusivities {
	    
		public:
        	class Initial : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		values[0] = 1;
                		if (x[0] > 0.2 && x[0] < 0.4 && x[1] > 0.2 && x[1] < 0.8)
                    			values[0] = 10;
                		else if (x[0] > 0.6 && x[0] < 0.8 && x[1] > 0.2 && x[1] < 0.8)
                    			values[0] = 10;
            		}
        	};

        	class DirichletBoundary : public dolfin::SubDomain {
            		bool inside(const dolfin::Array<double>& x, bool on_boundary) const
            		{
                		return false;
            		}
        	};

        	class Source : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		values[0] = 0;
                		// if(x[0]>1.0-(0.5+DOLFIN_EPS) && x[0]<1.0-(0.5-DOLFIN_EPS) &&
                		// x[1]>1.0-(0.5+DOLFIN_EPS) && x[1]<1.0-(0.5-DOLFIN_EPS))
                		// values[0] = 10;
            		}
        	};

        	class Neumann : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		values[0] = 0;
            		}
        	};

        	class Velocity : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		for (int i = 0; i < values.size(); i++) {
                    			values[0] = 0;
                		}
            		}

           		public:
            		Velocity(std::size_t dim) : dolfin::Expression(dim) {}
        	};

        	class Diffusivity : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		values[0] = 1;
                		if (x[0] > 0.2 && x[0] < 0.4 && x[1] > 0.2 && x[1] < 0.8)
                    			values[0] = 0.5;
                		else if (x[0] > 0.6 && x[0] < 0.8 && x[1] > 0.2 && x[1] < 0.8)
                    			values[0] = 0.75;
            		}
        	};
    	};

    	class TransferFunction {
	
		public:
        	class Initial : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		values[0] = 1;
                		if (x[0] > 0.4 && x[0] < 0.6 && x[1] > 0.4 && x[1] < 0.6)
                    			values[0] = 10;
            		}
        	};

        	class DirichletBoundary : public dolfin::SubDomain {
            		bool inside(const dolfin::Array<double>& x, bool on_boundary) const
            		{
                		return false;
            		}
        	};

        	class Source : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		values[0] = 0;
                		if (((x[0] > 0.4 && x[0] < 0.4 + DOLFIN_EPS) or
				(x[0] > 0.6 - DOLFIN_EPS && x[0] < 0.6)) and
                    		((x[1] > 0.4 && x[1] < 0.4 + DOLFIN_EPS) or
                     		(x[1] > 0.6 - DOLFIN_EPS && x[1] < 0.6))) {
                    			values[0] = -5;
                		} else if (((x[0] > 0.4 - DOLFIN_EPS && x[0] < 0.4) or
                           	(x[0] > 0.6 && x[0] < 0.6 + DOLFIN_EPS)) and
                           	((x[1] > 0.4 - DOLFIN_EPS && x[1] < 0.4) or
                            	(x[1] > 0.6 && x[1] < 0.6 + DOLFIN_EPS))) {
                    			values[0] = 5;
                		}
            		}
        	};

        	class Neumann : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		values[0] = 0;
            		}
        	};

        	class Velocity : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		for (int i = 0; i < values.size(); i++) {
                    		values[i] = 0;
                		}
            		}

           		public:
            		Velocity(std::size_t dim) : dolfin::Expression(dim) {}
        	};

        	class Diffusivity : public dolfin::Expression {
            		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            		{
                		values[0] = 1;
                		if (x[0] > 0.4 && x[0] < 0.6 && x[1] > 0.4 && x[1] < 0.6)
                    			values[0] = 0.5;
            		}
        	};
    	};
}

#endif  // DIFFUSION_FENICS_PDETESTEXAMPLES_H

//
// Created by js on 03.07.16.
//

#ifndef DIFFUSION_FENICS_PDETESTEXAMPLES_H
#define DIFFUSION_FENICS_PDETESTEXAMPLES_H

#include <dolfin.h>

class TestInitial : public dolfin::Expression {
    void eval(dolfin::Array<double> &values, const dolfin::Array<double> &x) const{
        values[0] = 1;
    }
};

//sourceTerm (right-hand side)
class TestSource : public dolfin::Expression {
    void eval(dolfin::Array<double> &values, const dolfin::Array<double> &x) const {
        values[0] = -1;
        if(x[0] > 1 - (0.5+DOLFIN_EPS) and x[0] < 1 - (0.5-DOLFIN_EPS) and x[1] > 1 - (0.5+DOLFIN_EPS) and x[1] < 1 - (0.5-DOLFIN_EPS))
            values[0] = 200;

    }
};

//Normal derivative (used for Neumann boundary condition)
class TestdUdN : public dolfin::Expression {
    void eval(dolfin::Array<double> &values, const dolfin::Array<double> &x) const {
        values[0] = x[0];
    }
};

//Sub domain for Dirichlet boundary condition
    class TestDirichletBoundary : public dolfin::SubDomain {
        bool inside(const dolfin::Array<double> &x, bool on_boundary) const {
            return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS or x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS;
        }
    };

class TestVelocity : public dolfin::Expression {

    void eval(dolfin::Array<double> &values, const dolfin::Array<double> &x) const {
        for(int i = 0; i< values.size(); i++){
            values[i] = 0;
        }
        
    }
    public:
    TestVelocity(std::size_t dim):dolfin::Expression(dim){};
};

class TestDiffusionCoefficient : public dolfin::Expression {
	
	void eval(dolfin::Array<double> &values, const dolfin::Array<double> &x) const {
		values [0] = 0.05;
	} 
};

namespace ConvectionDiffusion {
		
	class TestCases{
	public:
		std::shared_ptr<dolfin::Mesh> getMesh() {return this->mesh;}
		std::shared_ptr<dolfin::Expression> getInitial () {return this->initial;}
		std::shared_ptr<dolfin::SubDomain> getDirichletBoundary() {return this->dirichletBoundary;}
		std::shared_ptr<dolfin::Expression> getSource() {return this->source;}
		std::shared_ptr<dolfin::Expression> getNeumann() {return this->neumann;}
		std::shared_ptr<dolfin::Expression> getVelocity(){return this->velocity;}
		std::shared_ptr<dolfin::Expression> getDiffusivity() {return this->diffusivity;}

		void setMesh(dolfin::Mesh mesh)
			{this->mesh = std::make_shared<dolfin::Mesh>(mesh);}
		void setInitial(dolfin::Expression initial)
			{this->initial = std::make_shared<dolfin::Expression>(initial);}
		void setDirichletBoundary(dolfin::SubDomain dirichletBoundary)
			{this->dirichletBoundary = std::make_shared<dolfin::SubDomain>(dirichletBoundary);}
		void setSource(dolfin::Expression source)
			{this->source = std::make_shared<dolfin::Expression>(source);}
		void setNeumann(dolfin::Expression neumann)
			{this->neumann = std::make_shared<dolfin::Expression>(neumann);}
		void setVelocity(dolfin::Expression velocity)
			{this->velocity = std::make_shared<dolfin::Expression>(velocity);}
		void setDiffusivity(dolfin::Expression diffusivity)
			{this->diffusivity = std::make_shared<dolfin::Expression>(diffusivity);}

	protected:
		std::shared_ptr<dolfin::Mesh> mesh;
		std::shared_ptr<dolfin::Expression> initial;
		std::shared_ptr<dolfin::SubDomain> dirichletBoundary;
		std::shared_ptr<dolfin::Expression> source;
		std::shared_ptr<dolfin::Expression> neumann;
		std::shared_ptr<dolfin::Expression> velocity;
		std::shared_ptr<dolfin::Expression> diffusivity;
	};
	
	class CaseInOut : public TestCases {
		class Initial : public dolfin::Expression {
		
			void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const {
				values [0] = 1;
			}
		};		

		class DirichletBoundary : public dolfin::SubDomain{
			
			bool inside(const dolfin::Array<double>& x, bool on_boundary) const {
				//return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS or x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS;
				return false;
			}

			public:
			//DirichletBoundary(DirichletBoundary& copy) : dolfin::SubDomain(copy.map_tolerance){} 
		};

		class Source : public dolfin::Expression {
			void eval(dolfin::Array<double>&values, const dolfin::Array<double>& x) const {
				values[0] = 0;
				if(x[0] < DOLFIN_EPS)
					values[0] = 3;
				else if(x[0] > 1.0 - DOLFIN_EPS)
					values[0] = -2;
			}	
		};

		class Neumann : public dolfin::Expression {
			void eval(dolfin::Array<double>&values, const dolfin::Array<double>& x) const {
				values[0] = 0;
			}
		};		

		class Velocity : public dolfin::Expression {

			void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const {
				for(int i = 0; i<values.size(); i++){
					values[i] = 0;
				}
			}
			public:
			Velocity(std::size_t dim):dolfin::Expression(dim){}
		};

		class Diffusivity : public dolfin::Expression {
			void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const {
				values[0] = 0.5;
			}
		};

		public:
		CaseInOut(std::size_t dim){

			if(dim == 1) 
				setMesh(dolfin::UnitIntervalMesh(50));
			else if (dim == 2)
				setMesh(dolfin::UnitSquareMesh(50,50));
			else if (dim == 3)
				setMesh(dolfin::UnitCubeMesh(50,50,50));
			else
				setMesh(dolfin::UnitSquareMesh(50,50));
			setInitial(Initial());
			setDirichletBoundary(DirichletBoundary());
			setSource(Source());
			setNeumann(Neumann());
			setVelocity(Velocity(dim));
			setDiffusivity(Diffusivity());
										
		}

	};
} 
#endif //DIFFUSION_FENICS_PDETESTEXAMPLES_H

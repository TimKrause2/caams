#ifndef MBSYSTEM_H
#define MBSYSTEM_H

#include "caams.hpp"
#include "body.h"
#include "constraint.h"
#include "force_element.h"
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;
typedef std::vector< double > state_type;

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::VectorXd EVector;
typedef Eigen::Triplet<double> T;

class SolverParameters
{
public:
	SpMat A;
	EVector x;
	EVector y_lhs;
	EVector y_rhs;
    SolverParameters(long N_eqn):
		A(N_eqn,N_eqn),
		x(N_eqn),
		y_lhs(N_eqn),
		y_rhs(N_eqn)
	{
		y_lhs = Eigen::VectorXd::Zero(N_eqn);
		y_rhs = Eigen::VectorXd::Zero(N_eqn);

	};
    ~SolverParameters(void){}
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class System
{
public:
    std::vector<Body*> bodies;
    std::vector<Constraint*> constraints;
    std::vector<DependentConstraint*> dep_constraints;
    std::vector<ForceElement*> forces;
    std::vector<ForceModifierElement*> force_modifiers;

    long N_eqn;
    long p_index;
    long B_index;
    long N_nzs;
    SolverParameters *sp;
    System(void);
    ~System(void);

    void AddBody(Body *body);
    void AddConstraint(Constraint *constraint);
    void AddDependentConstraint(DependentConstraint *dep_constraint);
    void AddForce(ForceElement *force);
    void AddForceModifier(ForceModifierElement *force_modifier);
    void InitializeSolver(void);
    void rkStateFromVector(const state_type &x);
    void rkAccelerationVector(const state_type &x, state_type &dxdt);
    void StateToVector(state_type &x);
    void VectorToState(const state_type &x);
    void StateToRkState(void);
    void ConstraintReactions(void); // calculates the multipliers for the constraint reactions
    void rkSolve(void);
    void Integrate(double dt);

    void Draw(void);
};

class DxDt
{
    System &mbsystem;
public:
    void operator()(const state_type &x, state_type &dxdt, const double t);
    DxDt(System &mbsystem);
};

#endif

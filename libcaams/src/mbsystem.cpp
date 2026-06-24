#include "mbsystem.h"
#include <iostream>

#define STABILIZATION_ALPHA 8.0
#define STABILIZATION_BETA  8.0

System::System(void)
{
    sp = NULL;
}

System::~System(void)
{
    if(sp) delete sp;
}

void System::AddBody(Body *body)
{
    bodies.push_back(body);
}

void System::AddConstraint(Constraint *constraint)
{
    constraints.push_back(constraint);
}

void System::AddDependentConstraint(DependentConstraint *dep_constraint)
{
    dep_constraints.push_back(dep_constraint);
}

void System::AddForce(ForceElement *force)
{
    forces.push_back(force);
}

void System::AddForceModifier(ForceModifierElement *force_modifier)
{
    force_modifiers.push_back(force_modifier);
}

void System::InitializeSolver(void)
{
	long eqn_index = 0;
    for(auto body:bodies){
        body->eqn_index = eqn_index;
        eqn_index += 7;
    }
    p_index = eqn_index;
    for(auto body:bodies){
        body->p_index = eqn_index;
        eqn_index++;
    }
    B_index = eqn_index;
    for(auto constraint:constraints){
        constraint->eqn_index = eqn_index;
        eqn_index += constraint->N_eqn;
    }
    for(auto dconstraint:dep_constraints){
        dconstraint->eqn_index = eqn_index;
        eqn_index += dconstraint->N_eqn;
        dconstraint->theta_eqn_index = eqn_index;
        eqn_index += 1;
    }

    N_eqn = eqn_index;

    std::cout << "Number of eqations:" << N_eqn << std::endl;

    sp = new SolverParameters(N_eqn);

    // estimate the number of non-zero entries
    int N_nzs = 0;
    N_nzs += bodies.size()*(3+16+2*4);

    for(auto constraint:constraints){
        if(constraint->body1->eqn_index!=0){
            N_nzs += 2*7*constraint->N_eqn;
        }
        if(constraint->body2->eqn_index!=0){
            N_nzs += 2*7*constraint->N_eqn;
        }
    }

    for(auto dconstraint:dep_constraints){
        if(dconstraint->body1->eqn_index!=0){
            N_nzs += 2*7*dconstraint->N_eqn;
        }
        if(dconstraint->body2->eqn_index!=0){
            N_nzs += 2*7*dconstraint->N_eqn;
        }
        N_nzs += dconstraint->N_eqn;
    }

    std::cout << "Number of non-zeros:" << N_nzs << std::endl;
}

void sub(const Eigen::MatrixXd &M, long row, long col, std::vector<T> &nzs){
	for(long i_col=0;i_col<M.cols();i_col++){
		for(long i_row=0;i_row<M.rows();i_row++){
			if(M(i_row,i_col) != 0.0){
				nzs.push_back(T(row+i_row,col+i_col,M(i_row,i_col)));
            }
        }
	}
}

void System::rkSolve(void)
{
    std::vector<T> nzs;
    nzs.reserve(N_nzs);

    // build the non-zero triplet list
    for(auto body:bodies){
        sub(body->N(),body->eqn_index,body->eqn_index,nzs);
        sub(body->J_star(),body->eqn_index+3,body->eqn_index+3,nzs);
		sub(body->rk_p.transpose(),body->p_index,body->eqn_index+3,nzs);
        sub(body->rk_p,body->eqn_index+3,body->p_index,nzs);
		sp->y_lhs.segment<7>(body->eqn_index) = body->b_star();
		sp->y_lhs.segment<1>(body->p_index) = body->rk_p_dot.transpose()*body->rk_p_dot;
		sp->y_rhs.segment<7>(body->eqn_index) = Eigen::Matrix<double,7,1>::Zero();
    }

    for(auto force:forces){
        force->Apply(sp->y_rhs);
    }

    for(auto force:force_modifiers){
        force->Apply(sp->y_rhs);
    }

    for(auto constraint:constraints){
        if(constraint->body1->eqn_index>=0){
            sub(constraint->Body1ModifiedJacobian(),
                constraint->eqn_index,
                constraint->body1->eqn_index,
                nzs);
            sub(constraint->Body1ModifiedJacobian().transpose(),
                constraint->body1->eqn_index,
                constraint->eqn_index,
                nzs);
        }
        if(constraint->body2->eqn_index>=0){
            sub(constraint->Body2ModifiedJacobian(),
                constraint->eqn_index,
                constraint->body2->eqn_index,
                nzs);
            sub(constraint->Body2ModifiedJacobian().transpose(),
                constraint->body2->eqn_index,
                constraint->eqn_index,
                nzs);
        }
        Eigen::VectorXd gamma_t(constraint->ModifiedGamma());
        //std::cout << "gamma:" << std::endl << gamma_t << std::endl;
        gamma_t -= 2.0*STABILIZATION_ALPHA*constraint->dPHI();
        Eigen::VectorXd PHI(constraint->PHI());
        //PHI.print("PHI:");
        gamma_t -= STABILIZATION_BETA*STABILIZATION_BETA*PHI;
        sp->y_rhs.segment(constraint->eqn_index,constraint->N_eqn) = gamma_t;
    }

    for(auto constraint:dep_constraints){
        if(constraint->body1->eqn_index>=0){
            sub(constraint->Body1ModifiedJacobian(),
                constraint->eqn_index,
                constraint->body1->eqn_index,
                nzs);
            sub(constraint->Body1AccelerationEqn(),
                constraint->theta_eqn_index,
                constraint->body1->eqn_index,
                nzs);
            sub(constraint->Body1Jacobian().transpose(),
                constraint->body1->eqn_index,
                constraint->eqn_index,
                nzs);
        }
        if(constraint->body2->eqn_index>=0){
            sub(constraint->Body2ModifiedJacobian(),
                constraint->eqn_index,
                constraint->body2->eqn_index,
                nzs);
            sub(constraint->Body2AccelerationEqn(),
                constraint->theta_eqn_index,
                constraint->body2->eqn_index,
                nzs);
            sub(constraint->Body2Jacobian().transpose(),
                constraint->body2->eqn_index,
                constraint->eqn_index,
                nzs);
        }
        sub(constraint->ThetaModifiedJacobian(),
            constraint->eqn_index,
            constraint->theta_eqn_index,
            nzs);
        Eigen::VectorXd gamma_t(constraint->ModifiedGamma());
        //std::cout << "gamma:" << std::endl << gamma_t << std::endl;
        gamma_t -= 2.0*STABILIZATION_ALPHA*constraint->dPHI();
        Eigen::VectorXd PHI(constraint->PHI());
        //PHI.print("PHI:");
        gamma_t -= STABILIZATION_BETA*STABILIZATION_BETA*PHI;
        sp->y_rhs.segment(constraint->eqn_index,constraint->N_eqn) = gamma_t;
    }

    sp->A.setFromTriplets(nzs.begin(), nzs.end());

    Eigen::SparseLU<SpMat> solver;

	solver.analyzePattern(sp->A);
	solver.factorize(sp->A);

	sp->x = solver.solve(sp->y_rhs - sp->y_lhs);

}

void System::Integrate(double dt)
{
    if(!sp)return;

    // create the state vector
    state_type x(bodies.size()*14);

    // Initailize the state vector
    StateToVector(x);

    // Integrate using boost odeint
    DxDt dxdt(*this);
    runge_kutta4< state_type > rk4; // Stepper
    integrate_const(rk4, dxdt, x, 0.0, dt, dt);

    // Transfer the vector back to the state
    VectorToState(x);
}

void System::Draw(void)
{
    for(auto body:bodies){
        body->Draw();
    }
    for(auto constraint:constraints){
        constraint->Draw();
    }
    for(auto force:forces){
        force->Draw();
    }
}

void System::rkStateFromVector(const state_type &x_raw)
{
    // Map x into an Eigen vector
    Eigen::Map<const Eigen::VectorXd> x(x_raw.data(), x_raw.size());

    int b_index = 0;
    for(auto &body:bodies){
        int offset = b_index*14;
        body->rk_r = x.segment<3>(offset);
        body->rk_p = x.segment<4>(offset+3);
        body->rk_p.normalize();
        body->rk_r_dot = x.segment<3>(offset+7);
        body->rk_p_dot = x.segment<4>(offset+10);
        b_index++;
    }
}

void System::rkAccelerationVector(const state_type &x_raw, state_type &dxdt_raw)
{
    // Map the inputs into Eigen vectors
    Eigen::Map<const Eigen::VectorXd> x(x_raw.data(), x_raw.size());
    Eigen::Map<Eigen::VectorXd> dxdt(dxdt_raw.data(), dxdt_raw.size());
    int b_index=0;
    for(auto &body:bodies){
        int offset = b_index*14;
        // extract the accelerations from the solution
        Eigen::Matrix<double,7,1> q_ddot = sp->x.segment<7>(body->eqn_index);
        // extranct the velocity from x
        Eigen::Matrix<double,7,1> q_dot = x.segment<7>(offset+7);
        // update dxdt
        dxdt.segment<7>(offset) = q_dot;
        dxdt.segment<7>(offset+7) = q_ddot;
        b_index++;
    }
}

void System::StateToVector(state_type &x_raw)
{
    Eigen::Map<Eigen::VectorXd> x(x_raw.data(), x_raw.size());
    int b_index=0;
    for(auto &body:bodies){
        int offset = b_index*14;
        x.segment<3>(offset) = body->r;
        x.segment<4>(offset+3) = body->p;
        x.segment<3>(offset+7) = body->r_dot;
        x.segment<4>(offset+10) = body->p_dot;
        b_index++;
    }
}

void System::VectorToState(const state_type &x_raw)
{
    Eigen::Map<const Eigen::VectorXd> x(x_raw.data(), x_raw.size());
    int b_index = 0;
    for(auto &body:bodies){
        int offset = b_index*14;
        body->r = x.segment<3>(offset);
        body->p = x.segment<4>(offset+3);
        body->p.normalize();
        body->r_dot = x.segment<3>(offset+7);
        body->p_dot = x.segment<4>(offset+10);
        double sigma = body->p.dot(body->p_dot);
        body->p_dot -= sigma*body->p;
        b_index++;
    }
}

void System::StateToRkState(void)
{
    for(auto &body:bodies){
        body->rk_r = body->r;
        body->rk_p = body->p;
        body->rk_r_dot = body->r_dot;
        body->rk_p_dot = body->p_dot;
    }
}

void System::ConstraintReactions(void)
{
    // transfer the current state to the rk state
    StateToRkState();

    // solve for the multipliers
    rkSolve();

    // transfer the results to the constraints
    for(auto &constraint:constraints){
        constraint->lambda = sp->x.segment(
                    constraint->eqn_index,
                    constraint->N_eqn);
    }
}

DxDt::DxDt(System &mbsystem):
    mbsystem(mbsystem)
{
}

void DxDt::operator()(const state_type &x, state_type &dxdt, const double t)
{
    // Transfer the state from x to the rk variables
    mbsystem.rkStateFromVector(x);

    // Solve for the accelerations
    mbsystem.rkSolve();

    // Transfer the derivatives to dxdt
    mbsystem.rkAccelerationVector(x, dxdt);
}

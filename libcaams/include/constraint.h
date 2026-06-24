#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include "caams.hpp"
#include "body.h"

class Constraint
{
public:
    Body *body1;
    Body *body2;
    int N_eqn;
    int eqn_index;
    Eigen::VectorXd lambda;
    Constraint(
            Body *body1,
            Body *body2,
            int N_eqn);
	~Constraint(){}
    void Reaction(
            Eigen::Vector3d &body1_force,
            Eigen::Vector3d &body1_torque,
            Eigen::Vector3d &body2_force,
            Eigen::Vector3d &body2_torque);
    Eigen::VectorXd dPHI(void);
	Eigen::Vector3d h(const Eigen::Vector4d &p_dot,
					  const Eigen::Vector3d &s_p);
	virtual Eigen::MatrixXd Body1Jacobian(void)=0;
	virtual Eigen::MatrixXd Body2Jacobian(void)=0;
	virtual Eigen::MatrixXd Body1ModifiedJacobian(void)=0;
	virtual Eigen::MatrixXd Body2ModifiedJacobian(void)=0;
	virtual Eigen::VectorXd PHI(void)=0;
	virtual Eigen::VectorXd ModifiedGamma(void)=0;
    virtual void Draw(void)=0;
    virtual void DrawReaction(void)=0;
};

class DependentConstraint : public Constraint
{
public:
    int theta_eqn_index;
    double theta;
    double theta0;
    double theta_dot;
    double rk_theta;
    double rk_theta_dot;
    DependentConstraint(
            Body *body1, Body *body2, int N_eqn,
            double theta0, double theta_dot0);
    Eigen::Matrix<double,1,4> k_theta_dot;
    Eigen::Matrix<double,1,4> k_theta_ddot;
    virtual Eigen::VectorXd ThetaModifiedJacobian(void)=0;
    virtual Eigen::MatrixXd Body1AccelerationEqn(void)=0;
    virtual Eigen::MatrixXd Body2AccelerationEqn(void)=0;
};

class ScrewJoint : public DependentConstraint
{
public:
    Eigen::Vector3d s1_p; // point on axis in body1 frame
    Eigen::Vector3d s2_p; // point on axis in body2 frame
    // normal vectors - must be unit length!
    Eigen::Vector3d na_p; // axis in body one frame
    Eigen::Vector3d n1_p; // normal vector in body1 frame
    Eigen::Vector3d n2A_p; // normal vector in body2 frame. initial angle is theta0
    Eigen::Vector3d n2B_p; // bi-normal vector in body2 frame
    double alpha; // pitch rate (distance/radian)
    double d0; // initial distance. determined from the initial conditions
    ScrewJoint(
            Body *body1,
            Body *body2,
            Eigen::Vector3d s1_p,
            Eigen::Vector3d s2_p,
            Eigen::Vector3d na_p,
            Eigen::Vector3d n1_p,
            Eigen::Vector3d n2A_p,
            Eigen::Vector3d n2B_p,
            double alpha,
            double theta0,
            double theta_dot0);
    ~ScrewJoint(void) {}
private:
    bool theta_in_A_zone(void);
    Eigen::Vector3d d_vec(void);
    Eigen::Vector3d d_dot_vec(void);
public:
    virtual Eigen::MatrixXd Body1Jacobian(void);
    virtual Eigen::MatrixXd Body2Jacobian(void);
    virtual Eigen::MatrixXd Body1ModifiedJacobian(void);
    virtual Eigen::MatrixXd Body2ModifiedJacobian(void);
    virtual Eigen::VectorXd PHI(void);
    virtual Eigen::VectorXd ModifiedGamma(void);
    virtual void Draw(void);
    virtual void DrawReaction(void);
    virtual Eigen::VectorXd ThetaModifiedJacobian(void);
    virtual Eigen::MatrixXd Body1AccelerationEqn(void);
    virtual Eigen::MatrixXd Body2AccelerationEqn(void);
};

class SphericalJoint : public Constraint
{
public:
	Eigen::Vector3d s1_p;
	Eigen::Vector3d s2_p;
    SphericalJoint(
            Body *body1,
            Body *body2,
			Eigen::Vector3d s1_p,
			Eigen::Vector3d s2_p);
	~SphericalJoint(){}
	virtual Eigen::MatrixXd Body1Jacobian(void);
	virtual Eigen::MatrixXd Body2Jacobian(void);
	virtual Eigen::MatrixXd Body1ModifiedJacobian(void);
	virtual Eigen::MatrixXd Body2ModifiedJacobian(void);
	virtual Eigen::VectorXd PHI(void);
	virtual Eigen::VectorXd ModifiedGamma(void);
	virtual void Draw(void);
    virtual void DrawReaction(void);

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};



class SphericalSphericalJoint : public Constraint
{
public:
	Eigen::Vector3d s1_p;
	Eigen::Vector3d s2_p;
    double length;
    SphericalSphericalJoint(
            Body *body1,
            Body *body2,
			Eigen::Vector3d s1_p,
			Eigen::Vector3d s2_p,
            double length);
	~SphericalSphericalJoint(){}
	virtual Eigen::MatrixXd Body1Jacobian(void);
	virtual Eigen::MatrixXd Body2Jacobian(void);
	virtual Eigen::MatrixXd Body1ModifiedJacobian(void);
	virtual Eigen::MatrixXd Body2ModifiedJacobian(void);
	virtual Eigen::VectorXd PHI(void);
	virtual Eigen::VectorXd ModifiedGamma(void);
	virtual void Draw(void);
    virtual void DrawReaction(void);

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class Normal1_1 : public Constraint
{
public:
	Eigen::Vector3d s1_p;
	Eigen::Vector3d s2_p;
    Normal1_1(
            Body *body1,
            Body *body2,
			Eigen::Vector3d s1_p,
			Eigen::Vector3d s2_p);
    ~Normal1_1(){};
	virtual Eigen::MatrixXd Body1Jacobian(void);
	virtual Eigen::MatrixXd Body2Jacobian(void);
	virtual Eigen::MatrixXd Body1ModifiedJacobian(void);
	virtual Eigen::MatrixXd Body2ModifiedJacobian(void);
	virtual Eigen::VectorXd PHI(void);
	virtual Eigen::VectorXd ModifiedGamma(void);
    virtual void Draw(void);
    virtual void DrawReaction(void);

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class Normal2_1 : public Constraint
{
public:
	Eigen::Vector3d s1_p;
	Eigen::Vector3d s1B_p;
	Eigen::Vector3d s2B_p;
    Normal2_1(
            Body *body1,
            Body *body2,
			Eigen::Vector3d s1_p,
			Eigen::Vector3d s1B_p,
			Eigen::Vector3d s2B_p);
    ~Normal2_1(void){}
	virtual Eigen::MatrixXd Body1Jacobian(void);
	virtual Eigen::MatrixXd Body2Jacobian(void);
	virtual Eigen::MatrixXd Body1ModifiedJacobian(void);
	virtual Eigen::MatrixXd Body2ModifiedJacobian(void);
	virtual Eigen::VectorXd PHI(void);
	virtual Eigen::VectorXd ModifiedGamma(void);
    virtual void Draw(void);
    virtual void DrawReaction(void);

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class Parallel1_2 : public Constraint
{
public:
	Eigen::Vector3d s1_p;
	Eigen::Vector3d s2_p;
    Parallel1_2(
            Body *body1,
            Body *body2,
			Eigen::Vector3d s1_p,
			Eigen::Vector3d s2_p);
    ~Parallel1_2(void){}
	virtual Eigen::MatrixXd Body1Jacobian(void);
	virtual Eigen::MatrixXd Body2Jacobian(void);
	virtual Eigen::MatrixXd Body1ModifiedJacobian(void);
	virtual Eigen::MatrixXd Body2ModifiedJacobian(void);
	virtual Eigen::VectorXd PHI(void);
	virtual Eigen::VectorXd ModifiedGamma(void);
	virtual void Draw(void);
    virtual void DrawReaction(void);

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class Parallel2_2 : public Constraint
{
public:
	Eigen::Vector3d s1_p;
	Eigen::Vector3d s1B_p;
	Eigen::Vector3d s2B_p;
    Parallel2_2(
            Body *body1,
            Body *body2,
			Eigen::Vector3d s1_p,
			Eigen::Vector3d s1B_p,
			Eigen::Vector3d s2B_p);
    ~Parallel2_2(void){}
	virtual Eigen::MatrixXd Body1Jacobian(void);
	virtual Eigen::MatrixXd Body2Jacobian(void);
	virtual Eigen::MatrixXd Body1ModifiedJacobian(void);
	virtual Eigen::MatrixXd Body2ModifiedJacobian(void);
	virtual Eigen::VectorXd PHI(void);
	virtual Eigen::VectorXd ModifiedGamma(void);
	virtual void Draw(void);
    virtual void DrawReaction(void);

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class PlanarJoint : public Constraint
{
public:
    Eigen::Vector3d s1_p; // point on plane in body1 space
    Eigen::Vector3d s2_p; // point on plane in body2 space
    Eigen::Vector3d u2_p; // plane normal in body2 space
    PlanarJoint(
            Body *body1,
            Body *body2,
            Eigen::Vector3d s1_p,
            Eigen::Vector3d s2_p,
            Eigen::Vector3d u2_p);
    ~PlanarJoint(void){}
    virtual Eigen::MatrixXd Body1Jacobian(void);
    virtual Eigen::MatrixXd Body2Jacobian(void);
    virtual Eigen::MatrixXd Body1ModifiedJacobian(void);
    virtual Eigen::MatrixXd Body2ModifiedJacobian(void);
    virtual Eigen::VectorXd PHI(void);
    virtual Eigen::VectorXd ModifiedGamma(void);
    virtual void Draw(void);
    virtual void DrawReaction(void);

public:
};

class ScrewJoint_1 : public Constraint
{
public:
    Eigen::Vector3d s1_p; // point on axis in body1 space
    Eigen::Vector3d s2_p; // point on axis in body2 space
    Eigen::Vector3d u1_p; // normal in body1 space
    Eigen::Vector3d u2_p; // perpendicular normal in body2 space
    Eigen::Vector3d a1_p; // axis of rotation in body one space
    double beta;
    double theta0;
    ScrewJoint_1(
            Body *body1,
            Body *body2,
            Eigen::Vector3d s1_p,
            Eigen::Vector3d s2_p,
            Eigen::Vector3d u1_p,
            Eigen::Vector3d u2_p,
            Eigen::Vector3d a1_p,
            double beta);
    ~ScrewJoint_1(void){}
    virtual Eigen::MatrixXd Body1Jacobian(void);
    virtual Eigen::MatrixXd Body2Jacobian(void);
    virtual Eigen::MatrixXd Body1ModifiedJacobian(void);
    virtual Eigen::MatrixXd Body2ModifiedJacobian(void);
    virtual Eigen::VectorXd PHI(void);
    virtual Eigen::VectorXd ModifiedGamma(void);
    virtual void Draw(void);
    virtual void DrawReaction(void);

public:
};

//
// Composite constraints
//
class UniversalJoint : public Constraint
{
public:
    SphericalJoint sphericalJoint;
    Normal1_1 normalJoint;
    Eigen::Vector3d s1_p; // center of joint in body1 space
    Eigen::Vector3d s2_p; // center of joint in body2 space
    Eigen::Vector3d h1_p; // rotation axis in body1 space
    Eigen::Vector3d h2_p; // rotation axis in body2 space
    UniversalJoint(
            Body *body1,
            Body *body2,
            Eigen::Vector3d s1_p, // center of joint in body1 space
            Eigen::Vector3d s2_p, // center of joint in body2 space
            Eigen::Vector3d h1_p, // rotation axis in body1 space
            Eigen::Vector3d h2_p); // rotation axis in body2 space
    ~UniversalJoint(void){}
    virtual Eigen::MatrixXd Body1Jacobian(void);
    virtual Eigen::MatrixXd Body2Jacobian(void);
    virtual Eigen::MatrixXd Body1ModifiedJacobian(void);
    virtual Eigen::MatrixXd Body2ModifiedJacobian(void);
    virtual Eigen::VectorXd PHI(void);
    virtual Eigen::VectorXd ModifiedGamma(void);
    virtual void Draw(void);
    virtual void DrawReaction(void);

public:
};

class RevoluteJoint : public Constraint
{
public:
    SphericalJoint sphericalJoint;
    Parallel1_2    parallelJoint;
    Eigen::Vector3d s1_p; // center of joint in body1 space
    Eigen::Vector3d s2_p; // center of joint in body2 space
    Eigen::Vector3d h1_p; // rotation axis in body1 space
    Eigen::Vector3d h2_p; // rotation axis in body2 space
    RevoluteJoint(
            Body *body1,
            Body *body2,
            Eigen::Vector3d s1_p, // center of joint in body1 space
            Eigen::Vector3d s2_p, // center of joint in body2 space
            Eigen::Vector3d h1_p, // rotation axis in body1 space
            Eigen::Vector3d h2_p); // rotation axis in body2 space
    ~RevoluteJoint(void){}
    virtual Eigen::MatrixXd Body1Jacobian(void);
    virtual Eigen::MatrixXd Body2Jacobian(void);
    virtual Eigen::MatrixXd Body1ModifiedJacobian(void);
    virtual Eigen::MatrixXd Body2ModifiedJacobian(void);
    virtual Eigen::VectorXd PHI(void);
    virtual Eigen::VectorXd ModifiedGamma(void);
    virtual void Draw(void);
    virtual void DrawReaction(void);

public:
};

class TranslationalJoint : public Constraint
{
public:
    Parallel1_2 parallel1Joint;
    Parallel2_2 parallel2Joint;
    Normal1_1 normalJoint;
    Eigen::Vector3d a1_p; // rotation axis in body1 space
    Eigen::Vector3d a2_p; // rotation axis in body2 space
    Eigen::Vector3d s1B_p; // point on axis in body1 space
    Eigen::Vector3d s2B_p; // point on axis in body2 space
    Eigen::Vector3d h1_p; // normal axis in body1 space
    Eigen::Vector3d h2_p; // normal axis in body2 space
    TranslationalJoint(
            Body *body1,
            Body *body2,
            Eigen::Vector3d a1_p, // rotation axis in body1 space
            Eigen::Vector3d a2_p, // rotation axis in body2 space
            Eigen::Vector3d s1B_p, // point on axis in body1 space
            Eigen::Vector3d s2B_p, // point on axis in body2 space
            Eigen::Vector3d h1_p, // normal axis in body1 space
            Eigen::Vector3d h2_p); // normal axis in body2 space
    ~TranslationalJoint(void){}
    virtual Eigen::MatrixXd Body1Jacobian(void);
    virtual Eigen::MatrixXd Body2Jacobian(void);
    virtual Eigen::MatrixXd Body1ModifiedJacobian(void);
    virtual Eigen::MatrixXd Body2ModifiedJacobian(void);
    virtual Eigen::VectorXd PHI(void);
    virtual Eigen::VectorXd ModifiedGamma(void);
    virtual void Draw(void);
    virtual void DrawReaction(void);

public:
};

class CylindricalJoint : public Constraint
{
public:
    Parallel1_2 parallel1Joint;
    Parallel2_2 parallel2Joint;
    Eigen::Vector3d a1_p; // rotation axis in body1 space
    Eigen::Vector3d a2_p; // rotation axis in body2 space
    Eigen::Vector3d s1B_p; // point on axis in body1 space
    Eigen::Vector3d s2B_p; // point on axis in body2 space
    CylindricalJoint(
            Body *body1,
            Body *body2,
            Eigen::Vector3d a1_p, // rotation axis in body1 space
            Eigen::Vector3d a2_p, // rotation axis in body2 space
            Eigen::Vector3d s1B_p, // point on axis in body1 space
            Eigen::Vector3d s2B_p); // point on axis in body2 space
    ~CylindricalJoint(void){}
    virtual Eigen::MatrixXd Body1Jacobian(void);
    virtual Eigen::MatrixXd Body2Jacobian(void);
    virtual Eigen::MatrixXd Body1ModifiedJacobian(void);
    virtual Eigen::MatrixXd Body2ModifiedJacobian(void);
    virtual Eigen::VectorXd PHI(void);
    virtual Eigen::VectorXd ModifiedGamma(void);
    virtual void Draw(void);
    virtual void DrawReaction(void);

public:
};

#endif

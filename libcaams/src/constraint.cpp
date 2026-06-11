#include "constraint.h"
#include <stdio.h>
#include <iostream>
#include <GL/gl.h>
#include <GL/freeglut.h>
#define GLM_FORCE_RADIANS
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/type_ptr.hpp>

Constraint::Constraint(
        Body *body1,
        Body *body2,
        int N_eqn):
    body1(body1),
    body2(body2),
    N_eqn(N_eqn)
{
    k_lambda.resize(N_eqn, 4);
    lambda.resize(N_eqn);
}

void Constraint::Reaction(
        Eigen::Vector3d &body1_force,
        Eigen::Vector3d &body1_torque,
        Eigen::Vector3d &body2_force,
        Eigen::Vector3d &body2_torque)
{
    Eigen::Matrix<double, 7, 1> g1 = Body1Jacobian().transpose()*lambda;
    Eigen::Matrix<double, 7, 1> g2 = Body2Jacobian().transpose()*lambda;
    body1_force = g1.head<3>();
    body1_torque = 0.5*caams::G(body1->p)*g1.tail<4>();
    body2_force = g2.head<3>();
    body2_torque = 0.5*caams::G(body2->p)*g2.tail<4>();
}

Eigen::VectorXd Constraint::dPHI(void)
{
	Eigen::Matrix<double,7,1> q1_dot;
	Eigen::Matrix<double,7,1> q2_dot;

	q1_dot.head<3>() = body1->rk_r_dot;
	q1_dot.tail<4>() = body1->rk_p_dot;
	q2_dot.head<3>() = body2->rk_r_dot;
	q2_dot.tail<4>() = body2->rk_p_dot;

    return Body1Jacobian()*q1_dot + Body2Jacobian()*q2_dot;
}

Eigen::Vector3d Constraint::h(const Eigen::Vector4d &p_dot,
							  const Eigen::Vector3d &s_p)
{
	return -2.0*caams::G(p_dot)*caams::L(p_dot).transpose()*s_p;
}

SphericalJoint::SphericalJoint(
        Body *body1,
        Body *body2,
		Eigen::Vector3d s1_p,
		Eigen::Vector3d s2_p):
    Constraint(body1,body2,3),
    s1_p(s1_p),
    s2_p(s2_p)
{
}

Eigen::MatrixXd SphericalJoint::Body1Jacobian(void)
{
	Eigen::Matrix<double,3,7> r;
	r.block<3,3>(0,0) = Eigen::Matrix3d::Identity();
	r.block<3,4>(0,3) = 2.0*(caams::G(body1->rk_p)*caams::a_minus(s1_p)+s1_p*body1->rk_p.transpose());
    return r;
}

Eigen::MatrixXd SphericalJoint::Body2Jacobian(void)
{
	Eigen::Matrix<double,3,7> r;
	r.block<3,3>(0,0) = -Eigen::Matrix3d::Identity();
	r.block<3,4>(0,3) =-2.0*(caams::G(body2->rk_p)*caams::a_minus(s2_p)+s2_p*body2->rk_p.transpose());
	return r;
}

Eigen::MatrixXd SphericalJoint::Body1ModifiedJacobian(void)
{
	Eigen::Matrix<double,3,7> r;
	r.block<3,3>(0,0) = Eigen::Matrix3d::Identity();
	r.block<3,4>(0,3) = 2.0*caams::G(body1->rk_p)*caams::a_minus(s1_p);
	return r;
}

Eigen::MatrixXd SphericalJoint::Body2ModifiedJacobian(void)
{
	Eigen::Matrix<double,3,7> r;
	r.block<3,3>(0,0) = -Eigen::Matrix3d::Identity();
	r.block<3,4>(0,3) =-2.0*caams::G(body2->rk_p)*caams::a_minus(s2_p);
	return r;
}

Eigen::VectorXd SphericalJoint::PHI(void)
{
    return body1->rk_r + caams::Ap(body1->rk_p)*s1_p
           - body2->rk_r - caams::Ap(body2->rk_p)*s2_p;
}

Eigen::VectorXd SphericalJoint::ModifiedGamma(void)
{
    return h(body1->rk_p_dot,s1_p) - h(body2->rk_p_dot,s2_p);
}

void SphericalJoint::Draw(void)
{
	Eigen::Vector3d s1 = body1->r + caams::Ap(body1->p)*s1_p;
	Eigen::Vector3d s2 = body2->r + caams::Ap(body2->p)*s2_p;

    glBegin(GL_LINES);
        glColor3f(1.0f,1.0f,0.0f);
		glVertex3dv(body1->r.data());
		glVertex3dv(s1.data());
        glColor3f(0.0f,1.0f,1.0f);
		glVertex3dv(body2->r.data());
		glVertex3dv(s2.data());
    glEnd();
}

void SphericalJoint::DrawReaction(void)
{
    Eigen::Vector3d s1 = body1->r + caams::Ap(body1->p)*s1_p;
    Eigen::Vector3d s2 = body2->r + caams::Ap(body2->p)*s2_p;
    Eigen::Vector3d body1_force;
    Eigen::Vector3d body1_torque;
    Eigen::Vector3d body2_force;
    Eigen::Vector3d body2_torque;

    Reaction(body1_force,body1_torque, body2_force, body2_torque);
    body1_force.normalize();
    body1_torque.normalize();
    body2_force.normalize();
    body2_torque.normalize();
    body1_force += s1;
    body1_torque += s1;
    body2_force += s2;
    body2_torque += s2;
    glBegin(GL_LINES);
        glColor3d(1.0,1.0,0.0);
        glVertex3dv(s1.data());
        glVertex3dv(body1_force.data());
        glColor3d(1.0,0.0,1.0);
        glVertex3dv(s2.data());
        glVertex3dv(body2_force.data());
        glColor3d(0.0,1.0,1.0);
        glVertex3dv(s1.data());
        glVertex3dv(body1_torque.data());
        glColor3d(1.0,0.5,0.1);
        glVertex3dv(s2.data());
        glVertex3dv(body2_torque.data());
    glEnd();
}

SphericalSphericalJoint::SphericalSphericalJoint(
        Body *body1,
        Body *body2,
		Eigen::Vector3d s1_p,
		Eigen::Vector3d s2_p,
        double length):
    Constraint(body1,body2,1),
    s1_p(s1_p),
    s2_p(s2_p),
    length(length)
{
}

Eigen::MatrixXd SphericalSphericalJoint::Body1Jacobian(void)
{
	Eigen::Vector3d d =
                body2->rk_r+caams::Ap(body2->rk_p)*s2_p
				-body1->rk_r-caams::Ap(body1->rk_p)*s1_p;
	Eigen::Matrix<double,1,7> r;
	r.head<3>() = -2.0*d.transpose();
	Eigen::Matrix3x4d B1 = 2.0*(caams::G(body1->rk_p)*caams::a_minus(s1_p)
					 + s1_p*body1->rk_p.transpose());
	r.tail<4>() = -2.0*d.transpose()*B1;
    return r;
}

Eigen::MatrixXd SphericalSphericalJoint::Body2Jacobian(void)
{
	Eigen::Vector3d d =
				body2->rk_r+caams::Ap(body2->rk_p)*s2_p
				-body1->rk_r-caams::Ap(body1->rk_p)*s1_p;
	Eigen::Matrix<double,1,7> r;
	r.head<3>() = 2.0*d.transpose();
	Eigen::Matrix3x4d B2 = 2.0*(caams::G(body2->rk_p)*caams::a_minus(s2_p)
					 + s2_p*body2->rk_p.transpose());
	r.tail<4>() = 2.0*d.transpose()*B2;
	return r;
}

Eigen::MatrixXd SphericalSphericalJoint::Body1ModifiedJacobian(void)
{
	Eigen::Vector3d d =
				body2->rk_r+caams::Ap(body2->rk_p)*s2_p
				-body1->rk_r-caams::Ap(body1->rk_p)*s1_p;
	Eigen::Matrix<double,1,7> r;
	r.head<3>() = -2.0*d.transpose();
	Eigen::Matrix3x4d B1 = 2.0*caams::G(body1->rk_p)*caams::a_minus(s1_p);
	r.tail<4>() = -2.0*d.transpose()*B1;
	return r;
}

Eigen::MatrixXd SphericalSphericalJoint::Body2ModifiedJacobian(void)
{
	Eigen::Vector3d d =
				body2->rk_r+caams::Ap(body2->rk_p)*s2_p
				-body1->rk_r-caams::Ap(body1->rk_p)*s1_p;
	Eigen::Matrix<double,1,7> r;
	r.head<3>() = 2.0*d.transpose();
	Eigen::Matrix3x4d B2 = 2.0*caams::G(body2->rk_p)*caams::a_minus(s2_p);
	r.tail<4>() = 2.0*d.transpose()*B2;
	return r;
}

Eigen::VectorXd SphericalSphericalJoint::PHI(void)
{
	Eigen::Vector3d d =
				body2->rk_r+caams::Ap(body2->rk_p)*s2_p
				-body1->rk_r-caams::Ap(body1->rk_p)*s1_p;
	return d.transpose()*d - Eigen::Matrix<double,1,1>(length*length);
}

Eigen::VectorXd SphericalSphericalJoint::ModifiedGamma(void)
{
	Eigen::Vector3d d =
				body2->rk_r+caams::Ap(body2->rk_p)*s2_p
				-body1->rk_r-caams::Ap(body1->rk_p)*s1_p;
	Eigen::Vector3d omega1 = 2.0*caams::G(body1->rk_p)*body1->rk_p_dot;
	Eigen::Vector3d omega2 = 2.0*caams::G(body2->rk_p)*body2->rk_p_dot;
	Eigen::Vector3d d_dot(
                body2->rk_r_dot + caams::SS(omega2)*caams::Ap(body2->rk_p)*s2_p
                -body1->rk_r_dot - caams::SS(omega1)*caams::Ap(body1->rk_p)*s1_p);

	return 2.0*(d.transpose()*(h(body2->rk_p_dot,s2_p)-h(body1->rk_p_dot,s1_p))-d_dot.transpose()*d_dot);
}

void SphericalSphericalJoint::Draw(void)
{
	Eigen::Vector3d s1 = body1->r + caams::Ap(body1->p)*s1_p;
	Eigen::Vector3d s2 = body2->r + caams::Ap(body2->p)*s2_p;
	//Eigen::Vector3d d(s2-s1);
	//double l = d.norm();
	//std::cout << "SphericalSphericalJoint::Draw l:" << l << std::endl;

    glBegin(GL_LINES);
        glColor3f(1.0f,0.5f,0.5f);
		glVertex3dv(body1->r.data());
		glVertex3dv(s1.data());
        glColor3f(0.5f,0.5f,0.5f);
		glVertex3dv(s1.data());
		glVertex3dv(s2.data());
        glColor3f(0.5f,1.0f,0.5f);
		glVertex3dv(s2.data());
		glVertex3dv(body2->r.data());
    glEnd();
}

void SphericalSphericalJoint::DrawReaction(void)
{
}

Normal1_1::Normal1_1(
                Body *body1,
                Body *body2,
				Eigen::Vector3d s1_p,
				Eigen::Vector3d s2_p):
    Constraint(body1,body2,1),
    s1_p(s1_p),
    s2_p(s2_p)
{
}

Eigen::MatrixXd Normal1_1::Body1Jacobian(void)
{
	Eigen::Vector3d s2 = caams::Ap(body2->rk_p)*s2_p;
	Eigen::Matrix3x4d C1 = 2.0*(caams::G(body1->rk_p)*caams::a_minus(s1_p)+s1_p*body1->rk_p.transpose());
	Eigen::Matrix<double,1,7> r;
	r.head<3>() = Eigen::RowVector3d::Zero();
	r.tail<4>() = s2.transpose()*C1;
    return r;
}

Eigen::MatrixXd Normal1_1::Body2Jacobian(void)
{
	Eigen::Vector3d s1 = caams::Ap(body1->rk_p)*s1_p;
	Eigen::Matrix3x4d C2 = 2.0*(caams::G(body2->rk_p)*caams::a_minus(s2_p)+s2_p*body2->rk_p.transpose());
	Eigen::Matrix<double,1,7> r;
	r.head<3>() = Eigen::RowVector3d::Zero();
	r.tail<4>() = s1.transpose()*C2;
	return r;
}

Eigen::MatrixXd Normal1_1::Body1ModifiedJacobian(void)
{
	Eigen::Vector3d s2 = caams::Ap(body2->rk_p)*s2_p;
	Eigen::Matrix<double,1,7> r;
	r.head<3>() = Eigen::RowVector3d::Zero();
	r.tail<4>() = s2.transpose()*caams::G(body1->rk_p)*caams::a_minus(s1_p);
    return r;
}

Eigen::MatrixXd Normal1_1::Body2ModifiedJacobian(void)
{
	Eigen::Vector3d s1 = caams::Ap(body1->rk_p)*s1_p;
	Eigen::Matrix<double,1,7> r;
	r.head<3>() = Eigen::RowVector3d::Zero();
	r.tail<4>() = s1.transpose()*caams::G(body2->rk_p)*caams::a_minus(s2_p);
	return r;
}

Eigen::VectorXd Normal1_1::PHI(void)
{
	Eigen::Vector3d s1 = caams::Ap(body1->rk_p)*s1_p;
	Eigen::Vector3d s2 = caams::Ap(body2->rk_p)*s2_p;
	return s1.transpose()*s2;
}

Eigen::VectorXd Normal1_1::ModifiedGamma(void)
{
	Eigen::Vector3d s1 = caams::Ap(body1->rk_p)*s1_p;
	Eigen::Vector3d s2 = caams::Ap(body2->rk_p)*s2_p;
	Eigen::Vector3d omega1 = 2.0*caams::G(body1->rk_p)*body1->rk_p_dot;
	Eigen::Vector3d omega2 = 2.0*caams::G(body2->rk_p)*body2->rk_p_dot;
	Eigen::Vector3d s1_dot = caams::SS(omega1)*s1;
	Eigen::Vector3d s2_dot = caams::SS(omega2)*s2;
	return s1.transpose()*h(body2->rk_p_dot,s2_p) + s2.transpose()*h(body1->rk_p_dot,s1_p)
			- 2.0*s1_dot.transpose()*s2_dot;
}


void Normal1_1::Draw(void)
{

}

void Normal1_1::DrawReaction(void)
{
}

Normal2_1::Normal2_1(
        Body *body1,
        Body *body2,
		Eigen::Vector3d s1_p,
		Eigen::Vector3d s1B_p,
		Eigen::Vector3d s2B_p):
    Constraint(body1,body2,1),
    s1_p(s1_p),
    s1B_p(s1B_p),
    s2B_p(s2B_p)
{
}

Eigen::MatrixXd Normal2_1::Body1Jacobian(void)
{
	Eigen::Vector3d s1 = caams::Ap(body1->rk_p)*s1_p;
	Eigen::Vector3d d =
				body2->rk_r+caams::Ap(body2->rk_p)*s2B_p
				-body1->rk_r-caams::Ap(body1->rk_p)*s1B_p;
	Eigen::Matrix3x4d B1 = 2.0*(caams::G(body1->rk_p)*caams::a_minus(s1B_p)+s1B_p*body1->rk_p.transpose());
	Eigen::Matrix3x4d C1 = 2.0*(caams::G(body1->rk_p)*caams::a_minus(s1_p)+s1_p*body1->rk_p.transpose());
	Eigen::Matrix<double,1,7> r;
	r.head<3>() = -s1.transpose();
	r.tail<4>() = -s1.transpose()*B1 + d.transpose()*C1;
    return r;
}

Eigen::MatrixXd Normal2_1::Body2Jacobian(void)
{
	Eigen::Vector3d s1 = caams::Ap(body1->rk_p)*s1_p;
	Eigen::Matrix3x4d B2 = 2.0*(caams::G(body2->rk_p)*caams::a_minus(s2B_p)+s2B_p*body2->rk_p.transpose());
	Eigen::Matrix<double,1,7> r;
	r.head<3>() = s1.transpose();
	r.tail<4>() = s1.transpose()*B2;
	return r;
}

Eigen::MatrixXd Normal2_1::Body1ModifiedJacobian(void)
{
	Eigen::Vector3d s1 = caams::Ap(body1->rk_p)*s1_p;
	Eigen::Vector3d d =
				body2->rk_r+caams::Ap(body2->rk_p)*s2B_p
				-body1->rk_r-caams::Ap(body1->rk_p)*s1B_p;
	Eigen::Matrix3x4d B1m = 2.0*caams::G(body1->rk_p)*caams::a_minus(s1B_p);
	Eigen::Matrix3x4d C1m = 2.0*caams::G(body1->rk_p)*caams::a_minus(s1_p);
	Eigen::Matrix<double,1,7> r;
	r.head<3>() = -s1.transpose();
	r.tail<4>() = -s1.transpose()*B1m + d.transpose()*C1m;
	return r;
}

Eigen::MatrixXd Normal2_1::Body2ModifiedJacobian(void)
{
	Eigen::Vector3d s1 = caams::Ap(body1->rk_p)*s1_p;
	Eigen::Matrix3x4d B2m = 2.0*caams::G(body2->rk_p)*caams::a_minus(s2B_p);
	Eigen::Matrix<double,1,7> r;
	r.head<3>() = s1.transpose();
	r.tail<4>() = s1.transpose()*B2m;
	return r;
}

Eigen::VectorXd Normal2_1::PHI(void)
{
	Eigen::Vector3d s1 = caams::Ap(body1->rk_p)*s1_p;
	Eigen::Vector3d d =
				body2->rk_r+caams::Ap(body2->rk_p)*s2B_p
				-body1->rk_r-caams::Ap(body1->rk_p)*s1B_p;
	return s1.transpose()*d;
}

Eigen::VectorXd Normal2_1::ModifiedGamma(void)
{
	Eigen::Vector3d s1 = caams::Ap(body1->rk_p)*s1_p;
	Eigen::Vector3d d =
				body2->rk_r+caams::Ap(body2->rk_p)*s2B_p
				-body1->rk_r-caams::Ap(body1->rk_p)*s1B_p;
	Eigen::Vector3d omega1 = 2.0*caams::G(body1->rk_p)*body1->rk_p_dot;
	Eigen::Vector3d omega2 = 2.0*caams::G(body2->rk_p)*body2->rk_p_dot;
	Eigen::Vector3d d_dot(
				body2->rk_r_dot + caams::SS(omega2)*caams::Ap(body2->rk_p)*s2B_p
				-body1->rk_r_dot - caams::SS(omega1)*caams::Ap(body1->rk_p)*s1B_p);
	Eigen::Vector3d s1_dot = caams::SS(omega1)*s1;
	return s1.transpose()*(h(body2->rk_p_dot,s2B_p)-h(body1->rk_p_dot,s1B_p))
			+ d.transpose()*h(body1->rk_p_dot,s1_p)
			- 2.0*s1_dot.transpose()*d_dot;
}

void Normal2_1::Draw(void)
{

}

void Normal2_1::DrawReaction(void)
{

}




Parallel1_2::Parallel1_2(
        Body *body1,
        Body *body2,
		Eigen::Vector3d s1_p,
		Eigen::Vector3d s2_p):
    Constraint(body1,body2,2),
    s1_p(s1_p),
    s2_p(s2_p)
{
}

Eigen::MatrixXd Parallel1_2::Body1Jacobian(void)
{
	Eigen::Vector3d s1 = caams::Ap(body1->rk_p)*s1_p;
	Eigen::Vector3d s2 = caams::Ap(body2->rk_p)*s2_p;
	Eigen::Matrix3x4d C1 = 2.0*(caams::G(body1->rk_p)*caams::a_minus(s1_p)+s1_p*body1->rk_p.transpose());
	Eigen::Matrix<double,2,7> r;
	r.block<2,3>(0,0) = Eigen::Matrix<double,2,3>::Zero();
	Eigen::Matrix3x4d PHI_p(-caams::SS(s2)*C1);
	int i_max;
	s1.cwiseAbs().maxCoeff(&i_max);
	switch(i_max){
	case 0:
		r.block<2,4>(0,3) = PHI_p.block<2,4>(1,0);
        break;
	case 1:
		r.block<1,4>(0,3) = PHI_p.block<1,4>(0,0);
		r.block<1,4>(1,3) = PHI_p.block<1,4>(2,0);
        break;
	case 2:
		r.block<2,4>(0,3) = PHI_p.block<2,4>(0,0);
        break;
    }
    return r;
}

Eigen::MatrixXd Parallel1_2::Body2Jacobian(void)
{
	Eigen::Vector3d s1(caams::Ap(body1->rk_p)*s1_p);
	Eigen::Matrix3x4d C2 = 2.0*(caams::G(body2->rk_p)*caams::a_minus(s2_p) + s2_p*body2->rk_p.transpose());
	Eigen::Matrix<double,2,7> r;
	r.block<2,3>(0,0) = Eigen::Matrix<double,2,3>::Zero();
	Eigen::Matrix3x4d PHI_p(caams::SS(s1)*C2);
	int i_max;
	s1.cwiseAbs().maxCoeff(&i_max);
	switch(i_max){
	case 0:
		r.block<2,4>(0,3) = PHI_p.block<2,4>(1,0);
		break;
	case 1:
		r.block<1,4>(0,3) = PHI_p.block<1,4>(0,0);
		r.block<1,4>(1,3) = PHI_p.block<1,4>(2,0);
		break;
	case 2:
		r.block<2,4>(0,3) = PHI_p.block<2,4>(0,0);
		break;
	}
	return r;
}

Eigen::MatrixXd Parallel1_2::Body1ModifiedJacobian(void)
{
	Eigen::Vector3d s1 = caams::Ap(body1->rk_p)*s1_p;
	Eigen::Vector3d s2 = caams::Ap(body2->rk_p)*s2_p;
	Eigen::Matrix3x4d C1m = 2.0*caams::G(body1->rk_p)*caams::a_minus(s1_p);
	Eigen::Matrix<double,2,7> r;
	r.block<2,3>(0,0) = Eigen::Matrix<double,2,3>::Zero();
	Eigen::Matrix3x4d PHI_p(-caams::SS(s2)*C1m);
	int i_max;
	s1.cwiseAbs().maxCoeff(&i_max);
	switch(i_max){
	case 0:
		r.block<2,4>(0,3) = PHI_p.block<2,4>(1,0);
		break;
	case 1:
		r.block<1,4>(0,3) = PHI_p.block<1,4>(0,0);
		r.block<1,4>(1,3) = PHI_p.block<1,4>(2,0);
		break;
	case 2:
		r.block<2,4>(0,3) = PHI_p.block<2,4>(0,0);
		break;
	}
    return r;
}

Eigen::MatrixXd Parallel1_2::Body2ModifiedJacobian(void)
{
	Eigen::Vector3d s1(caams::Ap(body1->rk_p)*s1_p);
	Eigen::Matrix3x4d C2m = 2.0*caams::G(body2->rk_p)*caams::a_minus(s2_p);
	Eigen::Matrix<double,2,7> r;
	r.block<2,3>(0,0) = Eigen::Matrix<double,2,3>::Zero();
	Eigen::Matrix3x4d PHI_p(caams::SS(s1)*C2m);
	int i_max;
	s1.cwiseAbs().maxCoeff(&i_max);
	switch(i_max){
	case 0:
		r.block<2,4>(0,3) = PHI_p.block<2,4>(1,0);
		break;
	case 1:
		r.block<1,4>(0,3) = PHI_p.block<1,4>(0,0);
		r.block<1,4>(1,3) = PHI_p.block<1,4>(2,0);
		break;
	case 2:
		r.block<2,4>(0,3) = PHI_p.block<2,4>(0,0);
		break;
	}
    return r;
}

Eigen::VectorXd Parallel1_2::PHI(void)
{
	Eigen::Vector3d s1 = caams::Ap(body1->rk_p)*s1_p;
	Eigen::Vector3d s2 = caams::Ap(body2->rk_p)*s2_p;
	Eigen::Vector3d r_all = caams::SS(s1)*s2;
	Eigen::Vector2d r;
	int i_max;
	s1.cwiseAbs().maxCoeff(&i_max);
	switch(i_max){
	case 0:
		r = r_all.tail<2>();
		break;
	case 1:
		r(0) = r_all(0);
		r(1) = r_all(2);
		break;
	case 2:
		r = r_all.head<2>();
		break;
	}
	return r;
}

Eigen::VectorXd Parallel1_2::ModifiedGamma(void)
{
	Eigen::Vector3d s1 = caams::Ap(body1->rk_p)*s1_p;
	Eigen::Vector3d s2 = caams::Ap(body2->rk_p)*s2_p;
	Eigen::Vector3d omega1 = 2.0*caams::G(body1->rk_p)*body1->rk_p_dot;
	Eigen::Vector3d omega2 = 2.0*caams::G(body2->rk_p)*body2->rk_p_dot;
	Eigen::Vector3d s1_dot = caams::SS(omega1)*s1;
	Eigen::Vector3d s2_dot = caams::SS(omega2)*s2;
	Eigen::Vector3d r_all(caams::SS(s1)*h(body2->rk_p_dot,s2_p)
                        -caams::SS(s2)*h(body1->rk_p_dot,s1_p)
                        -2.0*caams::SS(s1_dot)*s2_dot);
	Eigen::Vector2d r;
	int i_max;
	s1.cwiseAbs().maxCoeff(&i_max);
	switch(i_max){
	case 0:
		r = r_all.tail<2>();
		break;
	case 1:
		r(0) = r_all(0);
		r(1) = r_all(2);
		break;
	case 2:
		r = r_all.head<2>();
		break;
	}
	return r;
}

void Parallel1_2::Draw(void)
{

}

void Parallel1_2::DrawReaction(void)
{

}

Parallel2_2::Parallel2_2(
        Body *body1,
        Body *body2,
		Eigen::Vector3d s1_p,
		Eigen::Vector3d s1B_p,
		Eigen::Vector3d s2B_p):
    Constraint(body1,body2,2),
    s1_p(s1_p),
    s1B_p(s1B_p),
    s2B_p(s2B_p)
{
}

Eigen::MatrixXd Parallel2_2::Body1Jacobian(void)
{
	Eigen::Vector3d s1(caams::Ap(body1->rk_p)*s1_p);
	Eigen::Vector3d d(
                body2->rk_r+caams::Ap(body2->rk_p)*s2B_p
                -body1->rk_r-caams::Ap(body1->rk_p)*s1B_p);
	Eigen::Matrix<double,2,7> r;
	Eigen::Matrix<double,3,7> PHI_q;
	PHI_q.block<3,3>(0,0) = -caams::SS(s1);
	Eigen::Matrix3x4d B1(2.0*(caams::G(body1->rk_p)*caams::a_minus(s1B_p) + s1B_p*body1->rk_p.transpose()));
	Eigen::Matrix3x4d C1(2.0*(caams::G(body1->rk_p)*caams::a_minus(s1_p) + s1_p*body1->rk_p.transpose()));
	PHI_q.block<3,4>(0,3) = -caams::SS(s1)*B1 - caams::SS(d)*C1;
	int i_max;
	s1.cwiseAbs().maxCoeff(&i_max);
	switch(i_max){
	case 0:
		r.block<2,7>(0,0) = PHI_q.block<2,7>(1,0);
		break;
	case 1:
		r.block<1,7>(0,0) = PHI_q.block<1,7>(0,0);
		r.block<1,7>(1,0) = PHI_q.block<1,7>(2,0);
		break;
	case 2:
		r.block<2,7>(0,0) = PHI_q.block<2,7>(0,0);
		break;
	}
	return r;
}

Eigen::MatrixXd Parallel2_2::Body2Jacobian(void)
{
	Eigen::Vector3d s1(caams::Ap(body1->rk_p)*s1_p);
	Eigen::Matrix<double,2,7> r;
	Eigen::Matrix<double,3,7> PHI_q;
	PHI_q.block<3,3>(0,0) = caams::SS(s1);
	Eigen::Matrix3x4d B2(2.0*(caams::G(body2->rk_p)*caams::a_minus(s2B_p) + s2B_p*body2->rk_p.transpose()));
	PHI_q.block<3,4>(0,3) = caams::SS(s1)*B2;
	int i_max;
	s1.cwiseAbs().maxCoeff(&i_max);
	switch(i_max){
	case 0:
		r.block<2,7>(0,0) = PHI_q.block<2,7>(1,0);
		break;
	case 1:
		r.block<1,7>(0,0) = PHI_q.block<1,7>(0,0);
		r.block<1,7>(1,0) = PHI_q.block<1,7>(2,0);
		break;
	case 2:
		r.block<2,7>(0,0) = PHI_q.block<2,7>(0,0);
		break;
	}
    return r;
}

Eigen::MatrixXd Parallel2_2::Body1ModifiedJacobian(void)
{
	Eigen::Vector3d s1(caams::Ap(body1->rk_p)*s1_p);
	Eigen::Vector3d d(
				body2->rk_r+caams::Ap(body2->rk_p)*s2B_p
				-body1->rk_r-caams::Ap(body1->rk_p)*s1B_p);
	Eigen::Matrix<double,2,7> r;
	Eigen::Matrix<double,3,7> PHI_q;
	PHI_q.block<3,3>(0,0) = -caams::SS(s1);
	Eigen::Matrix3x4d B1m(2.0*(caams::G(body1->rk_p)*caams::a_minus(s1B_p)));
	Eigen::Matrix3x4d C1m(2.0*(caams::G(body1->rk_p)*caams::a_minus(s1_p)));
	PHI_q.block<3,4>(0,3) = -caams::SS(s1)*B1m - caams::SS(d)*C1m;
	int i_max;
	s1.cwiseAbs().maxCoeff(&i_max);
	switch(i_max){
	case 0:
		r.block<2,7>(0,0) = PHI_q.block<2,7>(1,0);
		break;
	case 1:
		r.block<1,7>(0,0) = PHI_q.block<1,7>(0,0);
		r.block<1,7>(1,0) = PHI_q.block<1,7>(2,0);
		break;
	case 2:
		r.block<2,7>(0,0) = PHI_q.block<2,7>(0,0);
		break;
	}
	return r;
}

Eigen::MatrixXd Parallel2_2::Body2ModifiedJacobian(void)
{
	Eigen::Vector3d s1(caams::Ap(body1->rk_p)*s1_p);
	Eigen::Matrix<double,2,7> r;
	Eigen::Matrix<double,3,7> PHI_q;
	PHI_q.block<3,3>(0,0) = caams::SS(s1);
	Eigen::Matrix3x4d B2m(2.0*(caams::G(body2->rk_p)*caams::a_minus(s2B_p)));
	PHI_q.block<3,4>(0,3) = caams::SS(s1)*B2m;
	int i_max;
	s1.cwiseAbs().maxCoeff(&i_max);
	switch(i_max){
	case 0:
		r.block<2,7>(0,0) = PHI_q.block<2,7>(1,0);
		break;
	case 1:
		r.block<1,7>(0,0) = PHI_q.block<1,7>(0,0);
		r.block<1,7>(1,0) = PHI_q.block<1,7>(2,0);
		break;
	case 2:
		r.block<2,7>(0,0) = PHI_q.block<2,7>(0,0);
		break;
	}
	return r;
}

Eigen::VectorXd Parallel2_2::PHI(void)
{
	Eigen::Vector3d s1 = caams::Ap(body1->rk_p)*s1_p;
	Eigen::Vector3d d(
				body2->rk_r+caams::Ap(body2->rk_p)*s2B_p
				-body1->rk_r-caams::Ap(body1->rk_p)*s1B_p);
	Eigen::Vector3d r_all = caams::SS(s1)*d;
	Eigen::Vector2d r;
	int i_max;
	s1.cwiseAbs().maxCoeff(&i_max);
	switch(i_max){
	case 0:
		r = r_all.tail<2>();
		break;
	case 1:
		r(0) = r_all(0);
		r(1) = r_all(2);
		break;
	case 2:
		r = r_all.head<2>();
		break;
	}
	return r;
}

Eigen::VectorXd Parallel2_2::ModifiedGamma(void)
{
	Eigen::Vector3d s1(caams::Ap(body1->rk_p)*s1_p);
	Eigen::Vector3d omega1(2.0*caams::G(body1->rk_p)*body1->rk_p_dot);
	Eigen::Vector3d s1_dot(caams::SS(omega1)*s1);
	Eigen::Vector3d d(
                body2->rk_r+caams::Ap(body2->rk_p)*s2B_p
                -body1->rk_r-caams::Ap(body1->rk_p)*s1B_p);
	Eigen::Vector3d omega2(2.0*caams::G(body2->rk_p)*body2->rk_p_dot);
	Eigen::Vector3d d_dot(
                body2->rk_r_dot + caams::SS(omega2)*caams::Ap(body2->rk_p)*s2B_p
                -body1->rk_r_dot - caams::SS(omega1)*caams::Ap(body1->rk_p)*s1B_p);

	Eigen::Vector3d r_all(caams::SS(s1)*(h(body2->rk_p_dot,s2B_p)-h(body1->rk_p_dot,s1B_p))
                        -caams::SS(d)*h(body1->rk_p_dot,s1_p)
                        -2.0*caams::SS(s1_dot)*d_dot);
	Eigen::Vector2d r;
	int i_max;
	s1.cwiseAbs().maxCoeff(&i_max);
	switch(i_max){
	case 0:
		r = r_all.tail<2>();
		break;
	case 1:
		r(0) = r_all(0);
		r(1) = r_all(2);
		break;
	case 2:
		r = r_all.head<2>();
		break;
	}
	return r;
}

void Parallel2_2::Draw(void)
{

}

void Parallel2_2::DrawReaction(void)
{

}

#include "constraint.h"
#include <stdio.h>
#include <iostream>
#include <cmath>
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

DependentConstraint::DependentConstraint(
        Body *body1, Body *body2, int N_eqn,
        double theta0, double theta_dot0) :
    Constraint(body1, body2, N_eqn),
    theta(theta0),
    theta0(theta0),
    theta_dot(theta_dot0)
{

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

ScrewJoint::ScrewJoint(
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
        double theta_dot0):
    DependentConstraint(body1, body2, 2, theta0, theta_dot0),
    s1_p(s1_p),
    s2_p(s2_p),
    na_p(na_p),
    n1_p(n1_p),
    n2A_p(n2A_p),
    n2B_p(n2B_p),
    alpha(alpha)
{
    // Initialize d0 from the current state
    Eigen::Vector3d d =
            body1->r+caams::Ap(body1->p)*s1_p
            -body2->r-caams::Ap(body2->p)*s2_p;
    double d2 = d.dot(d);
    d0 = sqrt(d2);
}

double primary_angle(double theta_raw)
{
    double theta = fmod(theta_raw, 2.0*M_PI);
    if(theta < 0.0) theta += 2.0*M_PI;
    if(theta > M_PI) theta -= 2.0*M_PI;
    return theta;
}

bool ScrewJoint::theta_in_A_zone(void)
{
    double ptheta = primary_angle(rk_theta);
    if((ptheta<M_PI/4 && ptheta>-M_PI/4) || ptheta>3.0*M_PI/4 || ptheta<-3.0*M_PI/4.0){
        return false;
    }else{
        return true;
    }
}

Eigen::Vector3d ScrewJoint::d_vec(void)
{
    Eigen::Vector3d r =
            body2->rk_r+caams::Ap(body2->rk_p)*s1_p
            -body1->rk_r-caams::Ap(body1->rk_p)*s2_p;
    return r;
}

Eigen::Vector3d ScrewJoint::d_dot_vec(void)
{
    Eigen::Matrix3d omega1 = caams::SS(caams::omega_p_dot(body1->rk_p, body1->rk_p_dot));
    Eigen::Matrix3d omega2 = caams::SS(caams::omega_p_dot(body2->rk_p, body2->rk_p_dot));
    Eigen::Vector3d d_dot =
            body2->rk_r_dot + omega2*s2_p - body1->rk_r_dot - omega1*s1_p;
    return d_dot;
}

Eigen::MatrixXd ScrewJoint::Body1Jacobian(void)
{
    Eigen::Vector3d d = d_vec();
    double d2 = d.dot(d);
    double D = sqrt(d2);
    Eigen::Matrix<double,2,7> r;
    r.block<1,3>(0,0) = -1.0/D*d.transpose();
    r.block<1,4>(0,3) = -2.0/D*d.transpose()*(caams::G(body1->rk_p)*caams::a_minus(s1_p) + s1_p*body1->rk_p.transpose());
    r.block<1,3>(1,0) = Eigen::RowVector3d::Zero();
    Eigen::Matrix3x4d C1 = 2.0*(caams::G(body1->rk_p)*caams::a_minus(n1_p)+n1_p*body1->rk_p.transpose());
    if(theta_in_A_zone()){
        Eigen::Vector3d n2A = caams::Ap(body2->rk_p)*n2A_p;
        r.block<1,4>(1,3) = n2A.transpose()*C1;
    }else{
        Eigen::Vector3d n2B = caams::Ap(body2->rk_p)*n2B_p;
        r.block<1,4>(1,3) = n2B.transpose()*C1;
    }
    std::cout << "ScrewJoint::Body1Jacobian r:\n"
              << r << std::endl;
    return r;
}

Eigen::MatrixXd ScrewJoint::Body2Jacobian(void)
{
    Eigen::Vector3d d = d_vec();
    double d2 = d.dot(d);
    double D = sqrt(d2);
    double phi = atan(alpha);
    double beta = cos(phi)*sin(phi);
    Eigen::Matrix<double,2,7> r;
    Eigen::RowVector3d phi_r = 1.0/D*d.transpose();
    r.block<1,3>(0,0) = phi_r;
    r.block<1,4>(0,3) = 2.0/D*d.transpose()*(caams::G(body2->rk_p)*caams::a_minus(s2_p) + s2_p*body2->rk_p.transpose());

    Eigen::Vector3d n1 = caams::Ap(body1->rk_p)*n1_p;
    Eigen::Matrix3x4d C2;
    if(theta_in_A_zone()){
        C2 = 2.0*(caams::G(body2->rk_p)*caams::a_minus(n2A_p) + n2A_p*body2->rk_p.transpose());
    }else{
        C2 = 2.0*(caams::G(body2->rk_p)*caams::a_minus(n2B_p) + n2B_p*body2->rk_p.transpose());
    }
    Eigen::RowVector4d phi_p = n1.transpose()*C2;
    r.block<1,4>(1,3) = phi_p;
    r.block<1,3>(1,0) = Eigen::RowVector3d::Zero();
    std::cout << "ScrewJoint::Body2Jacobian r:\n"
              << r << std::endl;
    return r;
}

Eigen::MatrixXd ScrewJoint::Body1ModifiedJacobian(void)
{
    Eigen::Vector3d d = d_vec();
    double d2 = d.dot(d);
    double D = sqrt(d2);
    Eigen::Matrix<double,2,7> r;
    r.block<1,3>(0,0) = -1.0/D*d.transpose();
    r.block<1,4>(0,3) = -2.0/D*d.transpose()*caams::G(body1->rk_p)*caams::a_minus(s1_p);

    r.block<1,3>(1,0) = Eigen::RowVector3d::Zero();
    Eigen::Matrix3x4d C1 = 2.0*caams::G(body1->rk_p)*caams::a_minus(n1_p);
    if(theta_in_A_zone()){
        Eigen::Vector3d n2A = caams::Ap(body2->rk_p)*n2A_p;
        r.block<1,4>(1,3) = n2A.transpose()*C1;
    }else{
        Eigen::Vector3d n2B = caams::Ap(body2->rk_p)*n2B_p;
        r.block<1,4>(1,3) = n2B.transpose()*C1;
    }
    std::cout << "ScrewJoint::Body1ModifiedJacobian r:\n"
              << r << std::endl;
    return r;
}

Eigen::MatrixXd ScrewJoint::Body2ModifiedJacobian(void)
{
    Eigen::Vector3d d = d_vec();
    double d2 = d.dot(d);
    double D = sqrt(d2);
    Eigen::Matrix<double,2,7> r;
    r.block<1,3>(0,0) = 1.0/D*d.transpose();
    r.block<1,4>(0,3) = 2.0/D*d.transpose()*caams::G(body2->rk_p)*caams::a_minus(s2_p);

    r.block<1,3>(1,0) = Eigen::RowVector3d::Zero();
    Eigen::Vector3d n1 = caams::Ap(body1->rk_p)*n1_p;

    Eigen::Matrix3x4d C2;
    if(theta_in_A_zone()){
        C2 = 2.0*caams::G(body2->rk_p)*caams::a_minus(n2A_p);
    }else{
        C2 = 2.0*caams::G(body2->rk_p)*caams::a_minus(n2B_p);
    }
    r.block<1,4>(1,3) = n1.transpose()*C2;
    std::cout << "ScrewJoint::Body2ModifiedJacobian r:\n"
              << r << std::endl;
    return r;
}

Eigen::VectorXd ScrewJoint::PHI(void)
{
    Eigen::Vector3d d = d_vec();
    double d2 = d.dot(d);
    double D = sqrt(d2);
    Eigen::Vector2d r;
    r(0) = D - d0 - alpha*(theta - theta0);
    Eigen::Vector3d n1 = caams::Ap(body1->rk_p)*n1_p;
    Eigen::Vector3d n2;
    double phi_theta;
    if(theta_in_A_zone()){
        n2 = caams::Ap(body2->rk_p)*n2A_p;
        phi_theta = cos(rk_theta);
    }else{
        n2 = caams::Ap(body2->rk_p)*n2B_p;
        phi_theta = cos(rk_theta - M_PI/2.0);
    }
    r(1) = n1.dot(n2) - phi_theta;
    std::cout << "ScrewJoint::PHI r:\n"
              << r << std::endl;
    return r;
}

Eigen::VectorXd ScrewJoint::ModifiedGamma(void)
{
    Eigen::Vector3d d = d_vec();
    double d2 = d.dot(d);
    double D = sqrt(d2);
    Eigen::Vector3d d_dot = d_dot_vec();

    Eigen::Vector2d r;
    r(0) = 2.0/D*d.transpose()*(caams::G(body1->rk_p_dot)*caams::L(body1->rk_p_dot).transpose()*s1_p
                              -caams::G(body2->rk_p_dot)*caams::L(body2->rk_p_dot).transpose()*s2_p);
    r(0) -= d_dot.dot(d_dot)/D;
    double d_dot_d_dot = d.dot(d_dot);
    r(0) += d_dot_d_dot*d_dot_d_dot/d2/D;
    Eigen::Vector3d n1 = caams::Ap(body1->rk_p)*n1_p;
    Eigen::Vector3d n1_dot = caams::SS(caams::omega_p_dot(body1->rk_p, body1->rk_p_dot))*n1;
    Eigen::Vector3d C1 = 2.0*caams::G(body1->rk_p_dot)*caams::L(body1->rk_p_dot).transpose()*n1_p;
    Eigen::Vector3d n2;
    Eigen::Vector3d n2_dot;
    Eigen::Vector3d C2;
    double theta_gamma;
    if(theta_in_A_zone()){
        n2 = caams::Ap(body2->rk_p)*n2A_p;
        n2_dot = caams::SS(caams::omega_p_dot(body2->rk_p, body2->rk_p_dot))*n2;
        C2 = 2.0*caams::G(body2->rk_p_dot)*caams::L(body2->rk_p_dot).transpose()*n2A_p;
        theta_gamma = rk_theta_dot*rk_theta_dot*cos(rk_theta);
    }else{
        n2 = caams::Ap(body2->rk_p)*n2B_p;
        n2_dot = caams::SS(caams::omega_p_dot(body2->rk_p, body2->rk_p_dot))*n2;
        C2 = 2.0*caams::G(body2->rk_p_dot)*caams::L(body2->rk_p_dot).transpose()*n2B_p;
        theta_gamma = rk_theta_dot*rk_theta_dot*cos(rk_theta - M_PI/2.0);
    }
    r(1) = -n1.dot(C2) - n2.dot(C1) - 2.0*n1_dot.dot(n2_dot) - theta_gamma;
    std::cout << "ScrewJoint::ModifiedGamma r:\n"
              << r << std::endl;

    return r;
}

void ScrewJoint::Draw(void)
{
}

void ScrewJoint::DrawReaction(void)
{
    Eigen::Matrix<double,1,7> phi2_1 = Body2Jacobian().block<1,7>(0,0);
    Eigen::Matrix<double,7,1> g2_1 = phi2_1.transpose()*lambda(0);
    Eigen::Vector3d body2_1_force = g2_1.head<3>();
    Eigen::Vector3d body2_1_torque = 0.5*caams::G(body2->p)*g2_1.tail<4>();

    std::cout << "ScrewJoint::DrawReaction\n"
              << "lambda:\n" << lambda << "\n"
              << "body2_1_force:\n"
              << body2_1_force
              << "\nbody2_1_torque:\n"
              << body2_1_torque << std::endl;

    Eigen::Vector3d body1_force, body1_torque, body2_force, body2_torque;
    Reaction(body1_force, body1_torque, body2_force, body2_torque);
    std::cout << "body2_force:\n"
              << body2_force
              << "\nbody2_torqu:\n"
              << body2_torque
              << std::endl;
}

Eigen::VectorXd ScrewJoint::ThetaModifiedJacobian(void)
{
    Eigen::Vector2d r;
    r(0) = -alpha;
    if(theta_in_A_zone()){
        r(1) = sin(rk_theta);
    }else{
        r(1) = sin(rk_theta - M_PI/2.0);
    }
    std::cout << "ScrewJoint::ThetaModifiedJacobian r:\n"
              << r << std::endl;
    return r;
}

Eigen::MatrixXd ScrewJoint::Body1AccelerationEqn(void)
{
    Eigen::Matrix<double,1,7> r;
    Eigen::Vector3d na = caams::Ap(body1->rk_p)*na_p;
    r.head<3>() = -na.transpose();
    r.tail<4>() = 2.0*alpha*na.transpose()*caams::G(body1->rk_p);
    std::cout << "ScrewJoint::Body1AccelerationEqn r:\n"
              << r << std::endl;
    return r;
}

Eigen::MatrixXd ScrewJoint::Body2AccelerationEqn(void)
{
    Eigen::Matrix<double,1,7> r;
    Eigen::Vector3d na = caams::Ap(body1->rk_p)*na_p;
    r.head<3>() = na.transpose();
    r.tail<4>() = -2.0*alpha*na.transpose()*caams::G(body2->rk_p);
    std::cout << "ScrewJoint::Body2AccelerationEqn r:\n"
              << r << std::endl;
    return r;
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
    r.tail<4>() = 2.0*s2.transpose()*caams::G(body1->rk_p)*caams::a_minus(s1_p);
    return r;
}

Eigen::MatrixXd Normal1_1::Body2ModifiedJacobian(void)
{
	Eigen::Vector3d s1 = caams::Ap(body1->rk_p)*s1_p;
	Eigen::Matrix<double,1,7> r;
	r.head<3>() = Eigen::RowVector3d::Zero();
    r.tail<4>() = 2.0*s1.transpose()*caams::G(body2->rk_p)*caams::a_minus(s2_p);
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
    std::cout << "Parallel1_2::Body1Jacobian r:\n"
              << r << std::endl;
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
    std::cout << "Parallel1_2::Body2Jacobian r:\n"
              << r << std::endl;
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
    std::cout << "Parallel1_2::Body1ModifiedJacobian r:\n"
              << r << std::endl;
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
    std::cout << "Parallel1_2::Body2ModifiedJacobian r:\n"
              << r << std::endl;
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
    std::cout << "Parallel1_2::PHI r:\n"
              << r << std::endl;
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
    std::cout << "Parallel1_2::ModifiedGamma r:\n"
              << r << std::endl;
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
    std::cout << "Parallel2_2::Body1Jacobian r:\n"
              << r << std::endl;
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
    std::cout << "Parallel2_2::Body2Jacobian r:\n"
              << r << std::endl;
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
    std::cout << "Parallel2_2::Body1ModifiedJacobian r:\n"
              << r << std::endl;
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
    std::cout << "Parallel2_2::Body2ModifiedJacobian r:\n"
              << r << std::endl;
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
    std::cout << "Parallel2_2::PHI r:\n"
              << r << std::endl;
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
    std::cout << "Parallel2_2::ModifiedGamma r:\n"
              << r << std::endl;
    return r;
}

void Parallel2_2::Draw(void)
{

}

void Parallel2_2::DrawReaction(void)
{

}

PlanarJoint::PlanarJoint(
        Body *body1,
        Body *body2,
        Eigen::Vector3d s1_p,
        Eigen::Vector3d s2_p,
        Eigen::Vector3d u2_p):
    Constraint(body1,body2,1),
    s1_p(s1_p),
    s2_p(s2_p),
    u2_p(u2_p)
{
}

Eigen::MatrixXd PlanarJoint::Body1Jacobian(void)
{
    Eigen::Vector3d u2 = caams::Ap(body1->rk_p)*u2_p;
    Eigen::Matrix<double,1,7> r;
    r.head<3>() = u2.transpose();
    r.tail<4>() = 2.0*u2.transpose()*(caams::G(body1->rk_p)*caams::a_minus(s1_p) + s1_p*body1->rk_p.transpose());
    return r;
}

Eigen::MatrixXd PlanarJoint::Body2Jacobian(void)
{
    Eigen::Vector3d s1 = body1->rk_r + caams::Ap(body1->rk_p)*s1_p;
    Eigen::Vector3d s2 = body2->rk_r + caams::Ap(body2->rk_p)*s2_p;
    Eigen::Vector3d u2 = caams::Ap(body2->rk_p)*u2_p;
    Eigen::Matrix<double,1,7> r;
    r.head<3>() = -u2.transpose();
    r.tail<4>() =
            2*(s1-s2).transpose()*(caams::G(body2->rk_p)*caams::a_minus(u2_p) + u2_p*body2->rk_p.transpose())
            -2*u2.transpose()*(caams::G(body2->rk_p)*caams::a_minus(s2_p) + s2_p*body2->rk_p.transpose());
    return r;
}

Eigen::MatrixXd PlanarJoint::Body1ModifiedJacobian(void)
{
    Eigen::Vector3d u2 = caams::Ap(body1->rk_p)*u2_p;
    Eigen::Matrix<double,1,7> r;
    r.head<3>() = u2.transpose();
    r.tail<4>() = 2.0*u2.transpose()*caams::G(body1->rk_p)*caams::a_minus(s1_p);
    return r;
}

Eigen::MatrixXd PlanarJoint::Body2ModifiedJacobian(void)
{
    Eigen::Vector3d s1 = body1->rk_r + caams::Ap(body1->rk_p)*s1_p;
    Eigen::Vector3d s2 = body2->rk_r + caams::Ap(body2->rk_p)*s2_p;
    Eigen::Vector3d u2 = caams::Ap(body2->rk_p)*u2_p;
    Eigen::Matrix<double,1,7> r;
    r.head<3>() = -u2.transpose();
    r.tail<4>() =
            2*(s1-s2).transpose()*caams::G(body2->rk_p)*caams::a_minus(u2_p)
            -2*u2.transpose()*caams::G(body2->rk_p)*caams::a_minus(s2_p);
    return r;
}

Eigen::VectorXd PlanarJoint::PHI(void)
{
    Eigen::Vector3d s1 = body1->rk_r + caams::Ap(body1->rk_p)*s1_p;
    Eigen::Vector3d s2 = body2->rk_r + caams::Ap(body2->rk_p)*s2_p;
    Eigen::Vector3d u2 = caams::Ap(body2->rk_p)*u2_p;
    Eigen::VectorXd phi = (s1-s2).transpose()*u2;
    return phi;
}

Eigen::VectorXd PlanarJoint::ModifiedGamma(void)
{
    Eigen::Vector3d s1 = body1->rk_r + caams::Ap(body1->rk_p)*s1_p;
    Eigen::Vector3d s2 = body2->rk_r + caams::Ap(body2->rk_p)*s2_p;
    Eigen::Vector3d u2 = caams::Ap(body2->rk_p)*u2_p;
    Eigen::Matrix3d A1_dot = caams::A_dot(body1->rk_p, body1->rk_p_dot);
    Eigen::Matrix3d A2_dot = caams::A_dot(body2->rk_p, body2->rk_p_dot);
    Eigen::Vector3d s1_dot = body1->rk_r_dot + A1_dot*s1_p;
    Eigen::Vector3d s2_dot = body2->rk_r_dot + A2_dot*s2_p;
    Eigen::Vector3d u2_dot = A2_dot*u2_p;
    Eigen::VectorXd gamma =
            -2*u2.transpose()*(caams::G(body1->rk_p_dot)*caams::L(body1->rk_p_dot).transpose()*s1_p
                               -caams::G(body2->rk_p_dot)*caams::L(body2->rk_p_dot).transpose()*s2_p)
            -2*(s1-s2).transpose()*caams::G(body2->rk_p_dot)*caams::L(body2->rk_p_dot).transpose()*u2_p
            -2*(s1_dot - s2_dot).transpose()*u2_dot;
    return gamma;
}

void PlanarJoint::Draw(void)
{
}

void PlanarJoint::DrawReaction(void)
{
}

//
// Composite constraints
//

UniversalJoint::UniversalJoint(
        Body *body1,
        Body *body2,
        Eigen::Vector3d s1_p, // center of joint in body1 space
        Eigen::Vector3d s2_p, // center of joint in body2 space
        Eigen::Vector3d h1_p, // rotation acis in bocy1 space
        Eigen::Vector3d h2_p): // rotation axis in body2 space
    Constraint(body1, body2, 4),
    sphericalJoint(body1, body2, s1_p, s2_p),
    normalJoint(body1, body2, h1_p, h2_p),
    s1_p(s1_p),
    s2_p(s2_p),
    h1_p(h1_p),
    h2_p(h2_p)
{

}

Eigen::MatrixXd UniversalJoint::Body1Jacobian(void)
{
    Eigen::Matrix<double,4,7> r;
    r.block<3,7>(0,0) = sphericalJoint.Body1Jacobian();
    r.block<1,7>(3,0) = normalJoint.Body1Jacobian();
    return r;
}

Eigen::MatrixXd UniversalJoint::Body2Jacobian(void)
{
    Eigen::Matrix<double,4,7> r;
    r.block<3,7>(0,0) = sphericalJoint.Body2Jacobian();
    r.block<1,7>(3,0) = normalJoint.Body2Jacobian();
    return r;
}

Eigen::MatrixXd UniversalJoint::Body1ModifiedJacobian(void)
{
    Eigen::Matrix<double,4,7> r;
    r.block<3,7>(0,0) = sphericalJoint.Body1ModifiedJacobian();
    r.block<1,7>(3,0) = normalJoint.Body1ModifiedJacobian();
    return r;
}

Eigen::MatrixXd UniversalJoint::Body2ModifiedJacobian(void)
{
    Eigen::Matrix<double,4,7> r;
    r.block<3,7>(0,0) = sphericalJoint.Body2ModifiedJacobian();
    r.block<1,7>(3,0) = normalJoint.Body2ModifiedJacobian();
    return r;
}

Eigen::VectorXd UniversalJoint::PHI(void)
{
    Eigen::Vector4d r;
    r.head<3>() = sphericalJoint.PHI();
    r.tail<1>() = normalJoint.PHI();
    return r;
}

Eigen::VectorXd UniversalJoint::ModifiedGamma(void)
{
    Eigen::Vector4d r;
    r.head<3>() = sphericalJoint.ModifiedGamma();
    r.tail<1>() = normalJoint.ModifiedGamma();
    return r;
}

void UniversalJoint::Draw(void)
{
    Eigen::Matrix3d A_body1 = caams::Ap(body1->p);
    Eigen::Matrix3d A_body2 = caams::Ap(body2->p);
    Eigen::Vector3d s1 = body1->r + A_body1*s1_p;
    Eigen::Vector3d s2 = body2->r + A_body2*s2_p;

    glBegin(GL_LINES);
        glColor3f(1.0f,1.0f,0.0f);
        glVertex3dv(body1->r.data());
        glVertex3dv(s1.data());
        glColor3f(0.0f,1.0f,1.0f);
        glVertex3dv(body2->r.data());
        glVertex3dv(s2.data());
    glEnd();

    Eigen::Vector3d h1 = A_body1*h1_p;
    Eigen::Vector3d h2 = A_body2*h2_p;
    Eigen::Vector3d a1p = s1 + h1;
    Eigen::Vector3d a1m = s1 - h1;
    Eigen::Vector3d a2p = s2 + h2;
    Eigen::Vector3d a2m = s2 - h2;

    glBegin(GL_LINES);
        glColor3d(0.9, 0.3, 0.3);
        glVertex3dv(a1p.data());
        glVertex3dv(a1m.data());
        glColor3d(0.3, 0.3, 0.9);
        glVertex3dv(a2p.data());
        glVertex3dv(a2m.data());
    glEnd();
}

void UniversalJoint::DrawReaction(void)
{
}

RevoluteJoint::RevoluteJoint(
        Body *body1,
        Body *body2,
        Eigen::Vector3d s1_p, // center of joint in body1 space
        Eigen::Vector3d s2_p, // center of joint in body2 space
        Eigen::Vector3d h1_p, // rotation acis in bocy1 space
        Eigen::Vector3d h2_p): // rotation axis in body2 space
    Constraint(body1, body2, 5),
    sphericalJoint(body1, body2, s1_p, s2_p),
    parallelJoint(body1, body2, h1_p, h2_p),
    s1_p(s1_p),
    s2_p(s2_p),
    h1_p(h1_p),
    h2_p(h2_p)
{

}

Eigen::MatrixXd RevoluteJoint::Body1Jacobian(void)
{
    Eigen::Matrix<double,5,7> r;
    r.block<3,7>(0,0) = sphericalJoint.Body1Jacobian();
    r.block<2,7>(3,0) = parallelJoint.Body1Jacobian();
    return r;
}

Eigen::MatrixXd RevoluteJoint::Body2Jacobian(void)
{
    Eigen::Matrix<double,5,7> r;
    r.block<3,7>(0,0) = sphericalJoint.Body2Jacobian();
    r.block<2,7>(3,0) = parallelJoint.Body2Jacobian();
    return r;
}

Eigen::MatrixXd RevoluteJoint::Body1ModifiedJacobian(void)
{
    Eigen::Matrix<double,5,7> r;
    r.block<3,7>(0,0) = sphericalJoint.Body1ModifiedJacobian();
    r.block<2,7>(3,0) = parallelJoint.Body1ModifiedJacobian();
    return r;
}

Eigen::MatrixXd RevoluteJoint::Body2ModifiedJacobian(void)
{
    Eigen::Matrix<double,5,7> r;
    r.block<3,7>(0,0) = sphericalJoint.Body2ModifiedJacobian();
    r.block<2,7>(3,0) = parallelJoint.Body2ModifiedJacobian();
    return r;
}

Eigen::VectorXd RevoluteJoint::PHI(void)
{
    Eigen::Matrix<double,5,1> r;
    r.head<3>() = sphericalJoint.PHI();
    r.tail<2>() = parallelJoint.PHI();
    return r;
}

Eigen::VectorXd RevoluteJoint::ModifiedGamma(void)
{
    Eigen::Matrix<double,5,1> r;
    r.head<3>() = sphericalJoint.ModifiedGamma();
    r.tail<2>() = parallelJoint.ModifiedGamma();
    return r;
}

void RevoluteJoint::Draw(void)
{
    sphericalJoint.Draw();
}

void RevoluteJoint::DrawReaction(void)
{
}


TranslationalJoint::TranslationalJoint(
        Body *body1,
        Body *body2,
        Eigen::Vector3d a1_p, // rotation axis in body1 space
        Eigen::Vector3d a2_p, // rotation axis in body2 space
        Eigen::Vector3d s1B_p, // point on axis in body1 space
        Eigen::Vector3d s2B_p, // point on axis in body2 space
        Eigen::Vector3d h1_p, // normal axis in body1 space
        Eigen::Vector3d h2_p): // normal axis in body2 space
    Constraint(body1, body2, 5),
    parallel1Joint(body1, body2, a1_p, a2_p),
    parallel2Joint(body1, body2, a1_p, s1B_p, s2B_p),
    normalJoint(body1, body2, h1_p, h2_p),
    a1_p(a1_p),
    a2_p(a2_p),
    s1B_p(s1B_p),
    s2B_p(s2B_p),
    h1_p(h1_p),
    h2_p(h2_p)
{
}

Eigen::MatrixXd TranslationalJoint::Body1Jacobian(void)
{
    Eigen::Matrix<double,5,7> r;
    r.block<2,7>(0,0) = parallel1Joint.Body1Jacobian();
    r.block<2,7>(2,0) = parallel2Joint.Body1Jacobian();
    r.block<1,7>(4,0) = normalJoint.Body1Jacobian();
    return r;
}

Eigen::MatrixXd TranslationalJoint::Body2Jacobian(void)
{
    Eigen::Matrix<double,5,7> r;
    r.block<2,7>(0,0) = parallel1Joint.Body2Jacobian();
    r.block<2,7>(2,0) = parallel2Joint.Body2Jacobian();
    r.block<1,7>(4,0) = normalJoint.Body2Jacobian();
    return r;
}

Eigen::MatrixXd TranslationalJoint::Body1ModifiedJacobian(void)
{
    Eigen::Matrix<double,5,7> r;
    r.block<2,7>(0,0) = parallel1Joint.Body1ModifiedJacobian();
    r.block<2,7>(2,0) = parallel2Joint.Body1ModifiedJacobian();
    r.block<1,7>(4,0) = normalJoint.Body1ModifiedJacobian();
    return r;
}

Eigen::MatrixXd TranslationalJoint::Body2ModifiedJacobian(void)
{
    Eigen::Matrix<double,5,7> r;
    r.block<2,7>(0,0) = parallel1Joint.Body2ModifiedJacobian();
    r.block<2,7>(2,0) = parallel2Joint.Body2ModifiedJacobian();
    r.block<1,7>(4,0) = normalJoint.Body2ModifiedJacobian();
    return r;
}

Eigen::VectorXd TranslationalJoint::PHI(void)
{
    Eigen::Matrix<double,5,1> r;
    r.segment<2>(0) = parallel1Joint.PHI();
    r.segment<2>(2) = parallel2Joint.PHI();
    r.segment<1>(4) = normalJoint.PHI();
    return r;
}

Eigen::VectorXd TranslationalJoint::ModifiedGamma(void)
{
    Eigen::Matrix<double,5,1> r;
    r.segment<2>(0) = parallel1Joint.ModifiedGamma();
    r.segment<2>(2) = parallel2Joint.ModifiedGamma();
    r.segment<1>(4) = normalJoint.ModifiedGamma();
    return r;
}

void TranslationalJoint::Draw(void)
{

}

void TranslationalJoint::DrawReaction(void)
{

}

CylindricalJoint::CylindricalJoint(
        Body *body1,
        Body *body2,
        Eigen::Vector3d a1_p, // rotation axis in body1 space
        Eigen::Vector3d a2_p, // rotation axis in body2 space
        Eigen::Vector3d s1B_p, // point on axis in body1 space
        Eigen::Vector3d s2B_p): // point on axis in body2 space
    Constraint(body1, body2, 4),
    parallel1Joint(body1, body2, a1_p, a2_p),
    parallel2Joint(body1, body2, a1_p, s1B_p, s2B_p),
    a1_p(a1_p),
    a2_p(a2_p),
    s1B_p(s1B_p),
    s2B_p(s2B_p)
{
}

Eigen::MatrixXd CylindricalJoint::Body1Jacobian(void)
{
    Eigen::Matrix<double,4,7> r;
    r.block<2,7>(0,0) = parallel1Joint.Body1Jacobian();
    r.block<2,7>(2,0) = parallel2Joint.Body1Jacobian();
    return r;
}

Eigen::MatrixXd CylindricalJoint::Body2Jacobian(void)
{
    Eigen::Matrix<double,4,7> r;
    r.block<2,7>(0,0) = parallel1Joint.Body2Jacobian();
    r.block<2,7>(2,0) = parallel2Joint.Body2Jacobian();
    return r;
}

Eigen::MatrixXd CylindricalJoint::Body1ModifiedJacobian(void)
{
    Eigen::Matrix<double,4,7> r;
    r.block<2,7>(0,0) = parallel1Joint.Body1ModifiedJacobian();
    r.block<2,7>(2,0) = parallel2Joint.Body1ModifiedJacobian();
    return r;
}

Eigen::MatrixXd CylindricalJoint::Body2ModifiedJacobian(void)
{
    Eigen::Matrix<double,4,7> r;
    r.block<2,7>(0,0) = parallel1Joint.Body2ModifiedJacobian();
    r.block<2,7>(2,0) = parallel2Joint.Body2ModifiedJacobian();
    return r;
}

Eigen::VectorXd CylindricalJoint::PHI(void)
{
    Eigen::Matrix<double,4,1> r;
    r.segment<2>(0) = parallel1Joint.PHI();
    r.segment<2>(2) = parallel2Joint.PHI();
    return r;
}

Eigen::VectorXd CylindricalJoint::ModifiedGamma(void)
{
    Eigen::Matrix<double,4,1> r;
    r.segment<2>(0) = parallel1Joint.ModifiedGamma();
    r.segment<2>(2) = parallel2Joint.ModifiedGamma();
    return r;
}

void CylindricalJoint::Draw(void)
{
    parallel2Joint.Draw();
}

void CylindricalJoint::DrawReaction(void)
{

}


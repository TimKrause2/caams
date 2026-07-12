#include "mbsystem.h"
#include "forces.h"
#include <math.h>
#include <iostream>
#include <GL/gl.h>
#include <GL/freeglut.h>
#define GLM_FORCE_RADIANS
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/transform.hpp>

#define G_GRAVITY 9.81
#define AXIS_MASS 1.0
#define AXIS_RADIUS 0.2
#define AXIS_LENGTH 5.0
#define SLIDER_MASS 0.5
#define SLIDER_RADIUS 0.5
#define SLIDER_LENGTH 0.2
#define INITIAL_OFFSET 0.75
#define REST_DROP 1.0
#define K_SPRING (SLIDER_MASS*G_GRAVITY/REST_DROP)
#define OMEGA_SPRING sqrt(K_SPRING/SLIDER_MASS)
#define D_SPRING INITIAL_OFFSET
#define ZETA_DAMP 0.01
#define K_DAMP (2*ZETA_DAMP*OMEGA_SPRING)
#define N_ROTATION 2.0
#define BETA_SCREW ((N_ROTATION*2*M_PI)/REST_DROP)

Body *datumBody = new DatumBody;
Body *axis;
Body *slider;
Constraint *revoluteJoint;
Constraint *cylindricalJoint;
Constraint *screwJoint;
ForceElement *gravity;
ForceElement *spring;

System mbsystem;

void init_system(void)
{
    // Revolute joint body1:datum body2:axis
    Eigen::Vector3d s1_p_rev(0.0, 0.0, 0.0);
    Eigen::Vector3d s2_p_rev(0.0, 0.0, 0.0);
    Eigen::Vector3d h1_p_rev(0.0, 0.0, 1.0);
    Eigen::Vector3d h2_p_rev(0.0, 0.0, 1.0);

    // Cylindrical joint body1:axis body2:slider
    Eigen::Vector3d a1_p_cyl(-1.0, 0.0, 0.0);
    Eigen::Vector3d a2_p_cyl(-1.0, 0.0, 0.0);
    Eigen::Vector3d s1B_p_cyl(-AXIS_LENGTH/2, 0.25, 0.0);
    Eigen::Vector3d s2B_p_cyl(0.0,-0.2, 0.0);

    // Screw joint body1:axis body2:slider
    Eigen::Vector3d s1_p_screw = s1B_p_cyl;
    Eigen::Vector3d s2_p_screw = s2B_p_cyl;
    Eigen::Vector3d u1_p_screw(0.0, 1.0, 0.0);
    Eigen::Vector3d u2_p_screw(0.0, 0.0, 1.0);
    Eigen::Vector3d a1_p_screw = a1_p_cyl;

    // Inititialize the postions of the bodies
    Eigen::Vector4d p_axis(1.0, 0.0, 0.0, 0.0);
    Eigen::Vector4d p_dot_axis = Eigen::Vector4d::Zero();
    Eigen::Vector4d p_slider = p_axis;
    Eigen::Vector4d p_dot_slider = Eigen::Vector4d::Zero();
    Eigen::Matrix3d A = caams::Ap(p_axis);
    Eigen::Vector3d r_axis(0.0, 0.0, 0.0);
    Eigen::Vector3d r_dot_axis = Eigen::Vector3d::Zero();
    Eigen::Vector3d r_slider =
            r_axis
            + A*s1B_p_cyl
            + A*Eigen::Vector3d(INITIAL_OFFSET+AXIS_LENGTH/2, 0.0, 0.0)
            - A*s2B_p_cyl;
    Eigen::Vector3d r_dot_slider = Eigen::Vector3d::Zero();

    axis = new CylinderXaxis(
                r_axis,
                p_axis,
                r_dot_axis,
                p_dot_axis,
                AXIS_MASS,
                AXIS_RADIUS,
                AXIS_LENGTH);

    slider = new CylinderXaxis(
                r_slider,
                p_slider,
                r_dot_slider,
                p_dot_slider,
                SLIDER_MASS,
                SLIDER_RADIUS,
                SLIDER_LENGTH);

    revoluteJoint = new RevoluteJoint(
                datumBody, axis,
                s1_p_rev, s2_p_rev,
                h1_p_rev, h2_p_rev);

    cylindricalJoint = new CylindricalJoint(
                axis, slider,
                a1_p_cyl, a2_p_cyl,
                s1B_p_cyl, s2B_p_cyl);

    screwJoint = new ScrewJointLocal(
                axis, slider,
                s1_p_screw, s2_p_screw,
                u1_p_screw, u2_p_screw,
                a1_p_screw,
                BETA_SCREW);

	gravity = new SystemGravityForce(
                Eigen::Vector3d(0.0, -G_GRAVITY, 0.0),
				mbsystem);

    spring = new LinearSpringDamperForce(
                axis, slider, s1B_p_cyl, s2B_p_cyl,
                D_SPRING+AXIS_LENGTH/2, K_SPRING, K_DAMP);

    mbsystem.AddBody(axis);
    mbsystem.AddBody(slider);
    mbsystem.AddConstraint(revoluteJoint);
    mbsystem.AddConstraint(cylindricalJoint);
    mbsystem.AddConstraint(screwJoint);
	mbsystem.AddForce(gravity);
    mbsystem.AddForce(spring);
    mbsystem.InitializeSolver();
}

glm::dmat4 camera_rotation(1.0);
glm::dmat4 camera_translation(glm::translate(glm::dvec3(0.0,0.0,3.0)));

#define NEAR_PLANE 0.1
#define FAR_PLANE 30.0
#define FOV (M_PI*90.0/180.0)

void display(void){
    int window_width;
    int window_height;

    window_width = glutGet(GLUT_WINDOW_WIDTH);
    window_height = glutGet(GLUT_WINDOW_HEIGHT);

    glm::dmat4x4 m_projection = glm::perspective(
        FOV,(double)window_width/window_height,NEAR_PLANE,FAR_PLANE);
    glMatrixMode(GL_PROJECTION);
    glLoadMatrixd(glm::value_ptr(m_projection));

    glm::dmat4 camera_to_world = camera_translation*camera_rotation;
    glm::dmat4 world_to_camera = glm::inverse(camera_to_world);

    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixd(glm::value_ptr(world_to_camera));

    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glColor3f(1.0f,1.0f,1.0f);

    mbsystem.Draw();

    datumBody->Draw();

    mbsystem.Integrate(1/60.0);

    glutSwapBuffers();

}

#define TIMER_INTERVAL (1000.0/60.0)
double timer_interval=0.0;


void timerFunc( int value )
{
    glutPostRedisplay();
    timer_interval += TIMER_INTERVAL;
    int ti_int = (int)floor(timer_interval);
    glutTimerFunc( ti_int, timerFunc, 0 );
    timer_interval -= ti_int;
}

void keyboardFunc(unsigned char key, int x, int y)
{
    switch(key){
        case ' ':
            break;
        case 27:
            glutLeaveMainLoop();
            break;
        default:
            break;
    }
}

int main(int argc, char **argv){
    glutInit(&argc,argv);
    init_system();
    glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );
    glutInitWindowSize (1280, 720);
    glutInitWindowPosition (100, 100);
    glutCreateWindow ("mbtest");
    glutDisplayFunc(display);
    glutTimerFunc( TIMER_INTERVAL, timerFunc, 0 );
    glutKeyboardFunc(keyboardFunc);
    glClearColor(0.0,0.0,0.0,0.0);
    glDrawBuffer(GL_BACK);
    glColor3d(1.0,1.0,1.0);
    glutMainLoop();
    return 0;
}

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
#define MASS 1.0
#define INITIAL_DROP 0.5
#define RADIUS 0.5
#define LENGTH 0.2
#define REST_DROP 1.0
#define K_SPRING (MASS*G_GRAVITY/REST_DROP)
#define OMEGA_SPRING sqrt(K_SPRING/MASS)
#define D_SPRING INITIAL_DROP
#define K_CRITICAL (2*sqrt(MASS*K_SPRING))
#define ZETA_DAMP 0.1
#define K_DAMP (2*ZETA_DAMP*OMEGA_SPRING)
#define N_ROTATION 10.0
#define BETA_SCREW ((N_ROTATION*2*M_PI)/REST_DROP)


Body *datumBody = new DatumBody;
Body *massBody;
Constraint *cylindricalJoint;
Constraint *screwJoint;
ForceElement *gravity;
ForceElement *spring;

System mbsystem;

void init_system(void)
{
    Eigen::Vector3d a1_p(0.0, 1.0, 0.0); // datum space
    Eigen::Vector3d a2_p(0.0, 1.0, 0.0);
    Eigen::Vector3d s1B_p(0.0, INITIAL_DROP, 0.0);
    Eigen::Vector3d s2B_p(0.0, 0.0, 0.0);
    Eigen::Vector3d r_mass(0.0, 0.0, 0.0);
    Eigen::Vector4d p_mass(1.0, 0.0, 0.0, 0.0);
    Eigen::Vector3d r_dot_mass(0.0, 0.0, 0.0);
    Eigen::Vector4d p_dot_mass(0.0, 0.0, 0.0, 0.0);

    massBody = new CylinderYaxis(
                r_mass,
                p_mass,
                r_dot_mass,
                p_dot_mass,
                MASS,
                RADIUS,
                LENGTH);

    cylindricalJoint = new CylindricalJoint(
                datumBody, massBody, a1_p, a2_p, s1B_p, s2B_p);

    Eigen::Vector3d s1_p(0.0, 0.0, 0.0);
    Eigen::Vector3d s2_p(0.0, 0.0, 0.0);
    Eigen::Vector3d u1_p(1.0, 0.0, 0.0);
    Eigen::Vector3d u2_p(0.0, 0.0,-1.0);

    screwJoint = new ScrewJoint_1(
                datumBody, massBody,
                s1B_p, s2B_p, u1_p, u2_p, a1_p,
                BETA_SCREW);

	gravity = new SystemGravityForce(
                Eigen::Vector3d(0.0, -G_GRAVITY, 0.0),
				mbsystem);

    spring = new LinearSpringDamperForce(
                datumBody, massBody, s1B_p, s2B_p,
                D_SPRING, K_SPRING, K_DAMP);

    mbsystem.AddBody(massBody);
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
    mbsystem.Integrate(1/60.0);

    mbsystem.Draw();

    datumBody->Draw();

    screwJoint->DrawReaction();

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

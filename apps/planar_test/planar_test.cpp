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
#define RADIUS 0.5
#define SPRING_LENGTH 1.0
#define REST_LENGTH 0.8
#define THETA_INCLINATION (45.0/180.0*M_PI)
#define K_SPRING (MASS*G_GRAVITY*sin(THETA_INCLINATION)/REST_LENGTH)
#define K_DAMP 0.0001


Body *datumBody = new DatumBody;
Body *massBody;
Constraint *planarJoint;
ForceElement *spring;
ForceElement *gravity;

System mbsystem;

void init_system(void)
{
    Eigen::Vector3d s1_p_plane(0.0, 0.0, 0.0);
    Eigen::Vector3d s2_p_plane(0.0, 0.0, 0.0);
    Eigen::Matrix3d A_frame = caams::AAA(THETA_INCLINATION, Eigen::Vector3d(0.0, 0.0, 1.0));
    Eigen::Vector3d u2_p_plane = A_frame.col(1);
    Eigen::Vector3d s1_p_spring(0.0, 0.0, 0.0); // center of the mass(body1)
    Eigen::Vector3d s2_p_spring = A_frame*Eigen::Vector3d(SPRING_LENGTH, 0.0, SPRING_LENGTH);

    Eigen::Vector3d r_mass(0.0, 0.0, 0.0);
    Eigen::Vector4d p_mass(1.0, 0.0, 0.0, 0.0);
    Eigen::Vector3d r_dot_mass(0.0, 0.0, 0.0);
    Eigen::Vector4d p_dot_mass(0.0, 0.0, 0.0, 0.0);

    massBody = new Sphere(
                r_mass,
                p_mass,
                r_dot_mass,
                p_dot_mass,
                MASS,
                RADIUS);

    planarJoint = new PlanarJoint(
                massBody, datumBody,
                s1_p_plane, s2_p_plane, u2_p_plane);

	gravity = new SystemGravityForce(
				Eigen::Vector3d(0.0, -9.81, 0.0),
				mbsystem);

    spring = new LinearSpringDamperForce(
                massBody, datumBody,
                s1_p_spring, s2_p_spring,
                SPRING_LENGTH, K_SPRING, K_DAMP);


    mbsystem.AddBody(massBody);
    mbsystem.AddConstraint(planarJoint);
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

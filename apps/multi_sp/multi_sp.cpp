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

#define N_CONSTRAINTS 1
#define SPHERE_RADIUS 0.5
#define SPHERE_MASS 1.0
#define THETA_Z (0.0/180.0*M_PI)
#define THETA_CONE (45.0/180.0*M_PI)
#define D_SPHERE 1.0
#define D_ARM 2.0

Body *datumBody = new DatumBody;
Body *sphereBody;
Constraint *sphericalJoints[N_CONSTRAINTS];
ForceElement *gravity;

System mbsystem;

void init_system(void)
{
    // Initialize the arms in local coordinate space
    // The rotational axis is the y-axis
    // The 'cone' of the arms points toward the positive y-axis
    Eigen::Vector3d s1_p[N_CONSTRAINTS];
    for(int i=0;i<N_CONSTRAINTS;i++){
        Eigen::Matrix3d A_rot_z = caams::AAA(THETA_CONE,
                                             Eigen::Vector3d(0.0, 0.0, 1.0));
        double theta_y = (double)i*2.0*M_PI/N_CONSTRAINTS;
        Eigen::Matrix3d A_rot_y = caams::AAA(theta_y,
                                             Eigen::Vector3d(0.0, 1.0, 0.0));
        Eigen::Vector3d n_arm = A_rot_y*A_rot_z*Eigen::Vector3d(0.0, -1.0, 0.0);
        s1_p[i] = n_arm*D_ARM;
    }

    // Initialize the orientation and location of the sphere
    Eigen::Vector4d p_sphere = caams::pAA(THETA_Z,Eigen::Vector3d(0.0, 0.0, -1.0));
    Eigen::Matrix3d A_sphere = caams::Ap(p_sphere);
    Eigen::Vector3d r_sphere = A_sphere*Eigen::Vector3d(0.0, D_SPHERE, 0.0);
    Eigen::Vector3d r_dot_sphere = Eigen::Vector3d::Zero();
    Eigen::Vector4d p_dot_sphere = Eigen::Vector4d::Zero();

    // Initialize the constraint locations in global coordinates
    Eigen::Vector3d s2_p[N_CONSTRAINTS];
    for(int i=0;i<N_CONSTRAINTS;i++){
        s2_p[i] = r_sphere + A_sphere*s1_p[i];
    }

    sphereBody = new Sphere(
                r_sphere,
                p_sphere,
                r_dot_sphere,
                p_dot_sphere,
                SPHERE_MASS,
                SPHERE_RADIUS);

    for(int i=0;i<N_CONSTRAINTS;i++){
        sphericalJoints[i] = new SphericalJoint(sphereBody,datumBody,s1_p[i],s2_p[i]);
    }

	gravity = new SystemGravityForce(
				Eigen::Vector3d(0.0, -9.81, 0.0),
				mbsystem);

    mbsystem.AddBody(sphereBody);
    for(int i=0;i<N_CONSTRAINTS;i++){
        mbsystem.AddConstraint(sphericalJoints[i]);
    }
	mbsystem.AddForce(gravity);
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

    for(int i=0;i<N_CONSTRAINTS;i++){
        sphericalJoints[i]->DrawReaction();
    }

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

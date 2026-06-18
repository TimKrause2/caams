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

#define CAMERA_LATITUDE (-15.0/180.0*M_PI)
#define SHAFT_LENGTH 1.0
#define SHAFT_RADIUS 0.2
#define SHAFT_MASS 1.0
#define SHAFT_BEND (30.0/180.0*M_PI)
#define OMEGA_MAX (0.5*2.0*M_PI)
#define TORQUE_MAX 2.0
#define K_R (TORQUE_MAX/OMEGA_MAX)


Body *shaft1Body;
Body *shaft2Body;
Body *datumBody;
BodyLocalTorque *shaft1LocalTorque;
ForceElement *shaft2Friction;
Constraint *shaft1RevoluteJoint;
Constraint *shaft2ParallelJoint;
Constraint *universalJoint;

System mbsystem;

void init_system(void)
{
    // Initialize the orientations of the shafts
    Eigen::Vector4d p_shaft1(1.0, 0.0, 0.0, 0.0);
    Eigen::Vector4d p_shaft2 = caams::pAA(SHAFT_BEND, Eigen::Vector3d(0.0, 1.0, 0.0));
    Eigen::Matrix3d A_shaft1 = caams::Ap(p_shaft1);
    Eigen::Matrix3d A_shaft2 = caams::Ap(p_shaft2);

    Eigen::Vector3d s1_p_shaft1uni(0.0, 0.0, SHAFT_LENGTH/2.0);
    Eigen::Vector3d s1_p_shaft1sph(0.0, 0.0,-SHAFT_LENGTH/2.0);
    Eigen::Vector3d s2_p_shaft2uni(0.0, 0.0,-SHAFT_LENGTH/2.0);
    Eigen::Vector3d s2_p_shaft2sph(0.0, 0.0, SHAFT_LENGTH/2.0);

    Eigen::Vector3d r_shaft1 = -A_shaft1*s1_p_shaft1uni;
    Eigen::Vector3d r_shaft2 = -A_shaft2*s2_p_shaft2uni;

    Eigen::Vector3d r_dot_shaft1 = Eigen::Vector3d::Zero();
    Eigen::Vector3d r_dot_shaft2 = Eigen::Vector3d::Zero();
    Eigen::Vector4d p_dot_shaft1 = Eigen::Vector4d::Zero();
    Eigen::Vector4d p_dot_shaft2 = Eigen::Vector4d::Zero();

    datumBody = new DatumBody();

    shaft1Body = new CylinderZaxis(
                r_shaft1,
                p_shaft1,
                r_dot_shaft1,
                p_dot_shaft1,
                SHAFT_MASS,
                SHAFT_RADIUS,
                SHAFT_LENGTH);

    shaft2Body = new CylinderZaxis(
                r_shaft2,
                p_shaft2,
                r_dot_shaft2,
                p_dot_shaft2,
                SHAFT_MASS,
                SHAFT_RADIUS,
                SHAFT_LENGTH);

    shaft1RevoluteJoint = new RevoluteJoint(
                datumBody, shaft1Body,
                r_shaft1 + A_shaft1*s1_p_shaft1sph,
                s1_p_shaft1sph,
                Eigen::Vector3d(0.0, 0.0, 1.0),
                Eigen::Vector3d(0.0, 0.0, 1.0));

    shaft2ParallelJoint = new Parallel1_2(
                datumBody, shaft2Body,
                A_shaft2.col(2),
                Eigen::Vector3d(0.0, 0.0, 1.0));

    universalJoint = new UniversalJoint(
                shaft1Body, shaft2Body,
                s1_p_shaft1uni, s2_p_shaft2uni,
                Eigen::Vector3d(0.0, 1.0, 0.0),
                Eigen::Vector3d(1.0, 0.0, 0.0));

    shaft1LocalTorque = new BodyLocalTorque(shaft1Body,
                                        Eigen::Vector3d(0.0, 0.0, TORQUE_MAX));

    Eigen::Matrix3d k_t = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d k_r;
    k_r << K_R, 0.0, 0.0,
           0.0, K_R, 0.0,
           0.0, 0.0, K_R;
    shaft2Friction = new BodyDamping(shaft2Body, k_t, k_r);

    mbsystem.AddBody(shaft1Body);
    mbsystem.AddBody(shaft2Body);

    mbsystem.AddConstraint(shaft1RevoluteJoint);
    mbsystem.AddConstraint(shaft2ParallelJoint);
    mbsystem.AddConstraint(universalJoint);

    mbsystem.AddForce(shaft1LocalTorque);
    mbsystem.AddForce(shaft2Friction);

    mbsystem.InitializeSolver();
}

double camera_longitude = 0.0;
glm::dmat4 camera_rotation=glm::rotate(glm::dmat4(1.0), CAMERA_LATITUDE, glm::dvec3(1.0, 0.0, 0.0));
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

    glm::dmat4 pan = glm::rotate(glm::dmat4(1.0), camera_longitude, glm::dvec3(0.0, 1.0, 0.0));
    glm::dmat4 camera_to_world = pan * camera_rotation * camera_translation;
    glm::dmat4 world_to_camera = glm::inverse(camera_to_world);

    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixd(glm::value_ptr(world_to_camera));

    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glColor3f(1.0f,1.0f,1.0f);
    mbsystem.Draw();

    glutSwapBuffers();

    mbsystem.Integrate(1/60.0);
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

bool mousing = false;
int x_last;
int y_last;

void mouseFunc(int button, int state, int x, int y)
{
    if(button==GLUT_LEFT_BUTTON){
        if(state==GLUT_DOWN){
            mousing = true;
            x_last = x;
            y_last = y;
        }else{
            mousing = false;
        }
    }
}

void motionFunc(int x, int y)
{
    if(mousing){
        int delta_x = x - x_last;
        int delta_y = y - y_last;
        camera_longitude -= delta_x*0.25/180.0*M_PI;
        x_last = x;
        y_last = y;
    }
}

int main(int argc, char **argv){
    glutInit(&argc,argv);
    init_system();
    glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );
    glutInitWindowSize (1280, 720);
    glutInitWindowPosition (100, 100);
    glutCreateWindow ("Perpetual Tennis Racket");
    glutDisplayFunc(display);
    glutTimerFunc( TIMER_INTERVAL, timerFunc, 0 );
    glutKeyboardFunc(keyboardFunc);
    glutMouseFunc(mouseFunc);
    glutMotionFunc(motionFunc);
    glClearColor(0.0,0.0,0.0,0.0);
    glDrawBuffer(GL_BACK);
    glColor3d(1.0,1.0,1.0);
    glutMainLoop();
    return 0;
}

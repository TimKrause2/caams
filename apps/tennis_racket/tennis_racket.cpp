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

#define CUBOID_DX 1.0
#define CUBOID_DY 3.0
#define CUBOID_DZ 0.2
#define CUBOID_MASS 1.0
#define OMEGA_X  5*2*M_PI
#define OMEGA_Y  0.2*2*M_PI
#define OMEGA_Z  0.0
#define CAMERA_LATITUDE (-15.0/180.0*M_PI)

Body *cuboidBody;
DatumBody datumBody;
System mbsystem;

void init_system(void)
{
    Eigen::Vector4d p_cuboid(1.0, 0.0, 0.0, 0.0);
    Eigen::Vector3d r_cuboid(0.0, 0.0, 0.0);
    Eigen::Vector3d omega_p(OMEGA_X, OMEGA_Y, OMEGA_Z);
    Eigen::Vector3d r_dot_cuboid(0.0, 0.0, 0.0);
    Eigen::Vector4d p_dot_cuboid(0.5*caams::L(p_cuboid).transpose()*omega_p);

    cuboidBody = new Cuboid(
                r_cuboid,
                p_cuboid,
                r_dot_cuboid,
                p_dot_cuboid,
                CUBOID_MASS,
                CUBOID_DX,
                CUBOID_DY,
                CUBOID_DZ);

    mbsystem.AddBody(cuboidBody);
    mbsystem.InitializeSolver();
}

double camera_longitude = 0.0;
glm::dmat4 camera_rotation=glm::rotate(glm::dmat4(1.0), CAMERA_LATITUDE, glm::dvec3(1.0, 0.0, 0.0));
glm::dmat4 camera_translation(glm::translate(glm::dvec3(0.0,0.0,3.0)));

#define NEAR_PLANE 0.1
#define FAR_PLANE 30.0
#define FOV (M_PI*90.0/180.0)

void draw_angular_velocity(void)
{
    Eigen::Vector3d omega = caams::omega_p_dot(
                cuboidBody->p, cuboidBody->p_dot);
    Eigen::Vector3d omega_p = caams::omega_p_p_dot(
                cuboidBody->p, cuboidBody->p_dot);
    Eigen::Vector3d L_p = cuboidBody->J_p*omega_p;
    Eigen::Matrix3d A = caams::Ap(cuboidBody->p);
    Eigen::Vector3d L = A*L_p;

    L.normalize();
    glBegin(GL_LINES);
    glColor3d(0.0, 0.5, 1.0);
    glVertex3d(0.0, 0.0, 0.0);
    glVertex3dv(L.data());
    glEnd();

    omega.normalize();
    glBegin(GL_LINES);
    glColor3d(1.0, 0.0, 0.0);
    glVertex3d(0.0, 0.0, 0.0);
    glVertex3dv(omega.data());
    glEnd();
}



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

    datumBody.Draw();

    draw_angular_velocity();

    glutSwapBuffers();

    mbsystem.Integrate(1/600.0);
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

#include "primitives.h"
#include <GL/gl.h>
#include <math.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>


void WireCylinder(double radius, double height, int slices, int stacks)
{
	glm::dvec3 vertices[stacks+1][slices];

	for(int sl=0;sl<slices;sl++){
		double theta = 2.0*M_PI*sl/slices;
		vertices[0][sl] = glm::dvec3(radius*cos(theta),radius*sin(theta),0.0);
	}

	for(int st=1;st<(stacks+1);st++){
		glm::dvec3 z(0.0, 0.0, height*st/stacks);
		for(int sl=0;sl<slices;sl++){
			vertices[st][sl] = vertices[0][sl] + z;
		}
	}

	int sl0 = 0;
	int sl1 = 1;
	glBegin(GL_LINES);

	for(;sl0<slices;sl0++,sl1++){
		if(sl1==slices)sl1=0;
		int st;
		for(st=0;st<stacks;st++){
			glVertex3dv(glm::value_ptr(vertices[st][sl0]));
			glVertex3dv(glm::value_ptr(vertices[st][sl1]));
			glVertex3dv(glm::value_ptr(vertices[st][sl0]));
			glVertex3dv(glm::value_ptr(vertices[st+1][sl0]));
		}
		glVertex3dv(glm::value_ptr(vertices[st][sl0]));
		glVertex3dv(glm::value_ptr(vertices[st][sl1]));
	}

	glEnd();

}

void WireCube(double scale)
{
	scale *= 0.5;

	glBegin(GL_LINES);

	glVertex3d(scale,scale,scale);
	glVertex3d(scale,scale,-scale);

	glVertex3d(scale,scale,-scale);
	glVertex3d(scale,-scale,-scale);

	glVertex3d(scale,-scale,-scale);
	glVertex3d(scale,-scale,scale);

	glVertex3d(scale,-scale,scale);
	glVertex3d(scale,scale,scale);

	glVertex3d(-scale,scale,scale);
	glVertex3d(-scale,scale,-scale);

	glVertex3d(-scale,scale,-scale);
	glVertex3d(-scale,-scale,-scale);

	glVertex3d(-scale,-scale,-scale);
	glVertex3d(-scale,-scale,scale);

	glVertex3d(-scale,-scale,scale);
	glVertex3d(-scale,scale,scale);

	glVertex3d(scale,scale,scale);
	glVertex3d(-scale,scale,scale);

	glVertex3d(scale,scale,-scale);
	glVertex3d(-scale,scale,-scale);

	glVertex3d(scale,-scale,-scale);
	glVertex3d(-scale,-scale,-scale);

	glVertex3d(scale,-scale,scale);
	glVertex3d(-scale,-scale,scale);

	glEnd();
}

void WireSphere(double radius, int slices, int stacks)
{
	if(stacks<2)return;

	glm::dvec3 pole_p(0.0, 0.0, radius);
	glm::dvec3 pole_m(0.0, 0.0, -radius);

	glm::dvec3 vertices[stacks-1][slices];

	for(int st=0;st<(stacks-1);st++){
		double phi = M_PI*(st+1)/stacks;
		for(int sl=0;sl<slices;sl++){
			double theta = 2.0*M_PI*sl/slices;
			vertices[st][sl] = glm::dvec3(radius*sin(phi)*cos(theta),
										  radius*sin(phi)*sin(theta),
										  radius*cos(phi));
		}
	}

	glBegin(GL_LINES);

	for(int sl0=0,sl1=1;sl0<slices;sl0++,sl1++){
		if(sl1==slices)sl1=0;
		glVertex3dv(glm::value_ptr(pole_p));
		glVertex3dv(glm::value_ptr(vertices[0][sl0]));

		glVertex3dv(glm::value_ptr(pole_m));
		glVertex3dv(glm::value_ptr(vertices[stacks-2][sl0]));

		glVertex3dv(glm::value_ptr(vertices[0][sl0]));
		glVertex3dv(glm::value_ptr(vertices[0][sl1]));

		for(int st=1;st<stacks-1;st++){
			glVertex3dv(glm::value_ptr(vertices[st-1][sl0]));
			glVertex3dv(glm::value_ptr(vertices[st][sl0]));

			glVertex3dv(glm::value_ptr(vertices[st][sl0]));
			glVertex3dv(glm::value_ptr(vertices[st][sl1]));
		}
	}
	glEnd();
}



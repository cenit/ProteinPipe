#include <vector>
#include <cmath>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#define WINDING_NUMBER 30
#ifndef M_PI
#define M_PI 3.141519
#endif
#ifndef M_PI_2
#define M_PI_2 1.570759
#endif // !M_PI_2

class SolidSphere
{
protected:
  std::vector<GLfloat> vertices;
  std::vector<GLfloat> normals;
  std::vector<GLfloat> texcoords;
  std::vector<GLushort> indices;

public:
  SolidSphere() {};
  SolidSphere(const float &, const int &, const int &);
  void draw(const GLfloat &, 
            const GLfloat &, 
            const GLfloat &, 
            const GLfloat &, 
            const GLfloat &, 
            const GLfloat &
            );
};

SolidSphere sphere(.01f, 12, 24);

SolidSphere::SolidSphere(const float &radius, const int &rings, const int &sectors)
{
  float const R = 1.f / (float)(rings - 1);
  float const S = 1.f / (float)(sectors - 1);
  int r, s;

  vertices.resize(rings * sectors * 3);
  normals.resize(rings * sectors * 3);
  texcoords.resize(rings * sectors * 2);
  std::vector<GLfloat>::iterator v = vertices.begin();
  std::vector<GLfloat>::iterator n = normals.begin();
  std::vector<GLfloat>::iterator t = texcoords.begin();
  for (r = 0; r < rings; ++r) for (s = 0; s < sectors; ++s)
  {
    float const y = (float)std::sin(-M_PI_2 + M_PI * r * R);
    float const x = (float)std::cos(2 * M_PI * s * S) * (float)std::sin(M_PI * r * R);
    float const z = (float)std::sin(2 * M_PI * s * S) * (float)std::sin(M_PI * r * R);

    *t++ = s*S;
    *t++ = r*R;

    *v++ = x * radius;
    *v++ = y * radius;
    *v++ = z * radius;

    *n++ = x;
    *n++ = y;
    *n++ = z;
  }

  indices.resize(rings * sectors * 4);
  std::vector<GLushort>::iterator i = indices.begin();
  for (r = 0; r < rings - 1; ++r) for (s = 0; s < sectors - 1; ++s)
  {
    *i++ = r * sectors + s;
    *i++ = r * sectors + (s + 1);
    *i++ = (r + 1) * sectors + (s + 1);
    *i++ = (r + 1) * sectors + s;
  }
}

void SolidSphere::draw( const GLfloat &x, 
            const GLfloat &y, 
            const GLfloat &z, 
            const GLfloat &R, 
            const GLfloat &G, 
            const GLfloat &B
            )
{
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glColor3f(R, G, B);
  glTranslatef(x, y, z);

  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);

  glVertexPointer(3, GL_FLOAT, 0, &vertices[0]);
  glNormalPointer(GL_FLOAT, 0, &normals[0]);
  glTexCoordPointer(2, GL_FLOAT, 0, &texcoords[0]);
  glDrawElements(GL_QUADS, indices.size(), GL_UNSIGNED_SHORT, &indices[0]);
  glPopMatrix();
}



template<int N> void link(const GLfloat &radius, 
						  const GLfloat &x1, const GLfloat &y1, const GLfloat &z1, // first coordinates
						  const GLfloat &x2, const GLfloat &y2, const GLfloat &z2, // second coordinates
						  const GLfloat &R, const GLfloat &G, const GLfloat &B // colors
						  );

template<> void link<0>(  const GLfloat &radius, 
				          const GLfloat &x1, const GLfloat &y1, const GLfloat &z1, 
				          const GLfloat &x2, const GLfloat &y2, const GLfloat &z2, 
				          const GLfloat &R,  const GLfloat &G,  const GLfloat &B
				          ) // link as simple lines
{
	glBegin(GL_LINES);
	    glNormal3f(1.f, 1.f, 1.f);
	    glColor3f(R, G, B);
	    glVertex3f(x1, y1, z1);
	    glVertex3f(x2, y2, z2);
  	glEnd();
  	glPopMatrix();
}

template<> void link<1>(  const GLfloat &radius, 
				          const GLfloat &x1, const GLfloat &y1, const GLfloat &z1, 
				          const GLfloat &x2, const GLfloat &y2, const GLfloat &z2, 
				          const GLfloat &R,  const GLfloat &G,  const GLfloat &B
				          ) // link as springs
{
	// to fix with glRotatef
	glPushMatrix();
	//glTranslatef(x1, y1, z1); // Translate to the object's position.
	glBegin(GL_LINE_STRIP);
		glNormal3f(1.f, 1.f, 1.f);
		glColor3f(R, G, B);
		glVertex3f(x1, y1, z1);
		
		GLfloat ax = y2 - y1, ay = x1 - x2, az = 0.f;
		GLfloat norm = std::sqrt(ax*ax + ay*ay + az*az);
		if (norm != 0) { ax /= norm; ay /= norm; az /= norm; }
		GLfloat distance = std::sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) + (z2 - z1)*(z2 - z1)), ctheta = (z2 - z1) / distance, stheta = std::sqrt(1 - ctheta*ctheta);
		GLfloat r00 = ax*ax*(1 - ctheta) + ctheta,    r01 = ax*ay*(1 - ctheta),     r02 = -ay*stheta,
		    r10 = ax*ay*(1 - ctheta),       r11 = ay*ay*(1 - ctheta) + ctheta,  r12 = ax*stheta,
		    r20 = ay*stheta,            r21 = -ax*stheta,         r22 = ctheta;
		GLfloat x, y, z = z1, angle = 0.f, tmp_z = 0.f;
		while (angle < WINDING_NUMBER*2*M_PI)
		{
		  GLfloat tmp_x = radius * std::cos(angle), tmp_y = radius * std::sin(angle);
		  x = r00*tmp_x + r01*tmp_y + r02*tmp_z + x1;
		  y = r10*tmp_x + r11*tmp_y + r12*tmp_z + y1;
		  glVertex3f(x, y, z);
		  angle += .1f;
		  tmp_z = angle*distance / GLfloat(WINDING_NUMBER * 2 * M_PI);
		  z = r20*tmp_x + r21*tmp_y + r22*tmp_z + z1;
		}
		// Done drawing points
  	glEnd();
  	glPopMatrix();
}

template<> void link<2>(  const GLfloat &radius, 
				          const GLfloat &x1, const GLfloat &y1, const GLfloat &z1, 
				          const GLfloat &x2, const GLfloat &y2, const GLfloat &z2, 
				          const GLfloat &R,  const GLfloat &G,  const GLfloat &B
				          ) // link as cylinder
{
  glPushMatrix();
  glTranslatef(x1, y1, z1); // Translate to the object's position.
  GLfloat rx = x2 - x1 , ry = y2 - y1, rz = z2 - z1;
  GLfloat distance = std::sqrt(rx*rx + ry*ry + rz*rz);
  GLfloat ax = ry / distance, ay = -rx / distance;
  GLfloat norm = std::sqrt(ax*ax + ay*ay);
  GLfloat theta = (GLfloat)(- std::acos(rz / distance) / M_PI * 180.f);

  glRotatef(theta, ax / norm, ay / norm, 0.f); // Rotate the object.
  glColor3f(R, G, B);
  glNormal3f(1.f, 1.f, 1.f);
  glBegin(GL_POLYGON);
	GLUquadricObj *obj = gluNewQuadric();
	gluCylinder(obj, radius, radius, distance, 10, 10); // signatures function (obj, low_ray, up_ray, height, nquad, nquad)
  glEnd();

  glPopMatrix();
  return;
}

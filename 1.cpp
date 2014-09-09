#include <GL/glut.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h> 
#include <utility>

using namespace std;

double spx,spy,spz;
int viewx,viewy;

struct vector3d{
    double x,y,z;
};

class Sphere{
  vector3d c;
  double r;

};

class Plane{
  vector3d x,y,z;

};

class Object{
  int type;
  Sphere s;
  Plane p;
  double ka[3],kd[3],ks[3];
};




vector3d cross_prod(vector3d u,vector3d v)
{
    vector3d a;
    a.x=u.y*v.z-u.z*v.y;
    a.y=u.z*v.x-u.x*v.z;
    a.z=u.x*v.y-v.x*u.y;
    return a;
}

vector3d norm(vector3d c)
{
    vector3d p;
    double a=sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
    p.x=c.x/a;p.y=c.y/a;p.z=c.z/a;
    return p;
}

vector3d scalarmulti(double x,vector3d y)
{
    vector3d o;
    o.x=x*y.x;
    o.y=x*y.y;
    o.z=x*y.z;
    return o;
}

vector3d subvec(vector3d a,vector3d b)
{
    vector3d g;
    g.x=a.x-b.x;
    g.y=a.y-b.y;
    g.z=a.z-b.z;
    return g;
}
vector3d addvec(vector3d a,vector3d b)
{
    vector3d g;
    g.x=a.x+b.x;
    g.y=a.y+b.y;
    g.z=a.z+b.z;
    return g;
}
double dot(vector3d a,vector3d b)
{
  return a.x*b.x+a.y*b.y+a.z*b.z;
}


void specialKeys( int key, int x, int y ) {
 
  //  Right arrow - increase rotation by 5 degree
  // if (key == GLUT_KEY_RIGHT)
  //   rotate_y += 5;
 
  // //  Left arrow - decrease rotation by 5 degree
  // else if (key == GLUT_KEY_LEFT)
  //   rotate_y -= 5;
 
  // else if (key == GLUT_KEY_UP)
  //   rotate_x += 5;
 
  // else if (key == GLUT_KEY_DOWN)
  //   rotate_x -= 5;
  // else if(key==GLUT_KEY_F1)
  //   view=0;
  // else if(key==GLUT_KEY_F2)
  //   view=1;
  // else if(key==GLUT_KEY_F3)
  //   view=2;
  // else if(key==GLUT_KEY_F4)
  //   view=3;
  // else if(key==GLUT_KEY_F7)
  //   view=4;

 
  //  Request display update
  glutPostRedisplay();
 
}




void display()
{

     glMatrixMode(GL_MODELVIEW);
    // clear the drawing buffer.
    glClear(GL_COLOR_BUFFER_BIT);
    glFlush();        
   glutSwapBuffers();

}


void reshape(int x, int y)
{
    if (y == 0 || x == 0) return;   
    glMatrixMode(GL_PROJECTION);  
    glLoadIdentity(); 
    //gluPerspective(39.0,(GLdouble)x/(GLdouble)y,0.6,21.0);
    glOrtho(-10.0, 10.0, -10.0, 10.0, -100.0, 100.0);
    glMatrixMode(GL_MODELVIEW);
    glViewport(0,0,x,y);  //Use the whole window for rendering
} 



int main (int argc, char **argv)
{
    
    glutInit(&argc, argv); 
    viewx=350;viewy=350;
    glutInitWindowSize(viewx,viewy);
    glutCreateWindow("Solid Sphere");
    spx=0.0;spy=0.0;spz=0.0;
    glutSpecialFunc(specialKeys);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMainLoop();
    return 0;
}

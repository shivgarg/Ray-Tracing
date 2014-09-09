#include <GL/glut.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <stdio.h>
#include "math.h"
#include <stdlib.h> 
#include <utility>
#include "vec.cpp"
using namespace std;

double spx,spy,spz;
int viewx,viewy;



class Ray{
 public: 
  vector3d org,dir;
  Ray(vector3d a,vector3d b)
  {
    org=a;dir=b;
  }

};




class Sphere{
 public: 
  vector3d c;
  double r;
  Sphere(vector3d a,double rad)
  {
    c=a;r=rad;
  }

};

class Plane{
public:
  vector3d x,y,z;
  Plane(){};
  Plane(vector3d *a,vector3d *b, vector3d *c){
    x=*a;
    y=*b;
    z=*c;
  }


};

class Object{
  int type;
  Sphere s;
  Plane p;
  double ka[3],kd[3],ks[3];
  double exp;
};


 class lightsource{
  double intensity[3];
  vector3d direction;
  
 };



bool intersectPlane(Plane p, vector3d origin,vector3d direct , float &t)
{
    // assuming vectors are all normalized

    vector3d normal =  norm(cross_prod(subvec(p.x,p.y),subvec(p.x,p.z)));
    float denom = dot(normal, direct);
    if (denom > 1e-6) {
        vector3d a = subvec(p.x , origin);
        t = dot(a, normal) / denom; 
        return (t >= 0);
    }
    return false;
}


bool intersectsphere(Sphere * c,double * t,Ray * r)
{
      vector3d l=cross_prod(subvec(r->org,c->c),r->dir);
      double dist=mod(l);
      if(dist>c->r)
      {
        return false;
      }
      else
      {        
        double b=2.0*(r->dir.x*(r->org.x-c->c.x)+r->dir.y*(r->org.y-c->c.y)+r->dir.z*(r->org.z-c->c.z));
        double d=(r->org.x-c->c.x)*(r->org.x-c->c.x)+(r->org.y-c->c.y)*(r->org.y-c->c.y)+(r->org.y-c->c.z)*(r->org.z-c->c.z)-(c->r)*(c->r);
        double t0=(-b-sqrt(b*b-4*d))/2.0;
        double t1=(-b+sqrt(b*b-4*d))/2.0;
        if(t0<t1)
          *t=t0;
        else
          *t=t1;
        return true;

      }


}

struct rgbf {float r; float g; float b;};
//WBL 9 May 2007 Based on
//http://www.codeguru.com/cpp/w-d/dislog/commondialogs/article.php/c1861/
//Common.h
void toRGBf(const float h, const float s, const float v,
      rgbf* rgb)
{
  /*
RGBType rgb;
  if(!h  && !s)
  {
    rgb.r = rgb.g = rgb.b = v;
  }
  */
  //rgbf* rgb = (rgbf*) out;
double min,max,delta,hue;
  
  max = v;
  delta = max * s;
  min = max - delta;

  hue = h;
  if(h > 300 || h <= 60)
  {
    rgb->r = max;
    if(h > 300)
    {
      rgb->g = min;
      hue = (hue - 360.0)/60.0;
      rgb->b = ((hue * delta - min) * -1);
    }
    else
    {
      rgb->b = min;
      hue = hue / 60.0;
      rgb->g = (hue * delta + min);
    }
  }
  else
  if(h > 60 && h < 180)
  {
    rgb->g = max;
    if(h < 120)
    {
      rgb->b = min;
      hue = (hue/60.0 - 2.0 ) * delta;
      rgb->r = min - hue;
    }
    else
    {
      rgb->r = min;
      hue = (hue/60 - 2.0) * delta;
      rgb->b = (min + hue);
    }
  }
  else
  {
    rgb->b = max;
    if(h < 240)
    {
      rgb->r = min;
      hue = (hue/60.0 - 4.0 ) * delta;
      rgb->g = (min - hue);
    }
    else
    {
      rgb->g = min;
      hue = (hue/60 - 4.0) * delta;
      rgb->r = (min + hue);
    }
  }
}


//Convert a wide range of data values into nice colours 
void colour(const float data, float* out) {
  //convert data to angle
  const float a = atan2(data,1)/(2*atan2(1,1)); // -1 .. +1
  const float angle = (1+a)*180; //red=0 at -1,+1

  const float saturation = 1;

  const float h = (data<-1||data>1)? 1 : fabs(data);

  toRGBf(angle,saturation,h,(rgbf*)out);
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
    float* pixels = new float[viewx*viewy*3];
    for(int i=0;i<viewx*viewy;i++) {
    colour(10.0-((i*20.0)/(viewx*viewy)),&pixels[i*3]);
  } 

  //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  //http://msdn2.microsoft.com/en-us/library/ms537062.aspx
  //glDrawPixels writes a block of pixels to the framebuffer.

  glDrawPixels(viewx,viewy,GL_RGB,GL_FLOAT,pixels);

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
    viewx=512;viewy=512;
    glutInitWindowSize(viewx,viewy);
    glutCreateWindow("Solid Sphere");
    double t=0;
    Plane p(new vector3d(1.0,0.0,0.0),new vector3d(0.0,1.0,0.0),new vector3d(0.0,0.0,0.0));
    vector3d a(0.1,0.1,-0.2);
    vector3d b(0.0,0.0,1.0);
 //   p.x=a;
    //p.y=a;
    //p.z=a;
   // p.y(1.0,0.0,0.0);
    //p.z = new vector3d(1,1,1);
   // intersectPlane(p,a,b,t);
    Sphere sd(vector3d(0,0,0),5);
    Ray rt(vector3d(0,-10,0),vector3d(1,0,0));
    bool y=intersectsphere(&sd,&t,&rt);
    if(y)
      cout <<"here"<<t <<endl;
    else
      cout << "didnt intersect"<< endl;

    spx=0.0;spy=0.0;spz=0.0;
    glutSpecialFunc(specialKeys);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMainLoop();



    return 0;
}

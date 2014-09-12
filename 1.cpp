#include <GL/glut.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <algorithm>
#include "math.h"
#include <stdlib.h> 
#include <utility>
#include <climits>
#include "vec.cpp"
using namespace std;

//double spx,spy,spz;
int viewx,viewy;
vector3d vcs,u,v,n,eye;
double disteye;
int viewl,viewb;
//int currentpixelintensity[3]={0,0,0};
float ambientintensity[3]={0.1,0.1,0.1};



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
  Sphere(){};
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
    x=*a;y=*b;z=*c;
  }
};

class Object{
 public: 
  int type;
  Sphere s;
  Plane p;
  double ka[3],kd[3],ks[3];
  double exp;
  Object(int a, Sphere b){
    type = a;
    s = b;
  }
  Object(int a, Plane b){
    type = a;
    p = b;
  }

};


 class lightsource{
 public:
  double intensity[3];
  vector3d direction;
  lightsource(){};
  lightsource(double a,double b,double c,vector3d d){intensity[0]=a;intensity[1]=b;intensity[2]=c;direction=d;};
 };

vector<Object> obj;
vector<lightsource> lightsources;

double area(vector3d a, vector3d b, vector3d c)
{
   return mod(cross_prod(subvec(a,b),subvec(a,c)));
}

void ambient( float*f,Object o){
  for(int i=0;i<3;i++){
    ////cout<<"before "<<*(f+i)<<endl;
    *(f+i) += o.ka[i]*ambientintensity[i];
    ////cout<<"after "<<o.ka[i]*ambientintensity[i]<<endl;


//  currentpixelintensity[1]+= o.ka[1]*ambientintensity[1];
//  currentpixelintensity[2]+= o.ka[2]*ambientintensity[2];
  }

  ////cout<<"ambient "<<*(f)<<endl;
}

void diffuse(float*f,lightsource l, vector3d normal, Object o){
  for(int i=0;i<3;i++){
    double temp = o.kd[i]*(l.intensity[i])*dot(l.direction,normal);
    // if(temp<0){
    //   temp=0;
    // }
    *(f+i)+= max(-1*temp,0.0);
  }
  // currentpixelintensity[1]+= o.kd[0]*(l.intensity[0])*dot(l.direction,normal);
  // currentpixelintensity[2]+= o.kd[0]*(l.intensity[0])*dot(l.direction,normal);
  ////cout<<"ambient "<<*(f)<<endl;

}

void specular(float*f,lightsource l, vector3d normal, vector3d v,Object o){
  vector3d r = addvec(l.direction,scalarmulti(-2*dot(l.direction,normal),normal));
    //vector3d v = ;
  for(int i=0;i<3;i++){
    double temp = o.ks[i]*(l.intensity[i])*(pow(dot(r,v),o.exp));
    //if(o.exp)
    *(f+i)+= max(temp,0.0);
  }
  // currentpixelintensity[1]+= o.kd[0]*(l.intensity[0])*dot(l.direction,normal);
  // currentpixelintensity[2]+= o.kd[0]*(l.intensity[0])*dot(l.direction,normal);

}


bool intersectPlane(Plane p, vector3d origin,vector3d direct , double &t,int px,int py)
{
    // assuming vectors are all normalized

    vector3d normal =  norm(cross_prod(subvec(p.x,p.y),subvec(p.x,p.z)));
    float denom = dot(normal, direct);
    if (denom > 1e-6) {
        vector3d a = subvec(p.x , origin);
        t = dot(a, normal) / denom;
        vector3d point;
        point  = addvec(origin, scalarmulti(t,direct));
        if(t>=0){
          ////cout<<area(point,p.x,p.y)<<" "<<area(point,p.z,p.y)<<" "<<area(point,p.x,p.z)<<" "<<area(p.z,p.x,p.y)<<endl;
          return ((area(point,p.x,p.y)+area(point,p.z,p.y)+area(point,p.x,p.z)) == area(p.z,p.x,p.y));
        }
    }
    return false;
}


bool intersectsphere(Sphere * c,double *t,Ray * r,int px,int py)
{
      vector3d l=cross_prod(subvec(r->org,c->c),r->dir);

      double dist=mod(l);
      //cout << dist << "distance"<< endl;  
      if(dist>c->r)
      {
        return false;
      }
      else
      {        
        double b=2.0*(r->dir.x*(r->org.x-c->c.x)+r->dir.y*(r->org.y-c->c.y)+r->dir.z*(r->org.z-c->c.z));
        double d=(r->org.x-c->c.x)*(r->org.x-c->c.x)+(r->org.y-c->c.y)*(r->org.y-c->c.y)+(r->org.z-c->c.z)*(r->org.z-c->c.z)-(c->r)*(c->r);
        double t0=(-b-sqrt(b*b-4.0*d))/2.0;
        double t1=(-b+sqrt(b*b-4.0*d))/2.0;
        //cout << "b "<< b <<  " d "<< b << endl;
//      
        //cout << "to "<< t0 <<  " t1 "<< t1 << endl;
//        if(t0<t1)

          if((min(t0,t1))>0){
           *t=min(t0,t1);
            return true;

          }
          else{
            return false;
          }
      }


}






void display()
{
    glMatrixMode(GL_MODELVIEW);
    // clear the drawing buffer.
    glClear(GL_COLOR_BUFFER_BIT);
    float* pixels = new float[viewx*viewy*3];
    for(int i=0;i<viewx*viewy;i++) {
    //colour(10.0-((i*20.0)/(viewx*viewy)),&pixels[i*3]);
  } 

  
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


vector3d viewporttovcs(int x,int y)
{
  x=x-viewx/2.0;y=viewy/2.0-y;
  float x1  =  (float)(x)*((float)viewl)/((float)viewx);
  float y1 = (float)(y)*((float)viewb)/((float)viewy);
  vector3d f=addvec(scalarmulti(x1,u),scalarmulti(y1,v));
  f=addvec(f,vcs);
  return f;
}

// void rayTracerRecurse(Ray rt,)
// {
//    glMatrixMode(GL_MODELVIEW);
//   glClear(GL_COLOR_BUFFER_BIT);
//   float* pixels = new float[viewx*viewy*3];
//   for(int i=0;i<viewx*viewy;i++)
//   {
//     vector3d f=viewporttovcs(i%viewx,i/viewy);
//     Ray rt(eye,norm(subvec(f,eye)));
//     double t = INT_MAX;
//     int h=obj.size();
//     double k=t;
//     vector3d normal;
//     int intersectid;
//     for(int j=0;j<h;j++)
//     {
//       Object a=obj[j];
//       if(a.type==1)
//         intersectsphere(&(a.s),&k,&rt,i%viewx,i/viewy);
//       else
//         intersectPlane((a.p),rt.org,rt.dir,k,i%viewx,i/viewy);
//       if(k<t){
//         t=k;
//         intersectid=j;
//       } 
//     }
//     if(t!=INT_MAX){//shadow part pending
//       Object a = obj[intersectid];
//       if(a.type==1){
//         normal = norm(subvec(addvec(rt.org,scalarmulti(t,rt.dir)),a.s.c));
//       }
//       else{
//         normal = norm(cross_prod(subvec(a.p.x,a.p.y),subvec(a.p.x,a.p.z)));
//       }
//       pixels[3*i] = 0;
//       pixels[3*i+1] = 0;
//       pixels[3*i+2] = 0;
//       ambient(&pixels[3*i],a);
//       for(int r=0;r<lightsources.size();r++){
//         diffuse(&pixels[3*i],lightsources[r],normal,a);
//       }
//       for(int r=0;r<lightsources.size();r++){
//         specular(&pixels[3*i],lightsources[r],normal,a);
//       }
//     }
//   }
//     glDrawPixels(viewx,viewy,GL_RGB,GL_FLOAT,pixels);
//     glFlush();        
//     glutSwapBuffers();
// }



void raytracer()
{
  glMatrixMode(GL_MODELVIEW);
  glClear(GL_COLOR_BUFFER_BIT);
  float* pixels = new float[viewx*viewy*3];
  for(int i=0;i<viewx*viewy;i++)
  {
    vector3d f=viewporttovcs(i%viewx,i/viewx);

    Ray rt(eye,norm(subvec(f,eye)));
     //cout << rt.org.x<< " "<< rt.org.y << " "<<rt.org.z << " primary ray origin "<< endl;
     //cout << rt.dir.x<< " "<< rt.dir.y << " "<<rt.dir.z << " primary ray direction "<< endl;
    double t = INT_MAX;
    int h=obj.size();
    double k=t;
    vector3d normal;
    int intersectid;
    for(int j=0;j<h;j++)
    {
      Object a=obj[j];bool p;
      if(a.type==1)
        {
          //cout << "F "<< f.x<<" "<< f.y<<" "<< f.z << "vcs "<<endl;
          //cout << rt.org.x<< " "<< rt.org.y << " "<<rt.org.z << " primary ray origin "<< endl;
          //cout << rt.dir.x<< " "<< rt.dir.y << " "<<rt.dir.z << " primary ray direction "<< endl;
          //cout <<a.s.c.x<<" "<<a.s.c.y<<" "<<a.s.c.z<<" sphere origin  "<<a.s.r<<endl;

          p=intersectsphere(&(a.s),&k,&rt,i%viewx,i/viewy);
          if(!p)
          {
            //cout << "Not intersecting"<< endl;

          }
          else
          {
            vector3d asd=addvec(rt.org,scalarmulti(k,rt.dir));
            //cout << "point intersection "<< asd.x<<" "<< asd.y<<" "<< asd.z<< endl;
          }

        }
      else
        p=intersectPlane((a.p),rt.org,rt.dir,k,i%viewx,i/viewy);
      if(k<t){
        t=k;
        intersectid=j;
      } 
    }
    if(t!=INT_MAX){//shadow part pending
      //cout << "T in dsdsdsds  " << t << endl;
      Object a = obj[intersectid];
      if(a.type==1){
        normal = norm(subvec(addvec(rt.org,scalarmulti(t,rt.dir)),a.s.c));
      }
      else{
        normal = norm(cross_prod(subvec(a.p.x,a.p.y),subvec(a.p.x,a.p.z)));
      }
      pixels[3*i] = 0;
      pixels[3*i+1] = 0;
      pixels[3*i+2] = 0;
      ambient(&pixels[3*i],a);

      for(int r=0;r<lightsources.size();r++)
      {
        bool l =false;
        vector3d src=addvec(rt.org,scalarmulti(t,rt.dir));
        //cout << src.x<< " "<< src.y << " "<<src.z << " ray origin "<< endl;
        Ray ty(src,scalarmulti(-1,lightsources[r].direction));
        for(int j=0;j<h;j++)
        {
          if(j==intersectid){
            continue;
          }
          Object a=obj[j];
          if(a.type==1)
            {
            //cout << "herererererer   "<<endl;
               l=intersectsphere(&(a.s),&k,&ty,i%viewx,i/viewy);
               //cout<< ty.org.x<<endl;
               if(l)
                { //cout << "J "<< j << "R "<< intersectid<< endl;
              break;}

            }
          else
          {
            l=intersectPlane((a.p),ty.org,ty.dir,k,i%viewx,i/viewy);
            if(l)
              break;
          }
        }
        if(!l)
          {
            diffuse(&pixels[3*i],lightsources[r],normal,a);specular(&pixels[3*i],lightsources[r],normal,scalarmulti(-1,rt.dir),a);
          }
         
      }
      // for(int r=0;r<lightsources.size();r++)
      // {
      //   bool l =false;
      //   vector3d src=addvec(rt.org,scalarmulti(t,rt.dir);
      //   Ray ty(src,subvec(vector3d(0,0,0),lightsources[r].direction));
      //   for(int j=0;j<h;j++)
      //   {
      //     Object a=obj[j];
      //     if(a.type==1)
      //       {
      //          l=intersectsphere(&(a.s),&k,&ty,i%viewx,i/viewy);
      //          if(l)
      //           break;

      //       }
      //     else
      //     {
      //       l=intersectPlane((a.p),ty.org,ty.dir,k,i%viewx,i/viewy);
      //       if(l)
      //         break;
      //     }
      //   }
      //   if(!l)
      //     specular(&pixels[3*i],lightsources[r],normal,a);
      // }
    }
  }

    glDrawPixels(viewx,viewy,GL_RGB,GL_FLOAT,pixels);
    glFlush();        
    glutSwapBuffers();
}




int main (int argc, char **argv)
{
    
    glutInit(&argc, argv); 
    viewx=500;viewy=500;viewl = 50;
    viewb = 50;
    //vector3d f=viewporttovcs(0,viewy/2);
    ////cout <<"sdsdsds   "<< f.x <<" "<< f.y << " "<<f.z <<endl;
    glutInitWindowSize(viewx,viewy);
    glutCreateWindow("Solid Sphere");
    u = vector3d(1,0,0);
    v = vector3d(0,1,0);
    vector3d vcs = vector3d(0,0,0);
    disteye = 12;
    
    n=cross_prod(u,v);
    n=norm(n);
    

    eye=subvec(vcs,scalarmulti(disteye,n));
    //cout <<"normal "<< eye.x << " "<<eye.y<<" "<< eye.z<< endl;
    Object o(1, Sphere(vector3d(0.0,0.0,0.0),6));
    o.ka[0]=1;
    o.ka[1]=0;
    o.ka[2]=0;
    o.kd[0]=0.3;
    o.kd[1]=0.5;
    o.kd[2]=1.0;
    o.ks[0]=1.0;
    o.ks[1]=1.0;
    o.ks[2]=1.0;
    o.exp=50.0;
   obj.push_back(o);
    o = Object(1, Sphere(vector3d(0.0,0.0,-7.0),1.5));
    o.ka[0]=1;
    o.ka[1]=1;
    o.ka[2]=1;
    o.kd[0]=1;
    o.kd[1]=1;
    o.kd[2]=1;
    o.ks[0]=1;
    o.ks[1]=1;
    o.ks[2]=1;
    o.exp=250.0;
    obj.push_back(o);
    lightsource ll = lightsource( 1,1,1 ,vector3d(-0.707,0,0.707));
    //lightsource ll = lightsource( 1,1,1 ,vector3d(0.577,0.577,0.577));
    //lightsource ll = lightsource( 1,1,1 ,vector3d(0,0,1.0));
    lightsources.push_back(ll);
    //obj.push_back();

    double t=0;
 //    Plane p(new vector3d(1.0,0.0,0.0),new vector3d(0.0,1.0,0.0),new vector3d(0.0,0.0,0.0));
 //    vector3d a(0.9,0.9,-0.2);
 //    vector3d b(0.0,0.0,1.0);
 // //   p.x=a;
 //    //p.y=a;
 //    //p.z=a;
 //   // p.y(1.0,0.0,0.0);
 //    //p.z = new vector3d(1,1,1);
 //     bool z = intersectPlane(p,a,b,t,0,0);
     ////cout<< "here"<<t<<endl;
    // Sphere sd(vector3d(0,0,0),5);
    // Ray rt(vector3d(0,-10,0),vector3d(1,0,0));
    // bool y=intersectsphere(&sd,&t,&rt,0,0);
    // if(z)
    //   //cout <<"here"<<t <<endl;
    // else
    //   //cout << "didnt intersect"<< endl;



    //spx=0.0;spy=0.0;spz=0.0;
   
    glutDisplayFunc(raytracer);
    glutReshapeFunc(reshape);
    glutMainLoop();



    return 0;
}

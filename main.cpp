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

#define MAX_DEPTH 3

using namespace std;

int viewx,viewy;
vector3d vcs,u,v,n,eye;
double disteye;
int viewl,viewb;
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
                double kreflec,krefrac,n;
                double exp;
                Object(){};
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
                *(f+i) += o.ka[i]*ambientintensity[i];
        }
}

vector3d ambientrecur( Object o){
        vector3d ans(0,0,0);
        ans.x = o.ka[0]*ambientintensity[0];
        ans.y = o.ka[1]*ambientintensity[1];
        ans.z = o.ka[2]*ambientintensity[2];
        return ans;
}



void diffuse(float*f,lightsource l, vector3d normal, Object o){
        for(int i=0;i<3;i++){
                double temp = o.kd[i]*(l.intensity[i])*dot(l.direction,normal);
                *(f+i)+= max(-1*temp,0.0);
        }
}


vector3d diffuserecur(lightsource l, vector3d normal, Object o){
        vector3d ans(0,0,0);
        double temp = o.kd[0]*(l.intensity[0])*dot(l.direction,normal);
        ans.x = max(-1*temp,0.0);
        temp = o.kd[1]*(l.intensity[1])*dot(l.direction,normal);
        ans.y = max(-1*temp,0.0);
        temp = o.kd[2]*(l.intensity[2])*dot(l.direction,normal);
        ans.z = max(-1*temp,0.0);
        return ans;
}


vector3d specularrecur(lightsource l, vector3d normal, vector3d v,Object o){
        vector3d ans(0,0,0);
        vector3d r = addvec(l.direction,scalarmulti(-2*dot(l.direction,normal),normal));
        double temp = o.ks[0]*(l.intensity[0])*(pow(dot(r,v),o.exp));
        if(dot(r,v)>0){
                ans.x+= max(temp,0.0);
        }
        temp = o.ks[1]*(l.intensity[1])*(pow(dot(r,v),o.exp));
        if(dot(r,v)>0){
                ans.y+= max(temp,0.0);
        }
        temp = o.ks[2]*(l.intensity[2])*(pow(dot(r,v),o.exp));
        if(dot(r,v)>0){
                ans.z+= max(temp,0.0);
        }
        return ans;
}

void specular(lightsource l, vector3d normal, vector3d v,Object o){
        vector3d ans(0,0,0);
        vector3d r = addvec(l.direction,scalarmulti(-2*dot(l.direction,normal),normal));
        for(int i=0;i<3;i++){
                double temp = o.ks[i]*(l.intensity[i])*(pow(dot(r,v),o.exp));
                ans.x= max(temp,0.0);
        }
}


bool intersectPlane(Plane p, vector3d origin,vector3d direct , double &t)
{
        vector3d normal =  norm(cross_prod(subvec(p.x,p.y),subvec(p.x,p.z)));
        float denom = dot(normal, direct);
        if (fabs(denom) > 1e-3) {
                vector3d a = subvec(p.x , origin);
                t = dot(a, normal) / denom;
                vector3d point;
                point  = addvec(origin, scalarmulti(t,direct));
                if(t>1e-3){
                        if (fabs( (area(point,p.x,p.y)+area(point,p.z,p.y)+area(point,p.x,p.z)) - area(p.z,p.x,p.y))<(pow(10,-2))) {
                                return true;
                        }
                        else{
                                return false;
                        }
                }
        }
        return false;
}


bool intersectsphere(Sphere * c,double *t,Ray * r)
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
                double d=(r->org.x-c->c.x)*(r->org.x-c->c.x)+(r->org.y-c->c.y)*(r->org.y-c->c.y)+(r->org.z-c->c.z)*(r->org.z-c->c.z)-(c->r)*(c->r);
                double t0=(-b-sqrt(b*b-4.0*d))/2.0;
                double t1=(-b+sqrt(b*b-4.0*d))/2.0;
                if((min(t0,t1))>pow(10,-3)){
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
        glClear(GL_COLOR_BUFFER_BIT);
        float* pixels = new float[viewx*viewy*3];
        glDrawPixels(viewx,viewy,GL_RGB,GL_FLOAT,pixels);
        glFlush();        
        glutSwapBuffers();

}


void reshape(int x, int y)
{
        if (y == 0 || x == 0) return;   
        glMatrixMode(GL_PROJECTION);  
        glLoadIdentity(); 
        glOrtho(-10.0, 10.0, -10.0, 10.0, -100.0, 100.0);
        glMatrixMode(GL_MODELVIEW);
        glViewport(0,0,x,y);
} 


vector3d viewporttovcs(float x,float y)
{
        x=x-viewx/2.0;y=viewy/2.0-y;
        float x1  =  (float)(x)*((float)viewl)/((float)viewx);
        float y1 = (float)(y)*((float)viewb)/((float)viewy);
        vector3d f=addvec(scalarmulti(x1,u),scalarmulti(y1,v));
        f=addvec(f,vcs);
        return f;
}
vector3d reflected(vector3d ,vector3d );
vector3d refracted(vector3d ,double,vector3d );


vector3d rayTracerRecurse(Ray rt,  int level){
        double t = INT_MAX;
        int h=obj.size();
        double k=t;
        vector3d normal;vector3d retr(0,0,0);
        int intersectid;
        for(int j=0;j<h;j++)
        {
                Object a=obj[j];bool p;
                if(a.type==1)
                {
                        p=intersectsphere(&(a.s),&k,&rt);
                }
                else
                        p=intersectPlane((a.p),rt.org,rt.dir,k);
                if(k<t && k>0 && p==true){
                        t=k;
                        intersectid=j;
                } 
        }
        if(t!=INT_MAX){
                Object a = obj[intersectid];
                if(a.type==1){
                        normal = norm(subvec(addvec(rt.org,scalarmulti(t,rt.dir)),a.s.c));
                }
                else{
                        normal = norm(cross_prod(subvec(a.p.x,a.p.y),subvec(a.p.x,a.p.z)));
                }      

                //////////////// SHadow///////////////////////
                for(int r=0;r<lightsources.size();r++)
                {
                        bool l =false;
                        vector3d src=addvec(rt.org,scalarmulti(t,rt.dir));
                        Ray ty(src,scalarmulti(-1,lightsources[r].direction));
                        int ghj;
                        for(int j=0;j<h;j++)
                        {
                                if(j==intersectid){
                                        continue;
                                }
                                Object a=obj[j];
                                if(a.type==1)
                                {
                                        l=intersectsphere(&(a.s),&k,&ty);
                                        if(l)
                                        { 
                                                ghj=j;
                                                break;
                                        }

                                }
                                else
                                {
                                        l=intersectPlane((a.p),ty.org,ty.dir,k);
                                        if(l)
                                        {   ghj=j;break;}
                                }
                        }

                        if(!l)
                        {

                                retr = addvec( retr, scalarmulti(1-a.kreflec-a.krefrac,diffuserecur(lightsources[r],normal,a)));
                                retr = addvec(retr, scalarmulti(1-a.kreflec-a.krefrac,specularrecur(lightsources[r],normal,scalarmulti(-1,rt.dir),a)));
                        }
                        else{
                                retr = addvec( retr, scalarmulti((1-a.kreflec-a.krefrac)*(obj[ghj].krefrac),diffuserecur(lightsources[r],normal,a)));
                                retr = addvec(retr, scalarmulti((1-a.kreflec-a.krefrac)*(obj[ghj].krefrac),specularrecur(lightsources[r],normal,scalarmulti(-1,rt.dir),a)));
                        }

                }
                if(level<MAX_DEPTH){

                        retr=addvec(retr,scalarmulti(a.kreflec,rayTracerRecurse(Ray(addvec(rt.org,scalarmulti(t,rt.dir)),reflected(normal,rt.dir)),level+1)));
                        retr=addvec(retr,scalarmulti(a.krefrac,rayTracerRecurse(Ray(addvec(rt.org,scalarmulti(t,rt.dir)),refracted(normal,a.n,rt.dir)),level+1)));
                }
                retr=addvec(retr,scalarmulti(1-a.kreflec-a.krefrac,ambientrecur(a)));
                return retr;      
        }
        else
                return vector3d(ambientintensity[0],ambientintensity[1],ambientintensity[2]);
}

vector3d reflected(vector3d normal,vector3d l)
{
        return subvec(l,scalarmulti(2*dot(normal,l),normal));
}

vector3d refracted(vector3d normal,double n,vector3d inc)
{

        double cos1=dot(inc,normal);double n1,n2;
        if(cos1>0)
        {
                n1=n;n2=1;
        }
        else
        {
                n1=1;n2=n;cos1=-cos1;
        }
        double cos2=1.0-(n1*n1*(1-cos1*cos1)/(n2*n2));
        if(cos2<0)
        {
                return addvec(inc,scalarmulti(2*dot(inc,normal),normal));
        }
        else
        {
                cos2=sqrt(cos2);

                return subvec(scalarmulti(n1/n2,addvec(inc,scalarmulti(cos1,normal))),scalarmulti(cos2,normal));
        }

}

void raytracer()
{
        vector3d retr(0,0,0);
        glMatrixMode(GL_MODELVIEW);
        glClear(GL_COLOR_BUFFER_BIT);
        float* pixels = new float[viewx*viewy*3];
        for(int i=0;i<viewx*viewy;i++)
        {
                retr = vector3d(0,0,0);
                for(int op=-1;op<2;op++)
                {for(int yu=-1;yu<2;yu++)
                        {

                                vector3d f=viewporttovcs((float)(i%viewx)+op*0.33,(float)(i/viewx)+yu*0.33);
                                Ray rt(eye,norm(subvec(f,eye)));
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
                                                p=intersectsphere(&(a.s),&k,&rt);
                                        }
                                        else{
                                                p=intersectPlane((a.p),rt.org,rt.dir,k);        
                                        }
                                        if(k<t && k>0 && p){
                                                t=k;
                                                intersectid=j;
                                        } 
                                }
                                if(t!=INT_MAX){
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
                                        retr = addvec(retr,ambientrecur(a));

                                        //////////////// Shadow///////////////////////
                                        for(int r=0;r<lightsources.size();r++)
                                        {
                                                bool l =false;
                                                vector3d src=addvec(rt.org,scalarmulti(t,rt.dir));
                                                Ray ty(src,scalarmulti(-1,lightsources[r].direction));
                                                int ghj;
                                                for(int j=0;j<h;j++)
                                                {
                                                        if(j==intersectid){
                                                                continue;
                                                        }
                                                        Object a=obj[j];
                                                        if(a.type==1)
                                                        {
                                                                l=intersectsphere(&(a.s),&k,&ty);
                                                                if(l)
                                                                { ghj=j;
                                                                        break;}

                                                        }
                                                        else
                                                        {
                                                                l=intersectPlane((a.p),ty.org,ty.dir,k);
                                                                if(l)
                                                                {ghj=j;break;}
                                                        }
                                                }
                                                if(!l)
                                                {
                                                        retr = addvec( retr, scalarmulti(1-a.kreflec-a.krefrac,diffuserecur(lightsources[r],normal,a)));
                                                        retr = addvec(retr, scalarmulti(1-a.kreflec-a.krefrac,specularrecur(lightsources[r],normal,scalarmulti(-1,rt.dir),a)));
                                                }
                                                else{
                                                        retr = addvec( retr, scalarmulti((1-a.kreflec-a.krefrac)*obj[ghj].krefrac,diffuserecur(lightsources[r],normal,a)));
                                                        retr = addvec(retr, scalarmulti((1-a.kreflec-a.krefrac)*obj[ghj].krefrac,specularrecur(lightsources[r],normal,scalarmulti(-1,rt.dir),a)));
                                                }


                                        }
                                        retr=addvec(retr,scalarmulti(a.kreflec,rayTracerRecurse(Ray(addvec(rt.org,scalarmulti(t,rt.dir)),reflected(normal,rt.dir)),1)));

                                        retr=addvec(retr,scalarmulti(a.krefrac,rayTracerRecurse(Ray(addvec(rt.org,scalarmulti(t,rt.dir)),refracted(normal,a.n,rt.dir)),1)));
                                }}
                }
                retr=scalarmulti(1.0/9.0,retr);     
                pixels[3*i]=retr.x;pixels[3*i+1]=retr.y;pixels[3*i+2]=retr.z;

        }
        glDrawPixels(viewx,viewy,GL_RGB,GL_FLOAT,pixels);
        glFlush();        
        glutSwapBuffers();
}




int main (int argc, char **argv)
{

        glutInit(&argc, argv); 
        scanf("%d",&viewx);scanf("%d",&viewy);scanf("%d",&viewl);scanf("%d",&viewb);
        scanf("%lf",&u.x);scanf("%lf",&u.y);scanf("%lf",&u.z);
        scanf("%lf",&v.x);scanf("%lf",&v.y);scanf("%lf",&v.z);
        scanf("%lf",&vcs.x);scanf("%lf",&vcs.y);scanf("%lf",&vcs.z);
        scanf("%lf",&disteye);
        glutInitWindowSize(viewx,viewy);
        glutCreateWindow("Solid Sphere");

        n=cross_prod(u,v);
        n=norm(n);


        eye=subvec(vcs,scalarmulti(disteye,n));
        int no;
        scanf("%d",&no);
        for(int g=0;g<no;g++)
        {
                Object a;
                int q;
                scanf("%d",&q);
                if(q==1)
                {
                        a.type=1;
                        scanf("%lf %lf %lf %lf",&a.s.c.x,&a.s.c.y,&a.s.c.z,&a.s.r);
                }
                else
                {
                        a.type=2;
                        scanf("%lf %lf %lf",&a.p.x.x,&a.p.x.y,&a.p.x.z);
                        scanf("%lf %lf %lf",&a.p.y.x,&a.p.y.y,&a.p.y.z);
                        scanf("%lf %lf %lf",&a.p.z.x,&a.p.z.y,&a.p.z.z);  

                }
                scanf("%lf %lf %lf",&a.ka[0],&a.ka[1],&a.ka[2]);
                scanf("%lf %lf %lf",&a.kd[0],&a.kd[1],&a.kd[2]);
                scanf("%lf %lf %lf",&a.ks[0],&a.ks[1],&a.ks[2]);
                scanf("%lf %lf %lf %lf",&a.exp,&a.kreflec,&a.krefrac,&a.n);
                obj.push_back(a);
        }
        scanf("%lf %lf %lf",&ambientintensity[0],&ambientintensity[1],&ambientintensity[2]);
        scanf("%d",&no);
        for(int gh=0;gh<no;gh++)
        {
                lightsource ll;
                scanf("%lf %lf %lf",&ll.intensity[0],&ll.intensity[1],&ll.intensity[2]);
                scanf("%lf %lf %lf",&ll.direction.x,&ll.direction.y,&ll.direction.z);
                lightsources.push_back(ll);
        }
        glutDisplayFunc(raytracer);
        glutReshapeFunc(reshape);
        glutMainLoop();
        return 0;
}

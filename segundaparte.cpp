#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

double der1(double t, double y, double x); 
double der2(double t, double v, double x);
void euler(double dt, string name);
void rk(double dt, string name);

#define K 100.0
#define M  2.0;

int main ()
{
    double t1 = 0.01;
    euler(t1, "euler.dat");
    rk(t1, "rk.dat");
    
    cout<<"Solucion en terminos de seno y coseno  (oscilador armonico) ";
    return 0;
}

double der1(double t, double v, double x)
{
  return v;
}


double der2(double t, double v, double x)
{
  return (-K*x);
}


void euler(double delt, string name)
{
    int step = 1000;
    
    double x[step];
    double vx[step];
    ofstream outfile;
    outfile.open(name);
    
    
    x[0] = 1;
    vx[0] = 0;
    
    outfile<< x[0]<<" "<<vx[0]<<" "<<endl;
    
    double t2=0.0;
    
    for (int i=1; i<=20; i++)
    {
        x[i] = (vx[i-1]*delt)+x[i-1];
        vx[i] = (x[i-1]*delt)+vx[i-1]; 
        
        outfile<< x[i]<<" "<<vx[i]<<" "<<endl;
        t2 = t2 +delt;
    }
    outfile.close();
    
}

void rk(double dt, string name)
{
    double t;
    double xpres; 
    double vxpres;
    double ypres; 
    double vypres;
    double xfut;
    double vxfut;
    double k0x;
    double k0vx;
    double k1x;
    double k1vx;
    double k2x;
    double k2vx;
    double k3x;
    double k3vx;
    double kx;
    double kvx;
    double dvx;
    double dvx1;
  
    ofstream outfile;
    outfile.open(name);
    
    xpres = 1;
    vxpres = 0;
    t=0.0;
    int n = 20;
    
    for (int i=0; i<=n;i++)
    {
        dvx = (ypres*dt)+vypres;      
        
        k0x = vxpres;
        k0vx = dvx;
        
        xfut = xpres + ((0.5 * dt) * k0x);
        vxfut = vxpres + ((0.5 * dt) * k0vx);

        k1x = vxfut;
        k1vx = dvx;
        
        xfut = xpres + ((0.5 * dt) * k1x);
        vxfut = vxpres + ((0.5 * dt) * k1vx);
        k2x = vxfut;
        k2vx = dvx;
        
       
        xfut = xpres + (dt * k2x);
        vxfut = vxpres + (dt * k2vx);
        k3x = vxfut;
        k3vx = dvx;
        
        
        kx  = (k0x/6.0) + (k1x/3.0) + (k2x/3.0) + (k3x/6.0);
        kvx = (k0vx/6.0) + (k1vx/3.0) + (k2vx/3.0) + (k3vx/6.0);
        
        xfut = xpres + (dt * kx);
        vxfut = vxpres + (dt * kvx);
        t = t + dt;
        
        xpres = xfut;
        vxpres = vxfut;
     
    }
    
    outfile<<xfut<<" "<<vxfut<<" "<<endl;
    outfile.close();
}
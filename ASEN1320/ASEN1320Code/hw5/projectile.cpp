#include<iostream> 
#include<cmath>
using namespace std; 
  
// Set up global variables
const double g  = 9.81;  // gravitational acceleration of the Earth [m/s2]
const double PI = 3.141592653589793238463; 
double mass     = 50.0;              // mass of the projectile [kg]
double radius   = 0.25;              // radius of the projectile [m]
double Cd       = 0.5;               // drag coefficent for the projectile 
double rho      = 1.275;             // density of air at 15 deg [kg/m3]
double A        = PI*radius*radius;  // cross-sectional area of the projectile  
 
// Function prototypes 
double Fx(double vx, double vy); // call-by-value
double Fy(double vx, double vy); // call-by-value
double euler(double &t, double &x, double y, double vx, double vy, double DT); // call-by-reference

int main() 
{ 
    double v0 = 10.0, angle = 45.0, DT = 0.01, y0 = 1.0; // Initial conditions ** Change these for testing on Gradescope** 
    double t = 0.0, x = 0.0, vx = v0*cos(angle*PI/180.0), vy = v0*sin(angle*PI/180.0);
    
    cout << "Initial speed (m/s): " << v0 << endl;
    cout << "Initial angle (deg): " << angle << endl;
    cout << "Initial height (m): "  << y0    << endl;
    cout << "Time step (sec): "     << DT    << endl;
    
    double ymax = euler(t, x, y0, vx, vy, DT); // call-by-reference to euler 
    
    cout << "Peak height (m): "   << ymax << endl;  // ** Check this ymax number (4/16 pt)
    cout << "Range (m): " << x  << endl;   
    cout << "Travel time (sec): "       << t  << endl;
 
    return 0; 
}  
  
double Fx(double vx, double vy)   // Fxfunction  ** Check this via unit test (3/16 pt)
{
    double theta = atan2(vy,vx); 
    double D     = 0.5*rho*Cd*A*(vx*vx + vy*vy);
    return (-D*cos(theta)/mass);     // returning value
} 

double Fy(double vx, double vy)   // Fy function ** Check this via unit test (3/16 pt)
{ 
    double theta = atan2(vy,vx);
    double D     = 0.5*rho*Cd*A*(vx*vx + vy*vy);
    return (-D*sin(theta)/mass - g); // returning value
} 
  
double euler(double &t, double &x, double y, double vx, double vy, double DT) // euler integrator function  ** Check this via unit test (6/16 pt)
{ 
    double ymax = y;
    while (y > 0.0) {   // Iterating till hitting the ground (y (height) > 0)
        t += DT;
        x += vx*DT; 
        y += vy*DT;
       vx += Fx(vx,vy)*DT;  // call-by-value to Fx
       vy += Fy(vx,vy)*DT;  // call-by-value to Fy
       if (ymax < y ) ymax = y; 
    } 
   return ymax;
} 

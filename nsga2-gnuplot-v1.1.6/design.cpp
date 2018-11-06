// Design.cpp
// All the new and tested functions are here
#include <iostream>
#include <stdio.h>
#include <Eigen/Dense>
#include <cmath>
#include <gsl/gsl_poly.h>
#include <limits>

using namespace Eigen;
//This function checks if the assumed point is in the workspace are not
bool inWorkspace(float l1, float r1, float l0, float k){
  float d1=0;
  d1 = pow(pow(l0/2,2) + pow(k,2),0.5);
  if(d1<l1+r1 && d1>abs(l1-r1)){
    return true;
  }
  else{
    return false;
  }
}

// Computes the radius with respect to the workspace boundaries
MatrixXd workspacerad(float l1, float r1, float l0, float k){
  float d1 = pow(pow(l0/2,2) + pow(k,2),0.5);
  MatrixXd radwork(2, 1);
  radwork << l1+r1-d1, d1-fabs(l1-r1);
  return radwork;
}

// This function returns the IK solution of a chosen branch
MatrixXd ikcal(float l, float r, float l0, float x, float y){
  MatrixXd iksol(2,1);
    iksol << -acos((-pow(l,2) + pow(r,2) - pow(x,2) - pow(y,2))/(2.*pow(r,2)*pow((pow(l,2)*(pow(x,2) + pow(y,2)))/pow(r,4), 0.5))) + atan2(-((l*y)/pow(r,2)),-((l*x)/pow(r,2))),acos((-pow(l,2) - pow(l0,2) + pow(r,2) + 2*l0*x - pow(x,2) - pow(y,2))/(2.*pow(r,2)*pow((pow(l,2)*(pow(l0,2) - 2*l0*x + pow(x,2) + pow(y,2)))/pow(r,4), 0.5))) + atan2(-((l*y)/pow(r,2)),(l*(l0 - x))/pow(r,2));
    return iksol;
}

// Returns the mom val at the given point using ik
double momcal(float l, float r, float l0, float x, float y, float th1, float th2){
  double mom = pow(l,4)*pow(y*cos(th1) - x*sin(th1),2)*pow(y*cos(th2) + (l0 - x)*sin(th2),2)*pow(-(l0*y) + l*y*cos(th1) - l*y*cos(th2) + l*l0*sin(th1) - l*x*sin(th1) + pow(l,2)*sin(th1 - th2) + l*x*sin(th2),-2);
  return mom;
}

// Returns the kappainv vals at a given point using ik
double kinvcal(float l, float r, float l0, float x, float y, float th1, float th2){
  double ki = (pow(y*cos(th1) - x*sin(th1),-2)*pow(y*cos(th2) + (l0 - x)*sin(th2),-2)*(pow(l0 - x + l*cos(th2),2)*pow(y*cos(th1) - x*sin(th1),2) + pow(y*cos(th1) - x*sin(th1),2)*pow(y - l*sin(th2),2) + pow(x - l*cos(th1),2)*pow(y*cos(th2) + (l0 - x)*sin(th2),2) + pow(y - l*sin(th1),2)*pow(y*cos(th2) + (l0 - x)*sin(th2),2))*pow(-(l0*y) + l*y*cos(th1) - l*y*cos(th2) + l*l0*sin(th1) - l*x*sin(th1) + pow(l,2)*sin(th1 - th2) + l*x*sin(th2),-2)*(pow(l0,2)*pow(x,2)*pow(sin(th1),2) - 2*l0*pow(x,3)*pow(sin(th1),2) + pow(x,4)*pow(sin(th1),2) + pow(x,2)*pow(y,2)*pow(sin(th1),2) + pow(l0,2)*pow(x,2)*pow(sin(th2),2) - 2*l0*pow(x,3)*pow(sin(th2),2) + pow(x,4)*pow(sin(th2),2) - 2*l0*x*pow(y,2)*pow(sin(th2),2) + pow(l0,2)*pow(y,2)*pow(sin(th2),2) + pow(x,2)*pow(y,2)*pow(sin(th2),2) - 2*l0*x*pow(l,2)*pow(sin(th1),2)*pow(sin(th2),2) + pow(l,2)*pow(l0,2)*pow(sin(th1),2)*pow(sin(th2),2) + 2*pow(l,2)*pow(x,2)*pow(sin(th1),2)*pow(sin(th2),2) + 4*l*l0*x*y*pow(sin(th2),2)*sin(th1) - 2*l*y*pow(l0,2)*pow(sin(th2),2)*sin(th1) - 2*l*y*pow(x,2)*pow(sin(th2),2)*sin(th1) + pow(cos(th2),2)*(pow(y,2)*(pow(x,2) + pow(y,2)) + pow(l,2)*(pow(x,2) + pow(y,2))*pow(sin(th1),2) - 2*l*pow(y,3)*sin(th1)) - 2*l*y*pow(x,2)*pow(sin(th1),2)*sin(th2) - 2*x*cos(th1)*(l*pow(l0 - x,2)*pow(sin(th2),2) + l*y*pow(cos(th2),2)*(y + l*sin(th1)) + 2*l*(l0 - x)*y*cos(th2)*(sin(th1) + sin(th2)) + y*sin(th1)*(-2*l0*x + pow(l0,2) + pow(x,2) + pow(y,2) + pow(l,2)*pow(sin(th2),2) - 2*l*y*sin(th2))) - 2*cos(th2)*(x*y*(pow(x,2) + pow(y,2))*sin(th2) + 2*l*(l0 - x)*pow(y,2)*sin(th1)*sin(th2) + l*x*pow(sin(th1),2)*(x*(-l0 + x) + l*y*sin(th2))) + l0*y*pow(x,2)*sin(2*th2) + l0*pow(y,3)*sin(2*th2) + l0*y*pow(l,2)*pow(sin(th1),2)*sin(2*th2) + pow(cos(th1),2)*(2*pow(l,2)*pow(y,2)*pow(cos(th2),2) + pow(l,2)*(-2*l0*x + pow(l0,2) + pow(x,2) + pow(y,2))*pow(sin(th2),2) - 2*l*pow(y,3)*sin(th2) - 2*l*y*cos(th2)*((-l0 + x)*y + l*x*sin(th2)) + y*(y*(-2*l0*x + pow(l0,2) + pow(x,2) + pow(y,2)) + l0*pow(l,2)*sin(2*th2)))))/4.;
  return ki;
}

// Returns the radius for other cases like FSM other factors of the singularity manifold
MatrixXd radfsm(float l1, float r1, float l0, float k){
  float pen = std::numeric_limits<float>::max();
  MatrixXd radvals(25, 1);
  radvals << pow(pow(k,2) + pow(l0,2)/4.,0.5),pow(pow(l0,2) + pow(-2*k - 2*r1 + pow(-pow(l0,2) + 4*pow(l1,2),0.5),2),0.5)/2.,pow(pow(l0,2) + pow(2*k + 2*r1 + pow(-pow(l0,2) + 4*pow(l1,2),0.5),2),0.5)/2.,pow(pow(l0,2)/4. + pow(k - r1 + pow(-pow(l0,2) + 4*pow(l1,2),0.5)/2.,2),0.5),pow(pow(l0,2)/4. + pow(-k + r1 + pow(-pow(l0,2) + 4*pow(l1,2),0.5)/2.,2),0.5),pow(-(r1*(l0 + r1)) + pow(k,2) + pow(l1,2) + k*pow(-4*l0*r1 - pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5),0.5),pow(-(r1*(l0 + r1)) + pow(k,2) + pow(l1,2) - k*pow(-4*l0*r1 - pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5),0.5),pow((l0 - r1)*r1 + pow(k,2) + pow(l1,2) + k*pow(4*l0*r1 - pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5),0.5),pow((l0 - r1)*r1 + pow(k,2) + pow(l1,2) - k*pow(4*l0*r1 - pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5),0.5),pow(pow(k,2) + pow(l0 - pow(pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5),2)/4.,0.5),pow(pow(l0 - pow(pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5),2)/4. + pow(4*k + pow(2,0.5)*pow(-5*pow(l0,2) - 4*pow(l1,2) + 4*pow(r1,2) - pow(9*pow(l0,4) + 8*pow(l0,2)*(5*pow(l1,2) + 3*pow(r1,2)) + 16*pow(pow(l1,2) - pow(r1,2),2),0.5),0.5),2)/16.,0.5),pow(pow(l0 - pow(pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5),2)/4. + pow(k - (pow(2,-0.5)*pow(-5*pow(l0,2) - 4*pow(l1,2) + 4*pow(r1,2) - pow(9*pow(l0,4) + 8*pow(l0,2)*(5*pow(l1,2) + 3*pow(r1,2)) + 16*pow(pow(l1,2) - pow(r1,2),2),0.5),0.5))/2.,2),0.5),pow(pow(l0 - pow(pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5),2)/4. + pow(4*k + pow(2,0.5)*pow(-5*pow(l0,2) - 4*pow(l1,2) + 4*pow(r1,2) + pow(9*pow(l0,4) + 8*pow(l0,2)*(5*pow(l1,2) + 3*pow(r1,2)) + 16*pow(pow(l1,2) - pow(r1,2),2),0.5),0.5),2)/16.,0.5),pow(pow(l0 - pow(pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5),2)/4. + pow(k - pow((-5*pow(l0,2))/8. + (-4*pow(l1,2) + 4*pow(r1,2) + pow(9*pow(l0,4) + 8*pow(l0,2)*(5*pow(l1,2) + 3*pow(r1,2)) + 16*pow(pow(l1,2) - pow(r1,2),2),0.5))/8.,0.5),2),0.5),pow(pow(k,2) + pow(l0 + pow(pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5),2)/4.,0.5),pow(pow(l0 + pow(pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5),2)/4. + pow(4*k + pow(2,0.5)*pow(-5*pow(l0,2) - 4*pow(l1,2) + 4*pow(r1,2) - pow(9*pow(l0,4) + 8*pow(l0,2)*(5*pow(l1,2) + 3*pow(r1,2)) + 16*pow(pow(l1,2) - pow(r1,2),2),0.5),0.5),2)/16.,0.5),pow(pow(l0 + pow(pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5),2)/4. + pow(k - (pow(2,-0.5)*pow(-5*pow(l0,2) - 4*pow(l1,2) + 4*pow(r1,2) - pow(9*pow(l0,4) + 8*pow(l0,2)*(5*pow(l1,2) + 3*pow(r1,2)) + 16*pow(pow(l1,2) - pow(r1,2),2),0.5),0.5))/2.,2),0.5),pow(pow(l0 + pow(pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5),2)/4. + pow(4*k + pow(2,0.5)*pow(-5*pow(l0,2) - 4*pow(l1,2) + 4*pow(r1,2) + pow(9*pow(l0,4) + 8*pow(l0,2)*(5*pow(l1,2) + 3*pow(r1,2)) + 16*pow(pow(l1,2) - pow(r1,2),2),0.5),0.5),2)/16.,0.5),pow(pow(l0 + pow(pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5),2)/4. + pow(k - pow((-5*pow(l0,2))/8. + (-4*pow(l1,2) + 4*pow(r1,2) + pow(9*pow(l0,4) + 8*pow(l0,2)*(5*pow(l1,2) + 3*pow(r1,2)) + 16*pow(pow(l1,2) - pow(r1,2),2),0.5))/8.,0.5),2),0.5),pow(pow(k,2) + pow(l0 - pow(pow(l0,2) - 4*pow(l1,2) + 4*pow(r1,2),0.5),2)/4.,0.5),pow(pow(k + pow(-pow(l0,2) + 4*pow(l1,2),0.5),2) + pow(l0 - pow(pow(l0,2) - 4*pow(l1,2) + 4*pow(r1,2),0.5),2)/4.,0.5),pow(pow(k - pow(-pow(l0,2) + 4*pow(l1,2),0.5),2) + pow(l0 - pow(pow(l0,2) - 4*pow(l1,2) + 4*pow(r1,2),0.5),2)/4.,0.5),pow(pow(k,2) + pow(l0 + pow(pow(l0,2) - 4*pow(l1,2) + 4*pow(r1,2),0.5),2)/4.,0.5),pow(pow(k + pow(-pow(l0,2) + 4*pow(l1,2),0.5),2) + pow(l0 + pow(pow(l0,2) - 4*pow(l1,2) + 4*pow(r1,2),0.5),2)/4.,0.5),pow(pow(k - pow(-pow(l0,2) + 4*pow(l1,2),0.5),2) + pow(l0 + pow(pow(l0,2) - 4*pow(l1,2) + 4*pow(r1,2),0.5),2)/4.,0.5);
  // Send only the real ones by filtering the nans as 0's
  return radvals.array().isNaN().select(pen, radvals);
}

// Returns <0 implies real lines not possible hence acnode, 0 -> cusp
float doublerootcheck(float l1, float r1, float l0, float p, float q){
  float drc = -24*p*pow(l0,7) + pow(l0,8) + 232*p*pow(l0,5)*pow(l1,2) - 16*pow(l0,6)*pow(l1,2) - 544*p*pow(l0,3)*pow(l1,4) + 68*pow(l0,4)*pow(l1,4) + 256*l0*p*pow(l1,6) - 64*pow(l0,2)*pow(l1,6) + 16*pow(l1,8) + 204*pow(l0,6)*pow(p,2) - 1288*pow(l0,4)*pow(l1,2)*pow(p,2) + 1600*pow(l0,2)*pow(l1,4)*pow(p,2) - 256*pow(l1,6)*pow(p,2) - 936*pow(l0,5)*pow(p,3) + 3648*pow(l0,3)*pow(l1,2)*pow(p,3) - 2112*l0*pow(l1,4)*pow(p,3) + 2628*pow(l0,4)*pow(p,4) - 5664*pow(l0,2)*pow(l1,2)*pow(p,4) + 1056*pow(l1,4)*pow(p,4) - 4608*pow(l0,3)*pow(p,5) + 4608*l0*pow(l1,2)*pow(p,5) + 4896*pow(l0,2)*pow(p,6) - 1536*pow(l1,2)*pow(p,6) - 2880*l0*pow(p,7) + 720*pow(p,8) + 936*p*pow(l0,5)*pow(q,2) - 96*pow(l0,6)*pow(q,2) - 2304*p*pow(l0,3)*pow(l1,2)*pow(q,2) + 296*pow(l0,4)*pow(l1,2)*pow(q,2) + 960*l0*p*pow(l1,4)*pow(q,2) - 256*pow(l0,2)*pow(l1,4)*pow(q,2) - 256*pow(l1,6)*pow(q,2) - 4104*pow(l0,4)*pow(p,2)*pow(q,2) + 6912*pow(l0,2)*pow(l1,2)*pow(p,2)*pow(q,2) - 960*pow(l1,4)*pow(p,2)*pow(q,2) + 10368*pow(l0,3)*pow(p,3)*pow(q,2) - 9216*l0*pow(l1,2)*pow(p,3)*pow(q,2) - 15264*pow(l0,2)*pow(p,4)*pow(q,2) + 4608*pow(l1,2)*pow(p,4)*pow(q,2) + 12096*l0*pow(p,5)*pow(q,2) - 4032*pow(p,6)*pow(q,2) + 4608*p*pow(l0,3)*pow(q,4) - 540*pow(l0,4)*pow(q,4) - 4608*l0*p*pow(l1,2)*pow(q,4) + 1056*pow(l0,2)*pow(l1,2)*pow(q,4) + 1056*pow(l1,4)*pow(q,4) - 14112*pow(l0,2)*pow(p,2)*pow(q,4) + 4608*pow(l1,2)*pow(p,2)*pow(q,4) + 19008*l0*pow(p,3)*pow(q,4) - 9504*pow(p,4)*pow(q,4) + 4032*l0*p*pow(q,6) - 864*pow(l0,2)*pow(q,6) - 1536*pow(l1,2)*pow(q,6) - 4032*pow(p,2)*pow(q,6) + 720*pow(q,8) - 160*p*pow(l0,5)*pow(r1,2) + 12*pow(l0,6)*pow(r1,2) + 896*p*pow(l0,3)*pow(l1,2)*pow(r1,2) - 96*pow(l0,4)*pow(l1,2)*pow(r1,2) - 768*l0*p*pow(l1,4)*pow(r1,2) + 176*pow(l0,2)*pow(l1,4)*pow(r1,2) - 64*pow(l1,6)*pow(r1,2) + 976*pow(l0,4)*pow(p,2)*pow(r1,2) - 3008*pow(l0,2)*pow(l1,2)*pow(p,2)*pow(r1,2) + 768*pow(l1,4)*pow(p,2)*pow(r1,2) - 3168*pow(l0,3)*pow(p,3)*pow(r1,2) + 4224*l0*pow(l1,2)*pow(p,3)*pow(r1,2) + 5424*pow(l0,2)*pow(p,4)*pow(r1,2) - 2112*pow(l1,2)*pow(p,4)*pow(r1,2) - 4608*l0*pow(p,5)*pow(r1,2) + 1536*pow(p,6)*pow(r1,2) + 2592*p*pow(l0,3)*pow(q,2)*pow(r1,2) - 368*pow(l0,4)*pow(q,2)*pow(r1,2) - 1920*l0*p*pow(l1,2)*pow(q,2)*pow(r1,2) + 576*pow(l0,2)*pow(l1,2)*pow(q,2)*pow(r1,2) + 768*pow(l1,4)*pow(q,2)*pow(r1,2) - 7200*pow(l0,2)*pow(p,2)*pow(q,2)*pow(r1,2) + 1920*pow(l1,2)*pow(p,2)*pow(q,2)*pow(r1,2) + 9216*l0*pow(p,3)*pow(q,2)*pow(r1,2) - 4608*pow(p,4)*pow(q,2)*pow(r1,2) + 4608*l0*p*pow(q,4)*pow(r1,2) - 1104*pow(l0,2)*pow(q,4)*pow(r1,2) - 2112*pow(l1,2)*pow(q,4)*pow(r1,2) - 4608*pow(p,2)*pow(q,4)*pow(r1,2) + 1536*pow(q,6)*pow(r1,2) - 352*p*pow(l0,3)*pow(r1,4) + 28*pow(l0,4)*pow(r1,4) + 768*l0*p*pow(l1,2)*pow(r1,4) - 160*pow(l0,2)*pow(l1,2)*pow(r1,4) + 96*pow(l1,4)*pow(r1,4) + 1408*pow(l0,2)*pow(p,2)*pow(r1,4) - 768*pow(l1,2)*pow(p,2)*pow(r1,4) - 2112*l0*pow(p,3)*pow(r1,4) + 1056*pow(p,4)*pow(r1,4) + 960*l0*p*pow(q,2)*pow(r1,4) - 320*pow(l0,2)*pow(q,2)*pow(r1,4) - 768*pow(l1,2)*pow(q,2)*pow(r1,4) - 960*pow(p,2)*pow(q,2)*pow(r1,4) + 1056*pow(q,4)*pow(r1,4) - 256*l0*p*pow(r1,6) + 48*pow(l0,2)*pow(r1,6) - 64*pow(l1,2)*pow(r1,6) + 256*pow(p,2)*pow(r1,6) + 256*pow(q,2)*pow(r1,6) + 16*pow(r1,8);
  return drc;
}

//Solves the univariate in x
int polysolvex(float l1 , float r1 ,float l0 ,float k, double *u){
  /* coefficients of the univariate obtained */
  double a[9] = {64*pow(l0,10)*pow(l1,4)*pow(pow(l0,2) - 4*pow(r1,2),4) + pow(k,2)*pow(l0,4)*(81*pow(l0,16) - 288*pow(l0,14)*(2*pow(l1,2) - 3*pow(r1,2)) + 32*pow(l0,12)*(23*pow(l1,4) - 258*pow(l1,2)*pow(r1,2) + 54*pow(r1,4)) + 256*pow(l0,10)*(73*pow(l1,4)*pow(r1,2) - 82*pow(l1,2)*pow(r1,4) - 30*pow(r1,6)) + 256*pow(l0,8)*(pow(l1,8) - 36*pow(l1,6)*pow(r1,2) + 618*pow(l1,4)*pow(r1,4) + 196*pow(l1,2)*pow(r1,6) - 74*pow(r1,8)) - 4096*pow(l0,6)*pow(r1,2)*(pow(l1,8) + 44*pow(l1,6)*pow(r1,2) + 66*pow(l1,4)*pow(r1,4) - 16*pow(l1,2)*pow(r1,6) - 10*pow(r1,8)) + 8192*pow(l0,4)*pow(r1,4)*(3*pow(l1,8) - 28*pow(l1,6)*pow(r1,2) + 49*pow(l1,4)*pow(r1,4) - 30*pow(l1,2)*pow(r1,6) + 6*pow(r1,8)) - 65536*pow(l0,2)*(pow(l1,4) - 2*pow(l1,2)*pow(r1,2) + 2*pow(r1,4))*pow(r1,6)*pow(pow(l1,2) - pow(r1,2),2) + 65536*pow(r1,8)*pow(pow(l1,2) - pow(r1,2),4)) + 8*pow(k,4)*pow(l0,2)*(81*pow(l0,16) + pow(l0,14)*(-396*pow(l1,2) + 864*pow(r1,2)) + 48*pow(l0,12)*(8*pow(l1,4) - 149*pow(l1,2)*pow(r1,2) + 36*pow(r1,4)) - 64*pow(l0,10)*(5*pow(l1,6) - 289*pow(l1,4)*pow(r1,2) + 371*pow(l1,2)*pow(r1,4) + 120*pow(r1,6)) + 256*pow(l0,8)*(pow(l1,8) - 48*pow(l1,6)*pow(r1,2) + 502*pow(l1,4)*pow(r1,4) + 105*pow(l1,2)*pow(r1,6) - 74*pow(r1,8)) + 1024*pow(l0,6)*pow(r1,2)*(11*pow(l1,8) - 218*pow(l1,6)*pow(r1,2) - 58*pow(l1,4)*pow(r1,4) + 63*pow(l1,2)*pow(r1,6) + 40*pow(r1,8)) - 4096*pow(l0,4)*pow(r1,2)*(pow(l1,6) - 26*pow(l1,4)*pow(r1,2) + 13*pow(l1,2)*pow(r1,4) - 12*pow(r1,6))*pow(pow(l1,2) - pow(r1,2),2) - 16384*pow(l0,2)*pow(r1,4)*(-6*pow(l1,4) + pow(l1,2)*pow(r1,2) + 8*pow(r1,4))*pow(-pow(l1,2) + pow(r1,2),3) + 65536*pow(r1,6)*pow(-pow(l1,2) + pow(r1,2),5)) + 16*pow(k,6)*pow(9*pow(l0,8) - 12*pow(l0,6)*(pow(l1,2) - 4*pow(r1,2)) + 16*pow(l0,4)*(pow(l1,4) - 17*pow(l1,2)*pow(r1,2) - 2*pow(r1,4)) + 64*pow(l0,2)*(5*pow(l1,4)*pow(r1,2) - pow(l1,2)*pow(r1,4) - 4*pow(r1,6)) + 256*pow(r1,2)*pow(-pow(l1,2) + pow(r1,2),3),2),-64*(pow(k,2)*pow(l0,7)*(27*pow(l0,12) + pow(l0,10)*(-132*pow(l1,2) + 216*pow(r1,2)) + 16*pow(l0,8)*(5*pow(l1,4) - 100*pow(l1,2)*pow(r1,2) + 9*pow(r1,4)) + 128*pow(l0,6)*(23*pow(l1,4)*pow(r1,2) - 19*pow(l1,2)*pow(r1,4) - 14*pow(r1,6)) - 1024*pow(l0,2)*pow(r1,4)*(12*pow(l1,6) + 14*pow(l1,4)*pow(r1,2) + pow(l1,2)*pow(r1,4) - 6*pow(r1,6)) - 256*pow(l0,4)*(2*pow(l1,6)*pow(r1,2) - 56*pow(l1,4)*pow(r1,4) - 28*pow(l1,2)*pow(r1,6) + 3*pow(r1,8)) - 4096*(2*pow(l1,6)*pow(r1,6) - 3*pow(l1,4)*pow(r1,8) + pow(r1,12))) + 4*pow(l0,9)*pow(l1,4)*pow(pow(l0,2) - 4*pow(r1,2),4) + 8*pow(k,4)*pow(l0,5)*(27*pow(l0,12) + pow(l0,10)*(-93*pow(l1,2) + 216*pow(r1,2)) + pow(l0,8)*(65*pow(l1,4) - 1484*pow(l1,2)*pow(r1,2) + 144*pow(r1,4)) - 16*pow(l0,6)*(3*pow(l1,6) - 173*pow(l1,4)*pow(r1,2) + 202*pow(l1,2)*pow(r1,4) + 112*pow(r1,6)) + 16*pow(l0,4)*(pow(l1,8) - 86*pow(l1,6)*pow(r1,2) + 894*pow(l1,4)*pow(r1,4) + 328*pow(l1,2)*pow(r1,6) - 48*pow(r1,8)) + 128*pow(l0,2)*pow(r1,2)*(7*pow(l1,8) - 120*pow(l1,6)*pow(r1,2) - 30*pow(l1,4)*pow(r1,4) + 14*pow(l1,2)*pow(r1,6) + 48*pow(r1,8)) - 256*pow(r1,4)*(-17*pow(l1,4) + 28*pow(l1,2)*pow(r1,2) + 16*pow(r1,4))*pow(pow(l1,2) - pow(r1,2),2)) + 16*pow(k,6)*pow(l0,3)*(27*pow(l0,12) - 54*pow(l0,10)*(pow(l1,2) - 4*pow(r1,2)) + 18*pow(l0,8)*(5*pow(l1,4) - 76*pow(l1,2)*pow(r1,2) + 8*pow(r1,4)) - 8*pow(l0,6)*(7*pow(l1,6) - 276*pow(l1,4)*pow(r1,2) + 504*pow(l1,2)*pow(r1,4) + 224*pow(r1,6)) + 32*pow(l0,4)*(pow(l1,8) - 81*pow(l1,6)*pow(r1,2) + 486*pow(l1,4)*pow(r1,4) + 104*pow(l1,2)*pow(r1,6) - 24*pow(r1,8)) + 384*pow(l0,2)*pow(r1,2)*(3*pow(l1,8) - 43*pow(l1,6)*pow(r1,2) + 12*pow(l1,4)*pow(r1,4) + 12*pow(l1,2)*pow(r1,6) + 16*pow(r1,8)) - 512*pow(r1,2)*(pow(l1,4) - 20*pow(l1,2)*pow(r1,2) - 8*pow(r1,4))*pow(pow(l1,2) - pow(r1,2),3))),64*(pow(k,2)*pow(l0,6)*(243*pow(l0,12) + pow(l0,10)*(-772*pow(l1,2) + 1368*pow(r1,2)) + 16*pow(l0,8)*(13*pow(l1,4) - 500*pow(l1,2)*pow(r1,2) - 39*pow(r1,4)) + 128*pow(l0,6)*(79*pow(l1,4)*pow(r1,2) - 51*pow(l1,2)*pow(r1,4) - 62*pow(r1,6)) - 1024*pow(l0,2)*pow(r1,4)*(12*pow(l1,6) + 14*pow(l1,4)*pow(r1,2) + pow(l1,2)*pow(r1,4) - 6*pow(r1,6)) - 256*pow(l0,4)*(2*pow(l1,6)*pow(r1,2) - 128*pow(l1,4)*pow(r1,4) - 76*pow(l1,2)*pow(r1,6) - 21*pow(r1,8)) - 4096*(2*pow(l1,6)*pow(r1,6) - 3*pow(l1,4)*pow(r1,8) + pow(r1,12))) + 4*pow(l0,8)*pow(l1,4)*pow(pow(l0,2) - 4*pow(r1,2),4) + 8*pow(k,4)*pow(l0,4)*(243*pow(l0,12) + pow(l0,10)*(-557*pow(l1,2) + 1368*pow(r1,2)) + pow(l0,8)*(289*pow(l1,4) - 7948*pow(l1,2)*pow(r1,2) - 624*pow(r1,4)) - 16*pow(l0,6)*(11*pow(l1,6) - 589*pow(l1,4)*pow(r1,2) + 634*pow(l1,2)*pow(r1,4) + 496*pow(r1,6)) + 128*pow(l0,2)*pow(r1,2)*(7*pow(l1,8) - 120*pow(l1,6)*pow(r1,2) - 30*pow(l1,4)*pow(r1,4) + 14*pow(l1,2)*pow(r1,6) + 48*pow(r1,8)) + 16*pow(l0,4)*(pow(l1,8) - 214*pow(l1,6)*pow(r1,2) + 2398*pow(l1,4)*pow(r1,4) + 1160*pow(l1,2)*pow(r1,6) + 336*pow(r1,8)) - 256*pow(r1,4)*(-17*pow(l1,4) + 28*pow(l1,2)*pow(r1,2) + 16*pow(r1,4))*pow(pow(l1,2) - pow(r1,2),2)) + 16*pow(k,6)*pow(l0,2)*(243*pow(l0,12) - 342*pow(l0,10)*(pow(l1,2) - 4*pow(r1,2)) + 6*pow(l0,8)*(79*pow(l1,4) - 1316*pow(l1,2)*pow(r1,2) - 104*pow(r1,4)) - 8*pow(l0,6)*(23*pow(l1,6) - 980*pow(l1,4)*pow(r1,2) + 1720*pow(l1,2)*pow(r1,4) + 992*pow(r1,6)) + 384*pow(l0,2)*pow(r1,2)*(3*pow(l1,8) - 43*pow(l1,6)*pow(r1,2) + 12*pow(l1,4)*pow(r1,4) + 12*pow(l1,2)*pow(r1,6) + 16*pow(r1,8)) + 96*pow(l0,4)*(pow(l1,8) - 75*pow(l1,6)*pow(r1,2) + 482*pow(l1,4)*pow(r1,4) + 184*pow(l1,2)*pow(r1,6) + 56*pow(r1,8)) - 512*pow(r1,2)*(pow(l1,4) - 20*pow(l1,2)*pow(r1,2) - 8*pow(r1,4))*pow(pow(l1,2) - pow(r1,2),3))),-1024*(8*pow(k,4)*pow(l0,7)*(75*pow(l0,8) + pow(l0,6)*(-106*pow(l1,2) + 272*pow(r1,2)) + 4*pow(l0,4)*(11*pow(l1,4) - 346*pow(l1,2)*pow(r1,2) - 88*pow(r1,4)) - 16*pow(l0,2)*(pow(l1,6) - 52*pow(l1,4)*pow(r1,2) + 54*pow(l1,2)*pow(r1,4) + 48*pow(r1,6)) - 64*(4*pow(l1,6)*pow(r1,2) - 47*pow(l1,4)*pow(r1,4) - 26*pow(l1,2)*pow(r1,6) - 12*pow(r1,8))) + pow(k,2)*pow(l0,9)*(75*pow(l0,8) - 16*pow(l0,6)*(9*pow(l1,2) - 17*pow(r1,2)) + 16*pow(l0,4)*(pow(l1,4) - 82*pow(l1,2)*pow(r1,2) - 22*pow(r1,4)) + 128*pow(l0,2)*(7*pow(l1,4)*pow(r1,2) - 4*pow(l1,2)*pow(r1,4) - 6*pow(r1,6)) + 768*(3*pow(l1,4)*pow(r1,4) + 2*pow(l1,2)*pow(r1,6) + pow(r1,8))) + 16*pow(k,6)*pow(l0,5)*(75*pow(l0,8) - 68*pow(l0,6)*(pow(l1,2) - 4*pow(r1,2)) + 16*pow(l0,4)*(5*pow(l1,4) - 91*pow(l1,2)*pow(r1,2) - 22*pow(r1,4)) - 16*pow(l0,2)*(pow(l1,6) - 44*pow(l1,4)*pow(r1,2) + 76*pow(l1,2)*pow(r1,4) + 48*pow(r1,6)) + 8*(pow(l1,8) - 72*pow(l1,6)*pow(r1,2) + 480*pow(l1,4)*pow(r1,4) + 224*pow(l1,2)*pow(r1,6) + 96*pow(r1,8)))),512*(8*pow(k,4)*pow(l0,6)*(443*pow(l0,8) + pow(l0,6)*(-346*pow(l1,2) + 912*pow(r1,2)) + 4*pow(l0,4)*(31*pow(l1,4) - 1066*pow(l1,2)*pow(r1,2) - 408*pow(r1,4)) - 16*pow(l0,2)*(pow(l1,6) - 52*pow(l1,4)*pow(r1,2) + 54*pow(l1,2)*pow(r1,4) + 48*pow(r1,6)) - 64*(4*pow(l1,6)*pow(r1,2) - 47*pow(l1,4)*pow(r1,4) - 26*pow(l1,2)*pow(r1,6) - 12*pow(r1,8))) + pow(k,2)*pow(l0,8)*(443*pow(l0,8) + pow(l0,6)*(-464*pow(l1,2) + 912*pow(r1,2)) + 16*pow(l0,4)*(pow(l1,4) - 242*pow(l1,2)*pow(r1,2) - 102*pow(r1,4)) + 128*pow(l0,2)*(7*pow(l1,4)*pow(r1,2) - 4*pow(l1,2)*pow(r1,4) - 6*pow(r1,6)) + 768*(3*pow(l1,4)*pow(r1,4) + 2*pow(l1,2)*pow(r1,6) + pow(r1,8))) + 16*pow(k,6)*pow(l0,4)*(443*pow(l0,8) - 228*pow(l0,6)*(pow(l1,2) - 4*pow(r1,2)) + 48*pow(l0,4)*(5*pow(l1,4) - 97*pow(l1,2)*pow(r1,2) - 34*pow(r1,4)) - 16*pow(l0,2)*(pow(l1,6) - 44*pow(l1,4)*pow(r1,2) + 76*pow(l1,2)*pow(r1,4) + 48*pow(r1,6)) + 8*(pow(l1,8) - 72*pow(l1,6)*pow(r1,2) + 480*pow(l1,4)*pow(r1,4) + 224*pow(l1,2)*pow(r1,6) + 96*pow(r1,8)))),-16384*pow(k,2)*(4*pow(k,2) + pow(l0,2))*pow(l0,7)*(25*pow(l0,6) - 12*pow(l0,4)*(pow(l1,2) - 2*pow(r1,2)) + 4*pow(k,2)*(25*pow(l0,4) - 6*pow(l0,2)*(pow(l1,2) - 4*pow(r1,2)) + 6*(pow(l1,4) - 20*pow(l1,2)*pow(r1,2) - 8*pow(r1,4))) - 48*pow(l0,2)*(2*pow(l1,2)*pow(r1,2) + pow(r1,4))),16384*pow(k,2)*(4*pow(k,2) + pow(l0,2))*pow(l0,6)*(27*pow(l0,6) - 4*pow(l0,4)*(pow(l1,2) - 2*pow(r1,2)) + 4*pow(k,2)*(27*pow(l0,4) - 2*pow(l0,2)*(pow(l1,2) - 4*pow(r1,2)) + 2*(pow(l1,4) - 20*pow(l1,2)*pow(r1,2) - 8*pow(r1,4))) - 16*pow(l0,2)*(2*pow(l1,2)*pow(r1,2) + pow(r1,4))),-262144*pow(k,2)*pow(l0,9)*pow(4*pow(k,2) + pow(l0,2),2),65536*pow(k,2)*pow(l0,8)*pow(4*pow(k,2) + pow(l0,2),2)};

  gsl_poly_complex_workspace * w
      = gsl_poly_complex_workspace_alloc (9);
  gsl_poly_complex_solve (a, 9, w, u);
  gsl_poly_complex_workspace_free (w);
  return 0;
  }

  // Solves for y for given values of x
int polysolvey(float l1 , float r1 ,float l0 ,float k, double x, double *y){
   /* coefficients of the univariate obtained */
   double a[5] = {k*(l0 - 2*x)*((l0 - x)*x + pow(l1,2) - pow(r1,2))*(-6*l0*x + pow(l0,2) - 2*pow(l1,2) + 2*pow(r1,2) + 6*pow(x,2)),-((l0 - 2*x)*(l0 + 2*l1 - 2*x)*(l0 - 2*(l1 + x))*pow(l0,2))/2.,-(k*(l0 - 2*x)*(-8*(l1 - r1)*(l1 + r1) + 3*pow(l0 - 2*x,2))),-2*(l0 - 2*x)*pow(l0,2),-6*k*(l0 - 2*x)};

   gsl_poly_complex_workspace * w
       = gsl_poly_complex_workspace_alloc (5);
   gsl_poly_complex_solve (a, 5, w, y);
   gsl_poly_complex_workspace_free (w);
   return 0;
   }

// Checks given x and y if it is a solution of eq3 also
bool singchecksol(float l1, float r1, float l0, float x, float y){
 double tol = pow(10, -8);
 float eq3 = pow(l0,4)*(pow(x,2) + pow(y,2)) - 2*x*pow(l0,3)*(-pow(l1,2) + pow(r1,2) + 3*(pow(x,2) + pow(y,2))) - 4*l0*x*(-pow(l1,2) + pow(r1,2) + pow(x,2) + pow(y,2))*(-pow(l1,2) + pow(r1,2) + 3*(pow(x,2) + pow(y,2))) + pow(l0,2)*(pow(l1,4) + pow(r1,4) + 10*pow(r1,2)*pow(x,2) + 13*pow(x,4) + 2*(pow(r1,2) + 9*pow(x,2))*pow(y,2) - 2*pow(l1,2)*(pow(r1,2) + 5*pow(x,2) + 3*pow(y,2)) + 5*pow(y,4)) + 4*(pow(x,2) + pow(y,2))*pow(-pow(l1,2) + pow(r1,2) + pow(x,2) + pow(y,2),2);
 if(eq3<tol)
   return true;
 else
   return false;
}

// The feasible radius values because of the bigger factor
MatrixXd radgain(float l1, float r1, float l0, float k){
  double tol = pow(10, -7);
  float pen = std::numeric_limits<float>::max();
  int count=0, i=0, j=0;
  double xval[16], yval[8];
  // Since there can be at max 8 x-roots and 4 y-roots, creating a container of size 32
  MatrixXd radval = MatrixXd::Constant(32, 1, pen);
  polysolvex(l1, r1, l0, k, xval);
  std::cout<<"ok polysolvex"<<std::endl;
  std::cout<<xval[3]<<std::endl;
  for(int i=0; i<8;i++){
    std::cout<<"in for loop"<<std::endl;
    if(xval[2*i+1]<tol){
      std::cout<<"solving for y"<<std::endl;
      polysolvey(l1, r1, l0, k, xval[2*i], yval);
      std::cout<<"ok polysolvy"<<std::endl;
      for(int j=0;j<4;j++){
        if(fabs(yval[2*j+1])<tol){
          if(singchecksol(l1, r1, l0, xval[2*i], yval[2*j])){
            if(doublerootcheck(l1, r1, l0, xval[2*i], yval[2*j])>tol){
              radval(count) = pow(pow(xval[2*i]-(l0/2), 2)+pow(yval[2*j]-k, 2), 0.5);
              count++;
              }
            }
          }
        }
      }
    }
  return radval;
}

// Gives the minimum possible feasible radius
float minrad(float l1, float r1, float l0, float k){
  MatrixXd workrad(2,1);
  MatrixXd fsmrad(25, 1);
  MatrixXd gainrad(32, 1);
  workrad = workspacerad(l1, r1, l0, k);
  // std::cout<<"ok workrad"<<std::endl;
  fsmrad = radfsm(l1, r1, l0, k);
  // std::cout<<"ok radfsm"<<std::endl;
  gainrad = radgain(l1, r1, l0, k);
  // std::cout<<"ok radgain"<<std::endl;
  // std::cout << fsmrad << '\n';
  // std::cout << gainrad << '\n';
  MatrixXd allrads(workrad.rows()+fsmrad.rows()+gainrad.rows(), workrad.cols());
  allrads << workrad, fsmrad, gainrad;
  // std::cout << workrad.minCoeff() << '\n';
  // std::cout << allrads.minCoeff() << '\n';
  return allrads.minCoeff();
  return 0;
}

//Now that we know the min radius we can compute the required objective
MatrixXd objective(float l1, float r1, float l0, float k, float rad){
  //Now we know the center (0, k) and the radius of the circle i.e. rad, so discretise the workspace
  float pen = std::numeric_limits<float>::max();
  double g = 0.0, h = 0.0, minr = pen, disc = 0.0, momval = 0.0, ikappa = 0.0, area = 0.0, th1 = 0.0, th2 = 0.0;
  disc = pow(10, 1);
  MatrixXd ikval(2, 1);
  minr = minrad(l1, r1, l0, k);
  area = M_PI*minr*minr;
  for(float alpha=0;alpha<2*M_PI;alpha+=pow(disc, -1)){
    for(float radius=0;radius<minr;radius+=pow(disc, -1)){
      g = (l0/2)+radius*cos(alpha);
      h = k+radius*sin(alpha);
      ikval = ikcal(l1, r1, l0, g, h);
      th1 = ikval(0);
      th2 = ikval(1);
      momval+=momcal(l1, r1, l0, g, h, th1, th2);
      ikappa += kinvcal(l1, r1, l0, g, h, th1, th2);
    }
  }
MatrixXd obj(3, 1);
obj(0) = momval/area;
obj(1) = ikappa/area;
obj(2) = area;
return obj;
}

#include <iostream>
#include <stdio.h>
#include <Eigen/Dense>
#include <cmath>
#include <gsl/gsl_poly.h>
#include <limits>
#include "design.cpp"
using namespace Eigen;

bool inWorkspace(float l1, float r1, float l0, float k);
MatrixXd ikcal(float l, float r, float l0, float x, float y);
double momcal(float l, float r, float l0, float x, float y, float th1, float th2);
double kinvcal(float l, float r, float l0, float x, float y, float th1, float th2);
MatrixXd radfsm(float l1, float r1, float l0, float k);
MatrixXd workspacerad(float l1, float r1, float l0, float k);
float doublerootcheck(float l1, float r1, float l0, float p, float q);
double* polysolvex(float l1 , float r1 ,float l0 ,float k, double *u);
double* polysolvey(float l1 , float r1 ,float l0 ,float k, double x, double *y);
bool singchecksol(float l1, float r1, float l0, float x, float y);
MatrixXd radgain(float l1, float r1, float l0, float k);
float minrad(float l1, float r1, float l0, float k);
MatrixXd objective(float l1, float r1, float l0, float k, float rad);

int main()
{
  //This means that we need to minimise MoM and GCI and max Area
  double l1 = 1.771548, r1 = 0.851070, l0 = 0.652008, k = 1.887573, x = 0.8, y = 0.8;
  MatrixXd objs(3, 1);
  float minr = minrad(l1, r1, l0, k);
  std::cout<<minr<<std::endl;
  objs = objective(l1, r1, l0, k, minr);
  std::cout<<objs<<std::endl;
  return 0;
}

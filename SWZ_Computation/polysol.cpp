// polysol.cpp

# include <stdio.h>
// # include <iostream>
# include <gsl/gsl_poly.h>
# include <math.h>
# include <Eigen/Dense>
# include <limits>
// #include "matfun.cpp"

using namespace std;
using namespace Eigen;

// handling all imaginary roots verified on 09-07-2018
//
// int main(int argc, char const *argv[]) {
//   float l1 = 0.981747, r1 = 2.130107, l0 = 0.574475, k = 1.513348, rads = 0, finrad = 0;
//   bool cond;
//
//   cond = inWorkspace(l1, r1, l0, k);
//   cout<<"inWorkspace \n"<<cond<<endl;
//
//   MatrixXd sol1(3, 1);
//   sol1 = returnsol(l1, r1, l0, k);
//   cout<< "The sol1 is \n" <<sol1<<endl;
//   MatrixXd sol2(2, 15);
//   sol2 = returnsol2(l1, r1, l0, k);
//   cout<<"The sol2 is \n"<<sol2<<endl;
//   float sol3;
//   sol3 = fsmrad(l1, r1, l0, k);
//   cout<<"The sol3 is \n"<<sol3<<endl;
//   float sol4;
//   sol4 = workspacerad(l1, r1, l0, k);
//   cout<<"The sol4 is \n"<<sol4<<endl;
//
//   finrad = selectrad(sol1, sol2, sol3, sol4);
//
//   cout<<"The finaf rad is \n"<<finrad<<endl;
//   cout<< "getArea"<<getA(finrad)<<endl;
//
//   // VectorXd sol(3);
//   // sol = computeobj(l1, r1, l0, finrad, k);
//   // cout<<sol<<endl;
//
//   return 0;
// };


//Selects the smallest radius among all the other solutions
float selectrad(MatrixXd sol1, MatrixXd sol2, float sol3, float sol4){
  int i = 0, j = 0;
  MatrixXd finrads(4, 1);
  MatrixXd temp(1, 1);
  finrads << sol1(0,0), sol2(0,0), sol3, sol4;
  for(i=0;i<4;i++)
    for(j=i+1;j<4;j++){
      if(finrads(i, 0) > finrads(j, 0)){
        temp = finrads.row(i);
        finrads.row(i) = finrads.row(j);
        finrads.row(j) = temp;
      }
    }
    return finrads(0,0);
}

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

//We could do this, i.e. using a function to compute the area, but for a lot of runs this is inefficient as it adds to the number of jumps in the code. So compute area of the SWZ when and where needed
//returns the area of the SWZ
//float getA(float rad){
//  float A = M_PI*rad*rad;
//  return A;
//}


//Function to return the radius with the fsm circles
//Rewrite fsm1 again
float fsmrad(float l1, float r1, float l0, float k){
  float c = pow(pow(l1, 2)-pow(l0/2, 2), 0.5);
  float val = numeric_limits<float>::max();
  float rad1 = val, rad2 = val, rad3 = val, rad4 = val;
  if(abs(k-c)>=r1)
    rad1 = abs(k-c)-r1;
  else if(abs(k-c)<r1)
    rad2 = r1-abs(k-c);
  if(k+c>=r1)
    rad3 = k+c-r1;
  if(k+c<r1)
    rad4 = r1-(k+c);

    int i=0, j=0;
    MatrixXd min(4, 1);
    MatrixXd temp(1, 1);
    min << rad1, rad2, rad3, rad4;
    for(i=0;i<4;i++)
      for(j=i+1;j<4;j++){
        if(min(i, 0) > min(j, 0)){
          temp = min.row(i);
          min.row(i) = min.row(j);
          min.row(j) = temp;
        }
      }
      return min(0,0);
}


//Gets the radius with the tangency problem between the workspace and the circle chosen
float workspacerad(float l1, float r1, float l0, float k){
  float d1 = pow(pow(l0/2,2) + pow(k,2),0.5), rad1 = 0, rad2 = 0;
  rad1 = l1+r1-d1;
  rad2 = d1 - abs(l1-r1);
  if(rad1<=rad2){
    return rad1;
  }
  if(rad2<rad1){
    return rad2;
  }
}


//Gives the radius for the factored out parts of the equation
MatrixXd returnsol2(float l1, float r1 ,float l0 ,float k){
  int i=0, j=0, count = 0;
      MatrixXd a(1, 30);
      a.setZero();
      a << l0/2.,0,l0/2.,-pow(-4*l0*r1 - pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5)/2.,l0/2.,pow(-4*l0*r1 - pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5)/2.,l0/2.,-pow(4*l0*r1 - pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5)/2.,l0/2.,pow(4*l0*r1 - pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5)/2.,(l0 - pow(pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5))/2.,0,(l0 - pow(pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5))/2.,-(pow(2,-0.5)*pow(-5*pow(l0,2) - 4*pow(l1,2) + 4*pow(r1,2) - pow(9*pow(l0,4) + 8*pow(l0,2)*(5*pow(l1,2) + 3*pow(r1,2)) + 16*pow(pow(l1,2) - pow(r1,2),2),0.5),0.5))/2.,(l0 - pow(pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5))/2.,(pow(2,-0.5)*pow(-5*pow(l0,2) - 4*pow(l1,2) + 4*pow(r1,2) - pow(9*pow(l0,4) + 8*pow(l0,2)*(5*pow(l1,2) + 3*pow(r1,2)) + 16*pow(pow(l1,2) - pow(r1,2),2),0.5),0.5))/2.,(l0 - pow(pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5))/2.,-(pow(2,-0.5)*pow(-5*pow(l0,2) - 4*pow(l1,2) + 4*pow(r1,2) + pow(9*pow(l0,4) + 8*pow(l0,2)*(5*pow(l1,2) + 3*pow(r1,2)) + 16*pow(pow(l1,2) - pow(r1,2),2),0.5),0.5))/2.,(l0 - pow(pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5))/2.,pow((-5*pow(l0,2))/8. + (-4*pow(l1,2) + 4*pow(r1,2) + pow(9*pow(l0,4) + 8*pow(l0,2)*(5*pow(l1,2) + 3*pow(r1,2)) + 16*pow(pow(l1,2) - pow(r1,2),2),0.5))/8.,0.5),(l0 + pow(pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5))/2.,0,(l0 + pow(pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5))/2.,-(pow(2,-0.5)*pow(-5*pow(l0,2) - 4*pow(l1,2) + 4*pow(r1,2) - pow(9*pow(l0,4) + 8*pow(l0,2)*(5*pow(l1,2) + 3*pow(r1,2)) + 16*pow(pow(l1,2) - pow(r1,2),2),0.5),0.5))/2.,(l0 + pow(pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5))/2.,(pow(2,-0.5)*pow(-5*pow(l0,2) - 4*pow(l1,2) + 4*pow(r1,2) - pow(9*pow(l0,4) + 8*pow(l0,2)*(5*pow(l1,2) + 3*pow(r1,2)) + 16*pow(pow(l1,2) - pow(r1,2),2),0.5),0.5))/2.,(l0 + pow(pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5))/2.,-(pow(2,-0.5)*pow(-5*pow(l0,2) - 4*pow(l1,2) + 4*pow(r1,2) + pow(9*pow(l0,4) + 8*pow(l0,2)*(5*pow(l1,2) + 3*pow(r1,2)) + 16*pow(pow(l1,2) - pow(r1,2),2),0.5),0.5))/2.,(l0 + pow(pow(l0,2) + 4*pow(l1,2) - 4*pow(r1,2),0.5))/2.,pow((-5*pow(l0,2))/8. + (-4*pow(l1,2) + 4*pow(r1,2) + pow(9*pow(l0,4) + 8*pow(l0,2)*(5*pow(l1,2) + 3*pow(r1,2)) + 16*pow(pow(l1,2) - pow(r1,2),2),0.5))/8.,0.5);

      //cout<<a <<endl;
      Map<MatrixXd> b(a.data(), 2, 15);
      // cout<<"The reshaped matrix is \n"<<b<<endl;

      MatrixXd temp3(15, 3);
      temp3.setZero(15, 3);
      for(i=0;i<15;i++)
      //This statement removes all imaginary values
        if (b(0, i)==b(0, i) && b(1, i)==b(1, i)){
          temp3.block<1, 2>(count, 1) = b.col(i).transpose();
          temp3(count, 0) = pow(pow((b(0, i)-(l0/2)), 2) + pow((b(1, i)-k), 2), 0.5);
          count++;
        }

      //This returns the penality if there are no real solutions
      if(count==0){
        MatrixXd penality(3, 1);
        penality.setOnes();
        float val = numeric_limits<float>::max();
        penality = val * penality;
        return penality;
      }

      MatrixXd sols(count, 3);
      sols.setZero(count, 3);
      for(i=0;i<count;i++){
        sols.row(i) = temp3.row(i);
      }

      MatrixXd temp4(1, 3);
      temp4.setZero(1, 3);
      for(i=0;i<count;i++)
        for(j=i+1;j<count;j++){
          if(sols(i, 0) > sols(j, 0)){
            temp4 = sols.row(i);
            sols.row(i) = sols.row(j);
            sols.row(j) = temp4;
          }
        }

//Handling double point
  i = 0;
  float doubletangentcond = 1;
  while(doubletangentcond>pow(10, -4)&&i<count){
    //Condition for eliminating acnodes from the solutions
    float p = sols(i, 1), q = sols(i, 2);

    doubletangentcond = -24*p*pow(l0,7) + pow(l0,8) + 232*p*pow(l0,5)*pow(l1,2) - 16*pow(l0,6)*pow(l1,2) - 544*p*pow(l0,3)*pow(l1,4) + 68*pow(l0,4)*pow(l1,4) + 256*l0*p*pow(l1,6) - 64*pow(l0,2)*pow(l1,6) + 16*pow(l1,8) + 204*pow(l0,6)*pow(p,2) - 1288*pow(l0,4)*pow(l1,2)*pow(p,2) + 1600*pow(l0,2)*pow(l1,4)*pow(p,2) - 256*pow(l1,6)*pow(p,2) - 936*pow(l0,5)*pow(p,3) + 3648*pow(l0,3)*pow(l1,2)*pow(p,3) - 2112*l0*pow(l1,4)*pow(p,3) + 2628*pow(l0,4)*pow(p,4) - 5664*pow(l0,2)*pow(l1,2)*pow(p,4) + 1056*pow(l1,4)*pow(p,4) - 4608*pow(l0,3)*pow(p,5) + 4608*l0*pow(l1,2)*pow(p,5) + 4896*pow(l0,2)*pow(p,6) - 1536*pow(l1,2)*pow(p,6) - 2880*l0*pow(p,7) + 720*pow(p,8) + 936*p*pow(l0,5)*pow(q,2) - 96*pow(l0,6)*pow(q,2) - 2304*p*pow(l0,3)*pow(l1,2)*pow(q,2) + 296*pow(l0,4)*pow(l1,2)*pow(q,2) + 960*l0*p*pow(l1,4)*pow(q,2) - 256*pow(l0,2)*pow(l1,4)*pow(q,2) - 256*pow(l1,6)*pow(q,2) - 4104*pow(l0,4)*pow(p,2)*pow(q,2) + 6912*pow(l0,2)*pow(l1,2)*pow(p,2)*pow(q,2) - 960*pow(l1,4)*pow(p,2)*pow(q,2) + 10368*pow(l0,3)*pow(p,3)*pow(q,2) - 9216*l0*pow(l1,2)*pow(p,3)*pow(q,2) - 15264*pow(l0,2)*pow(p,4)*pow(q,2) + 4608*pow(l1,2)*pow(p,4)*pow(q,2) + 12096*l0*pow(p,5)*pow(q,2) - 4032*pow(p,6)*pow(q,2) + 4608*p*pow(l0,3)*pow(q,4) - 540*pow(l0,4)*pow(q,4) - 4608*l0*p*pow(l1,2)*pow(q,4) + 1056*pow(l0,2)*pow(l1,2)*pow(q,4) + 1056*pow(l1,4)*pow(q,4) - 14112*pow(l0,2)*pow(p,2)*pow(q,4) + 4608*pow(l1,2)*pow(p,2)*pow(q,4) + 19008*l0*pow(p,3)*pow(q,4) - 9504*pow(p,4)*pow(q,4) + 4032*l0*p*pow(q,6) - 864*pow(l0,2)*pow(q,6) - 1536*pow(l1,2)*pow(q,6) - 4032*pow(p,2)*pow(q,6) + 720*pow(q,8) - 160*p*pow(l0,5)*pow(r1,2) + 12*pow(l0,6)*pow(r1,2) + 896*p*pow(l0,3)*pow(l1,2)*pow(r1,2) - 96*pow(l0,4)*pow(l1,2)*pow(r1,2) - 768*l0*p*pow(l1,4)*pow(r1,2) + 176*pow(l0,2)*pow(l1,4)*pow(r1,2) - 64*pow(l1,6)*pow(r1,2) + 976*pow(l0,4)*pow(p,2)*pow(r1,2) - 3008*pow(l0,2)*pow(l1,2)*pow(p,2)*pow(r1,2) + 768*pow(l1,4)*pow(p,2)*pow(r1,2) - 3168*pow(l0,3)*pow(p,3)*pow(r1,2) + 4224*l0*pow(l1,2)*pow(p,3)*pow(r1,2) + 5424*pow(l0,2)*pow(p,4)*pow(r1,2) - 2112*pow(l1,2)*pow(p,4)*pow(r1,2) - 4608*l0*pow(p,5)*pow(r1,2) + 1536*pow(p,6)*pow(r1,2) + 2592*p*pow(l0,3)*pow(q,2)*pow(r1,2) - 368*pow(l0,4)*pow(q,2)*pow(r1,2) - 1920*l0*p*pow(l1,2)*pow(q,2)*pow(r1,2) + 576*pow(l0,2)*pow(l1,2)*pow(q,2)*pow(r1,2) + 768*pow(l1,4)*pow(q,2)*pow(r1,2) - 7200*pow(l0,2)*pow(p,2)*pow(q,2)*pow(r1,2) + 1920*pow(l1,2)*pow(p,2)*pow(q,2)*pow(r1,2) + 9216*l0*pow(p,3)*pow(q,2)*pow(r1,2) - 4608*pow(p,4)*pow(q,2)*pow(r1,2) + 4608*l0*p*pow(q,4)*pow(r1,2) - 1104*pow(l0,2)*pow(q,4)*pow(r1,2) - 2112*pow(l1,2)*pow(q,4)*pow(r1,2) - 4608*pow(p,2)*pow(q,4)*pow(r1,2) + 1536*pow(q,6)*pow(r1,2) - 352*p*pow(l0,3)*pow(r1,4) + 28*pow(l0,4)*pow(r1,4) + 768*l0*p*pow(l1,2)*pow(r1,4) - 160*pow(l0,2)*pow(l1,2)*pow(r1,4) + 96*pow(l1,4)*pow(r1,4) + 1408*pow(l0,2)*pow(p,2)*pow(r1,4) - 768*pow(l1,2)*pow(p,2)*pow(r1,4) - 2112*l0*pow(p,3)*pow(r1,4) + 1056*pow(p,4)*pow(r1,4) + 960*l0*p*pow(q,2)*pow(r1,4) - 320*pow(l0,2)*pow(q,2)*pow(r1,4) - 768*pow(l1,2)*pow(q,2)*pow(r1,4) - 960*pow(p,2)*pow(q,2)*pow(r1,4) + 1056*pow(q,4)*pow(r1,4) - 256*l0*p*pow(r1,6) + 48*pow(l0,2)*pow(r1,6) - 64*pow(l1,2)*pow(r1,6) + 256*pow(p,2)*pow(r1,6) + 256*pow(q,2)*pow(r1,6) + 16*pow(r1,8);

    i++;
  }

  //This returns the penality if there are no real solutions
  if(i==count){
    MatrixXd penality(3, 1);
    penality.setOnes();
    float val = numeric_limits<float>::max();
    penality = val * penality;
    return penality;
  }

      VectorXd sol2(3, 1);
      sol2.setZero(3, 1);
      sol2 = sols.row(i-1);
      return sol2.transpose();
}

//Uses the polynimial solution to find all the radii and choses the smallest solution which is not an acnode
MatrixXd returnsol(float l1, float r1 ,float l0 ,float k){
  int i=0, j=0;
  double u[16];
  polysolvex(l1, r1, l0, k, u);

  double n[12];
  int count = 0, count2 = 0;
  MatrixXd temp(50, 3);
  temp.setZero(50, 3);
  for(i=0;i<8;i++){
  if(abs(u[2*i+1])<=pow(10, -4)){
    polysolvey(l1, r1, l0, k, n, u[2*i]);
    for(j=0;j<6;j++){
      if(abs(n[2*j+1])<=pow(10, -4)){
        temp(count, 0) = pow(pow(u[2*i]-(0.5*l0), 2)+pow(n[2*j]-k, 2), 0.5);
        temp(count, 1) = u[2*i];
        temp(count, 2) = n[2*j];
        count++;
      }
    }
    count2++;
  }
}


if (count==0 || count2==0){
    MatrixXd penality(3,1);
    penality.setOnes();
    float val = numeric_limits<float>::max();
    penality = val * penality;
    return penality;
  }

  MatrixXd rads(count, 3);
  rads.setZero(count, 3);
  for(i=0; i<count; i++)
    for(j=0; j<3; j++)
      rads(i, j) = temp(i, j);


  MatrixXd temp2(1, 3);
  for(i=0;i<count;i++)
    for(j=i+1;j<count;j++){
      if(rads(i, 0) > rads(j, 0)){
        temp2 = rads.row(i);
        rads.row(i) = rads.row(j);
        rads.row(j) = temp2;
      }
    }

  i = 0;
  float doubletangentcond = 1;
  while(doubletangentcond>pow(10, -4)&&i<count){
    //Condition for eliminating acnodes from the solutions

    float p = rads(i, 1), q = rads(i, 2);

    doubletangentcond = -24*p*pow(l0,7) + pow(l0,8) + 232*p*pow(l0,5)*pow(l1,2) - 16*pow(l0,6)*pow(l1,2) - 544*p*pow(l0,3)*pow(l1,4) + 68*pow(l0,4)*pow(l1,4) + 256*l0*p*pow(l1,6) - 64*pow(l0,2)*pow(l1,6) + 16*pow(l1,8) + 204*pow(l0,6)*pow(p,2) - 1288*pow(l0,4)*pow(l1,2)*pow(p,2) + 1600*pow(l0,2)*pow(l1,4)*pow(p,2) - 256*pow(l1,6)*pow(p,2) - 936*pow(l0,5)*pow(p,3) + 3648*pow(l0,3)*pow(l1,2)*pow(p,3) - 2112*l0*pow(l1,4)*pow(p,3) + 2628*pow(l0,4)*pow(p,4) - 5664*pow(l0,2)*pow(l1,2)*pow(p,4) + 1056*pow(l1,4)*pow(p,4) - 4608*pow(l0,3)*pow(p,5) + 4608*l0*pow(l1,2)*pow(p,5) + 4896*pow(l0,2)*pow(p,6) - 1536*pow(l1,2)*pow(p,6) - 2880*l0*pow(p,7) + 720*pow(p,8) + 936*p*pow(l0,5)*pow(q,2) - 96*pow(l0,6)*pow(q,2) - 2304*p*pow(l0,3)*pow(l1,2)*pow(q,2) + 296*pow(l0,4)*pow(l1,2)*pow(q,2) + 960*l0*p*pow(l1,4)*pow(q,2) - 256*pow(l0,2)*pow(l1,4)*pow(q,2) - 256*pow(l1,6)*pow(q,2) - 4104*pow(l0,4)*pow(p,2)*pow(q,2) + 6912*pow(l0,2)*pow(l1,2)*pow(p,2)*pow(q,2) - 960*pow(l1,4)*pow(p,2)*pow(q,2) + 10368*pow(l0,3)*pow(p,3)*pow(q,2) - 9216*l0*pow(l1,2)*pow(p,3)*pow(q,2) - 15264*pow(l0,2)*pow(p,4)*pow(q,2) + 4608*pow(l1,2)*pow(p,4)*pow(q,2) + 12096*l0*pow(p,5)*pow(q,2) - 4032*pow(p,6)*pow(q,2) + 4608*p*pow(l0,3)*pow(q,4) - 540*pow(l0,4)*pow(q,4) - 4608*l0*p*pow(l1,2)*pow(q,4) + 1056*pow(l0,2)*pow(l1,2)*pow(q,4) + 1056*pow(l1,4)*pow(q,4) - 14112*pow(l0,2)*pow(p,2)*pow(q,4) + 4608*pow(l1,2)*pow(p,2)*pow(q,4) + 19008*l0*pow(p,3)*pow(q,4) - 9504*pow(p,4)*pow(q,4) + 4032*l0*p*pow(q,6) - 864*pow(l0,2)*pow(q,6) - 1536*pow(l1,2)*pow(q,6) - 4032*pow(p,2)*pow(q,6) + 720*pow(q,8) - 160*p*pow(l0,5)*pow(r1,2) + 12*pow(l0,6)*pow(r1,2) + 896*p*pow(l0,3)*pow(l1,2)*pow(r1,2) - 96*pow(l0,4)*pow(l1,2)*pow(r1,2) - 768*l0*p*pow(l1,4)*pow(r1,2) + 176*pow(l0,2)*pow(l1,4)*pow(r1,2) - 64*pow(l1,6)*pow(r1,2) + 976*pow(l0,4)*pow(p,2)*pow(r1,2) - 3008*pow(l0,2)*pow(l1,2)*pow(p,2)*pow(r1,2) + 768*pow(l1,4)*pow(p,2)*pow(r1,2) - 3168*pow(l0,3)*pow(p,3)*pow(r1,2) + 4224*l0*pow(l1,2)*pow(p,3)*pow(r1,2) + 5424*pow(l0,2)*pow(p,4)*pow(r1,2) - 2112*pow(l1,2)*pow(p,4)*pow(r1,2) - 4608*l0*pow(p,5)*pow(r1,2) + 1536*pow(p,6)*pow(r1,2) + 2592*p*pow(l0,3)*pow(q,2)*pow(r1,2) - 368*pow(l0,4)*pow(q,2)*pow(r1,2) - 1920*l0*p*pow(l1,2)*pow(q,2)*pow(r1,2) + 576*pow(l0,2)*pow(l1,2)*pow(q,2)*pow(r1,2) + 768*pow(l1,4)*pow(q,2)*pow(r1,2) - 7200*pow(l0,2)*pow(p,2)*pow(q,2)*pow(r1,2) + 1920*pow(l1,2)*pow(p,2)*pow(q,2)*pow(r1,2) + 9216*l0*pow(p,3)*pow(q,2)*pow(r1,2) - 4608*pow(p,4)*pow(q,2)*pow(r1,2) + 4608*l0*p*pow(q,4)*pow(r1,2) - 1104*pow(l0,2)*pow(q,4)*pow(r1,2) - 2112*pow(l1,2)*pow(q,4)*pow(r1,2) - 4608*pow(p,2)*pow(q,4)*pow(r1,2) + 1536*pow(q,6)*pow(r1,2) - 352*p*pow(l0,3)*pow(r1,4) + 28*pow(l0,4)*pow(r1,4) + 768*l0*p*pow(l1,2)*pow(r1,4) - 160*pow(l0,2)*pow(l1,2)*pow(r1,4) + 96*pow(l1,4)*pow(r1,4) + 1408*pow(l0,2)*pow(p,2)*pow(r1,4) - 768*pow(l1,2)*pow(p,2)*pow(r1,4) - 2112*l0*pow(p,3)*pow(r1,4) + 1056*pow(p,4)*pow(r1,4) + 960*l0*p*pow(q,2)*pow(r1,4) - 320*pow(l0,2)*pow(q,2)*pow(r1,4) - 768*pow(l1,2)*pow(q,2)*pow(r1,4) - 960*pow(p,2)*pow(q,2)*pow(r1,4) + 1056*pow(q,4)*pow(r1,4) - 256*l0*p*pow(r1,6) + 48*pow(l0,2)*pow(r1,6) - 64*pow(l1,2)*pow(r1,6) + 256*pow(p,2)*pow(r1,6) + 256*pow(q,2)*pow(r1,6) + 16*pow(r1,8);
    i++;
  }

  if (i==count){
      MatrixXd penality(3,1);
      penality.setOnes();
      float val = numeric_limits<float>::max();
      penality = val * penality;
      return penality;
    }


  MatrixXd sol(1,3);
  sol = rads.row(i-1);

  return sol.transpose();

}


//Solves the polynomial in y by substituting the x values into the gain type curve
int polysolvey(float l1 , float r1 ,float l0 ,float k, double *z, float x){
  /* coefficients of the univariate obtained */

  double a[7] = {2*x*pow(l0,3)*pow(l1,2) - 4*l0*x*pow(l1,4) + pow(l0,2)*pow(l1,4) - 2*x*pow(l0,3)*pow(r1,2) + 8*l0*x*pow(l1,2)*pow(r1,2) - 2*pow(l0,2)*pow(l1,2)*pow(r1,2) - 4*l0*x*pow(r1,4) + pow(l0,2)*pow(r1,4) + pow(l0,4)*pow(x,2) - 10*pow(l0,2)*pow(l1,2)*pow(x,2) + 4*pow(l1,4)*pow(x,2) + 10*pow(l0,2)*pow(r1,2)*pow(x,2) - 8*pow(l1,2)*pow(r1,2)*pow(x,2) + 4*pow(r1,4)*pow(x,2) - 6*pow(l0,3)*pow(x,3) + 16*l0*pow(l1,2)*pow(x,3) - 16*l0*pow(r1,2)*pow(x,3) + 13*pow(l0,2)*pow(x,4) - 8*pow(l1,2)*pow(x,4) + 8*pow(r1,2)*pow(x,4) - 12*l0*pow(x,5) + 4*pow(x,6),0,-6*x*pow(l0,3) + pow(l0,4) + 16*l0*x*pow(l1,2) - 6*pow(l0,2)*pow(l1,2) + 4*pow(l1,4) - 16*l0*x*pow(r1,2) + 2*pow(l0,2)*pow(r1,2) - 8*pow(l1,2)*pow(r1,2) + 4*pow(r1,4) + 18*pow(l0,2)*pow(x,2) - 16*pow(l1,2)*pow(x,2) + 16*pow(r1,2)*pow(x,2) - 24*l0*pow(x,3) + 12*pow(x,4),0,-12*l0*x + 5*pow(l0,2) - 8*pow(l1,2) + 8*pow(r1,2) + 12*pow(x,2),0,4};

  gsl_poly_complex_workspace * w
      = gsl_poly_complex_workspace_alloc (7);

  gsl_poly_complex_solve (a, 7, w, z);

  gsl_poly_complex_workspace_free (w);

return 0;
}


//Solves the univariate in x
int polysolvex(float l1 , float r1 ,float l0 ,float k, double *u){
  /* coefficients of the univariate obtained */

  double a[9] = {1296*pow(k,6)*pow(l0,16) + 648*pow(k,4)*pow(l0,18) + 81*pow(k,2)*pow(l0,20) - 3456*pow(k,6)*pow(l0,14)*pow(l1,2) - 3168*pow(k,4)*pow(l0,16)*pow(l1,2) - 576*pow(k,2)*pow(l0,18)*pow(l1,2) + 6912*pow(k,6)*pow(l0,12)*pow(l1,4) + 3072*pow(k,4)*pow(l0,14)*pow(l1,4) + 736*pow(k,2)*pow(l0,16)*pow(l1,4) + 64*pow(l0,18)*pow(l1,4) - 6144*pow(k,6)*pow(l0,10)*pow(l1,6) - 2560*pow(k,4)*pow(l0,12)*pow(l1,6) + 4096*pow(k,6)*pow(l0,8)*pow(l1,8) + 2048*pow(k,4)*pow(l0,10)*pow(l1,8) + 256*pow(k,2)*pow(l0,12)*pow(l1,8) + 13824*pow(k,6)*pow(l0,14)*pow(r1,2) + 6912*pow(k,4)*pow(l0,16)*pow(r1,2) + 864*pow(k,2)*pow(l0,18)*pow(r1,2) - 96768*pow(k,6)*pow(l0,12)*pow(l1,2)*pow(r1,2) - 57216*pow(k,4)*pow(l0,14)*pow(l1,2)*pow(r1,2) - 8256*pow(k,2)*pow(l0,16)*pow(l1,2)*pow(r1,2) + 221184*pow(k,6)*pow(l0,10)*pow(l1,4)*pow(r1,2) + 147968*pow(k,4)*pow(l0,12)*pow(l1,4)*pow(r1,2) + 18688*pow(k,2)*pow(l0,14)*pow(l1,4)*pow(r1,2) - 1024*pow(l0,16)*pow(l1,4)*pow(r1,2) - 335872*pow(k,6)*pow(l0,8)*pow(l1,6)*pow(r1,2) - 98304*pow(k,4)*pow(l0,10)*pow(l1,6)*pow(r1,2) - 9216*pow(k,2)*pow(l0,12)*pow(l1,6)*pow(r1,2) + 262144*pow(k,6)*pow(l0,6)*pow(l1,8)*pow(r1,2) + 90112*pow(k,4)*pow(l0,8)*pow(l1,8)*pow(r1,2) - 4096*pow(k,2)*pow(l0,10)*pow(l1,8)*pow(r1,2) - 131072*pow(k,6)*pow(l0,4)*pow(l1,10)*pow(r1,2) - 32768*pow(k,4)*pow(l0,6)*pow(l1,10)*pow(r1,2) + 27648*pow(k,6)*pow(l0,12)*pow(r1,4) + 13824*pow(k,4)*pow(l0,14)*pow(r1,4) + 1728*pow(k,2)*pow(l0,16)*pow(r1,4) - 423936*pow(k,6)*pow(l0,10)*pow(l1,2)*pow(r1,4) - 189952*pow(k,4)*pow(l0,12)*pow(l1,2)*pow(r1,4) - 20992*pow(k,2)*pow(l0,14)*pow(l1,2)*pow(r1,4) + 1904640*pow(k,6)*pow(l0,8)*pow(l1,4)*pow(r1,4) + 1028096*pow(k,4)*pow(l0,10)*pow(l1,4)*pow(r1,4) + 158208*pow(k,2)*pow(l0,12)*pow(l1,4)*pow(r1,4) + 6144*pow(l0,14)*pow(l1,4)*pow(r1,4) - 3506176*pow(k,6)*pow(l0,6)*pow(l1,6)*pow(r1,4) - 1785856*pow(k,4)*pow(l0,8)*pow(l1,6)*pow(r1,4) - 180224*pow(k,2)*pow(l0,10)*pow(l1,6)*pow(r1,4) + 4259840*pow(k,6)*pow(l0,4)*pow(l1,8)*pow(r1,4) + 917504*pow(k,4)*pow(l0,6)*pow(l1,8)*pow(r1,4) + 24576*pow(k,2)*pow(l0,8)*pow(l1,8)*pow(r1,4) - 2621440*pow(k,6)*pow(l0,2)*pow(l1,10)*pow(r1,4) - 786432*pow(k,4)*pow(l0,4)*pow(l1,10)*pow(r1,4) + 1048576*pow(k,6)*pow(l1,12)*pow(r1,4) - 122880*pow(k,6)*pow(l0,10)*pow(r1,6) - 61440*pow(k,4)*pow(l0,12)*pow(r1,6) - 7680*pow(k,2)*pow(l0,14)*pow(r1,6) + 57344*pow(k,6)*pow(l0,8)*pow(l1,2)*pow(r1,6) + 215040*pow(k,4)*pow(l0,10)*pow(l1,2)*pow(r1,6) + 50176*pow(k,2)*pow(l0,12)*pow(l1,2)*pow(r1,6) + 1572864*pow(k,6)*pow(l0,6)*pow(l1,4)*pow(r1,6) - 475136*pow(k,4)*pow(l0,8)*pow(l1,4)*pow(r1,6) - 270336*pow(k,2)*pow(l0,10)*pow(l1,4)*pow(r1,6) - 16384*pow(l0,12)*pow(l1,4)*pow(r1,6) - 7471104*pow(k,6)*pow(l0,4)*pow(l1,6)*pow(r1,6) - 2162688*pow(k,4)*pow(l0,6)*pow(l1,6)*pow(r1,6) - 229376*pow(k,2)*pow(l0,8)*pow(l1,6)*pow(r1,6) + 8388608*pow(k,6)*pow(l0,2)*pow(l1,8)*pow(r1,6) + 2490368*pow(k,4)*pow(l0,4)*pow(l1,8)*pow(r1,6) - 65536*pow(k,2)*pow(l0,6)*pow(l1,8)*pow(r1,6) - 6291456*pow(k,6)*pow(l1,10)*pow(r1,6) - 524288*pow(k,4)*pow(l0,2)*pow(l1,10)*pow(r1,6) - 303104*pow(k,6)*pow(l0,8)*pow(r1,8) - 151552*pow(k,4)*pow(l0,10)*pow(r1,8) - 18944*pow(k,2)*pow(l0,12)*pow(r1,8) + 1015808*pow(k,6)*pow(l0,6)*pow(l1,2)*pow(r1,8) + 516096*pow(k,4)*pow(l0,8)*pow(l1,2)*pow(r1,8) + 65536*pow(k,2)*pow(l0,10)*pow(l1,2)*pow(r1,8) + 3473408*pow(k,6)*pow(l0,4)*pow(l1,4)*pow(r1,8) + 2097152*pow(k,4)*pow(l0,6)*pow(l1,4)*pow(r1,8) + 401408*pow(k,2)*pow(l0,8)*pow(l1,4)*pow(r1,8) + 16384*pow(l0,10)*pow(l1,4)*pow(r1,8) - 7340032*pow(k,6)*pow(l0,2)*pow(l1,6)*pow(r1,8) - 1703936*pow(k,4)*pow(l0,4)*pow(l1,6)*pow(r1,8) + 262144*pow(k,2)*pow(l0,6)*pow(l1,6)*pow(r1,8) + 15728640*pow(k,6)*pow(l1,8)*pow(r1,8) + 2621440*pow(k,4)*pow(l0,2)*pow(l1,8)*pow(r1,8) + 65536*pow(k,2)*pow(l0,4)*pow(l1,8)*pow(r1,8) + 655360*pow(k,6)*pow(l0,6)*pow(r1,10) + 327680*pow(k,4)*pow(l0,8)*pow(r1,10) + 40960*pow(k,2)*pow(l0,10)*pow(r1,10) - 917504*pow(k,6)*pow(l0,4)*pow(l1,2)*pow(r1,10) - 1212416*pow(k,4)*pow(l0,6)*pow(l1,2)*pow(r1,10) - 245760*pow(k,2)*pow(l0,8)*pow(l1,2)*pow(r1,10) - 2097152*pow(k,6)*pow(l0,2)*pow(l1,4)*pow(r1,10) - 1966080*pow(k,4)*pow(l0,4)*pow(l1,4)*pow(r1,10) - 458752*pow(k,2)*pow(l0,6)*pow(l1,4)*pow(r1,10) - 20971520*pow(k,6)*pow(l1,6)*pow(r1,10) - 5242880*pow(k,4)*pow(l0,2)*pow(l1,6)*pow(r1,10) - 262144*pow(k,2)*pow(l0,4)*pow(l1,6)*pow(r1,10) + 786432*pow(k,6)*pow(l0,4)*pow(r1,12) + 393216*pow(k,4)*pow(l0,6)*pow(r1,12) + 49152*pow(k,2)*pow(l0,8)*pow(r1,12) + 5767168*pow(k,6)*pow(l0,2)*pow(l1,2)*pow(r1,12) + 3014656*pow(k,4)*pow(l0,4)*pow(l1,2)*pow(r1,12) + 393216*pow(k,2)*pow(l0,6)*pow(l1,2)*pow(r1,12) + 15728640*pow(k,6)*pow(l1,4)*pow(r1,12) + 5242880*pow(k,4)*pow(l0,2)*pow(l1,4)*pow(r1,12) + 393216*pow(k,2)*pow(l0,4)*pow(l1,4)*pow(r1,12) - 2097152*pow(k,6)*pow(l0,2)*pow(r1,14) - 1048576*pow(k,4)*pow(l0,4)*pow(r1,14) - 131072*pow(k,2)*pow(l0,6)*pow(r1,14) - 6291456*pow(k,6)*pow(l1,2)*pow(r1,14) - 2621440*pow(k,4)*pow(l0,2)*pow(l1,2)*pow(r1,14) - 262144*pow(k,2)*pow(l0,4)*pow(l1,2)*pow(r1,14) + 1048576*pow(k,6)*pow(r1,16) + 524288*pow(k,4)*pow(l0,2)*pow(r1,16) + 65536*pow(k,2)*pow(l0,4)*pow(r1,16),-27648*pow(k,6)*pow(l0,15) - 13824*pow(k,4)*pow(l0,17) - 1728*pow(k,2)*pow(l0,19) + 55296*pow(k,6)*pow(l0,13)*pow(l1,2) + 47616*pow(k,4)*pow(l0,15)*pow(l1,2) + 8448*pow(k,2)*pow(l0,17)*pow(l1,2) - 92160*pow(k,6)*pow(l0,11)*pow(l1,4) - 33280*pow(k,4)*pow(l0,13)*pow(l1,4) - 5120*pow(k,2)*pow(l0,15)*pow(l1,4) - 256*pow(l0,17)*pow(l1,4) + 57344*pow(k,6)*pow(l0,9)*pow(l1,6) + 24576*pow(k,4)*pow(l0,11)*pow(l1,6) - 32768*pow(k,6)*pow(l0,7)*pow(l1,8) - 8192*pow(k,4)*pow(l0,9)*pow(l1,8) - 221184*pow(k,6)*pow(l0,13)*pow(r1,2) - 110592*pow(k,4)*pow(l0,15)*pow(r1,2) - 13824*pow(k,2)*pow(l0,17)*pow(r1,2) + 1400832*pow(k,6)*pow(l0,11)*pow(l1,2)*pow(r1,2) + 759808*pow(k,4)*pow(l0,13)*pow(l1,2)*pow(r1,2) + 102400*pow(k,2)*pow(l0,15)*pow(l1,2)*pow(r1,2) - 2260992*pow(k,6)*pow(l0,9)*pow(l1,4)*pow(r1,2) - 1417216*pow(k,4)*pow(l0,11)*pow(l1,4)*pow(r1,2) - 188416*pow(k,2)*pow(l0,13)*pow(l1,4)*pow(r1,2) + 4096*pow(l0,15)*pow(l1,4)*pow(r1,2) + 2654208*pow(k,6)*pow(l0,7)*pow(l1,6)*pow(r1,2) + 704512*pow(k,4)*pow(l0,9)*pow(l1,6)*pow(r1,2) + 32768*pow(k,2)*pow(l0,11)*pow(l1,6)*pow(r1,2) - 1179648*pow(k,6)*pow(l0,5)*pow(l1,8)*pow(r1,2) - 458752*pow(k,4)*pow(l0,7)*pow(l1,8)*pow(r1,2) + 524288*pow(k,6)*pow(l0,3)*pow(l1,10)*pow(r1,2) - 147456*pow(k,6)*pow(l0,11)*pow(r1,4) - 73728*pow(k,4)*pow(l0,13)*pow(r1,4) - 9216*pow(k,2)*pow(l0,15)*pow(r1,4) + 4128768*pow(k,6)*pow(l0,9)*pow(l1,2)*pow(r1,4) + 1654784*pow(k,4)*pow(l0,11)*pow(l1,2)*pow(r1,4) + 155648*pow(k,2)*pow(l0,13)*pow(l1,2)*pow(r1,4) - 15925248*pow(k,6)*pow(l0,7)*pow(l1,4)*pow(r1,4) - 7323648*pow(k,4)*pow(l0,9)*pow(l1,4)*pow(r1,4) - 917504*pow(k,2)*pow(l0,11)*pow(l1,4)*pow(r1,4) - 24576*pow(l0,13)*pow(l1,4)*pow(r1,4) + 16908288*pow(k,6)*pow(l0,5)*pow(l1,6)*pow(r1,4) + 7864320*pow(k,4)*pow(l0,7)*pow(l1,6)*pow(r1,4) + 786432*pow(k,2)*pow(l0,9)*pow(l1,6)*pow(r1,4) - 12058624*pow(k,6)*pow(l0,3)*pow(l1,8)*pow(r1,4) - 2228224*pow(k,4)*pow(l0,5)*pow(l1,8)*pow(r1,4) + 1835008*pow(k,6)*pow(l0,9)*pow(r1,6) + 917504*pow(k,4)*pow(l0,11)*pow(r1,6) + 114688*pow(k,2)*pow(l0,13)*pow(r1,6) - 3407872*pow(k,6)*pow(l0,7)*pow(l1,2)*pow(r1,6) - 2686976*pow(k,4)*pow(l0,9)*pow(l1,2)*pow(r1,6) - 458752*pow(k,2)*pow(l0,11)*pow(l1,2)*pow(r1,6) - 4718592*pow(k,6)*pow(l0,5)*pow(l1,4)*pow(r1,6) + 1966080*pow(k,4)*pow(l0,7)*pow(l1,4)*pow(r1,6) + 917504*pow(k,2)*pow(l0,9)*pow(l1,4)*pow(r1,6) + 65536*pow(l0,11)*pow(l1,4)*pow(r1,6) + 28835840*pow(k,6)*pow(l0,3)*pow(l1,6)*pow(r1,6) + 8126464*pow(k,4)*pow(l0,5)*pow(l1,6)*pow(r1,6) + 524288*pow(k,2)*pow(l0,7)*pow(l1,6)*pow(r1,6) + 786432*pow(k,6)*pow(l0,7)*pow(r1,8) + 393216*pow(k,4)*pow(l0,9)*pow(r1,8) + 49152*pow(k,2)*pow(l0,11)*pow(r1,8) - 4718592*pow(k,6)*pow(l0,5)*pow(l1,2)*pow(r1,8) - 917504*pow(k,4)*pow(l0,7)*pow(l1,2)*pow(r1,8) + 65536*pow(k,2)*pow(l0,9)*pow(l1,2)*pow(r1,8) - 19398656*pow(k,6)*pow(l0,3)*pow(l1,4)*pow(r1,8) - 7471104*pow(k,4)*pow(l0,5)*pow(l1,4)*pow(r1,8) - 786432*pow(k,2)*pow(l0,7)*pow(l1,4)*pow(r1,8) - 65536*pow(l0,9)*pow(l1,4)*pow(r1,8) - 6291456*pow(k,6)*pow(l0,5)*pow(r1,10) - 3145728*pow(k,4)*pow(l0,7)*pow(r1,10) - 393216*pow(k,2)*pow(l0,9)*pow(r1,10) - 2097152*pow(k,6)*pow(l0,3)*pow(l1,2)*pow(r1,10) - 524288*pow(k,4)*pow(l0,5)*pow(l1,2)*pow(r1,10) + 4194304*pow(k,6)*pow(l0,3)*pow(r1,12) + 2097152*pow(k,4)*pow(l0,5)*pow(r1,12) + 262144*pow(k,2)*pow(l0,7)*pow(r1,12),248832*pow(k,6)*pow(l0,14) + 124416*pow(k,4)*pow(l0,16) + 15552*pow(k,2)*pow(l0,18) - 350208*pow(k,6)*pow(l0,12)*pow(l1,2) - 285184*pow(k,4)*pow(l0,14)*pow(l1,2) - 49408*pow(k,2)*pow(l0,16)*pow(l1,2) + 485376*pow(k,6)*pow(l0,10)*pow(l1,4) + 147968*pow(k,4)*pow(l0,12)*pow(l1,4) + 13312*pow(k,2)*pow(l0,14)*pow(l1,4) + 256*pow(l0,16)*pow(l1,4) - 188416*pow(k,6)*pow(l0,8)*pow(l1,6) - 90112*pow(k,4)*pow(l0,10)*pow(l1,6) + 98304*pow(k,6)*pow(l0,6)*pow(l1,8) + 8192*pow(k,4)*pow(l0,8)*pow(l1,8) + 1400832*pow(k,6)*pow(l0,12)*pow(r1,2) + 700416*pow(k,4)*pow(l0,14)*pow(r1,2) + 87552*pow(k,2)*pow(l0,16)*pow(r1,2) - 8085504*pow(k,6)*pow(l0,10)*pow(l1,2)*pow(r1,2) - 4069376*pow(k,4)*pow(l0,12)*pow(l1,2)*pow(r1,2) - 512000*pow(k,2)*pow(l0,14)*pow(l1,2)*pow(r1,2) + 8028160*pow(k,6)*pow(l0,8)*pow(l1,4)*pow(r1,2) + 4825088*pow(k,4)*pow(l0,10)*pow(l1,4)*pow(r1,2) + 647168*pow(k,2)*pow(l0,12)*pow(l1,4)*pow(r1,2) - 4096*pow(l0,14)*pow(l1,4)*pow(r1,2) - 7372800*pow(k,6)*pow(l0,6)*pow(l1,6)*pow(r1,2) - 1753088*pow(k,4)*pow(l0,8)*pow(l1,6)*pow(r1,2) - 32768*pow(k,2)*pow(l0,10)*pow(l1,6)*pow(r1,2) + 1179648*pow(k,6)*pow(l0,4)*pow(l1,8)*pow(r1,2) + 458752*pow(k,4)*pow(l0,6)*pow(l1,8)*pow(r1,2) - 524288*pow(k,6)*pow(l0,2)*pow(l1,10)*pow(r1,2) - 638976*pow(k,6)*pow(l0,10)*pow(r1,4) - 319488*pow(k,4)*pow(l0,12)*pow(r1,4) - 39936*pow(k,2)*pow(l0,14)*pow(r1,4) - 14090240*pow(k,6)*pow(l0,8)*pow(l1,2)*pow(r1,4) - 5193728*pow(k,4)*pow(l0,10)*pow(l1,2)*pow(r1,4) - 417792*pow(k,2)*pow(l0,12)*pow(l1,2)*pow(r1,4) + 47382528*pow(k,6)*pow(l0,6)*pow(l1,4)*pow(r1,4) + 19644416*pow(k,4)*pow(l0,8)*pow(l1,4)*pow(r1,4) + 2097152*pow(k,2)*pow(l0,10)*pow(l1,4)*pow(r1,4) + 24576*pow(l0,12)*pow(l1,4)*pow(r1,4) - 16908288*pow(k,6)*pow(l0,4)*pow(l1,6)*pow(r1,4) - 7864320*pow(k,4)*pow(l0,6)*pow(l1,6)*pow(r1,4) - 786432*pow(k,2)*pow(l0,8)*pow(l1,6)*pow(r1,4) + 12058624*pow(k,6)*pow(l0,2)*pow(l1,8)*pow(r1,4) + 2228224*pow(k,4)*pow(l0,4)*pow(l1,8)*pow(r1,4) - 8126464*pow(k,6)*pow(l0,8)*pow(r1,6) - 4063232*pow(k,4)*pow(l0,10)*pow(r1,6) - 507904*pow(k,2)*pow(l0,12)*pow(r1,6) + 18087936*pow(k,6)*pow(l0,6)*pow(l1,2)*pow(r1,6) + 9502720*pow(k,4)*pow(l0,8)*pow(l1,2)*pow(r1,6) + 1245184*pow(k,2)*pow(l0,10)*pow(l1,2)*pow(r1,6) + 4718592*pow(k,6)*pow(l0,4)*pow(l1,4)*pow(r1,6) - 1966080*pow(k,4)*pow(l0,6)*pow(l1,4)*pow(r1,6) - 917504*pow(k,2)*pow(l0,8)*pow(l1,4)*pow(r1,6) - 65536*pow(l0,10)*pow(l1,4)*pow(r1,6) - 28835840*pow(k,6)*pow(l0,2)*pow(l1,6)*pow(r1,6) - 8126464*pow(k,4)*pow(l0,4)*pow(l1,6)*pow(r1,6) - 524288*pow(k,2)*pow(l0,6)*pow(l1,6)*pow(r1,6) + 5505024*pow(k,6)*pow(l0,6)*pow(r1,8) + 2752512*pow(k,4)*pow(l0,8)*pow(r1,8) + 344064*pow(k,2)*pow(l0,10)*pow(r1,8) + 4718592*pow(k,6)*pow(l0,4)*pow(l1,2)*pow(r1,8) + 917504*pow(k,4)*pow(l0,6)*pow(l1,2)*pow(r1,8) - 65536*pow(k,2)*pow(l0,8)*pow(l1,2)*pow(r1,8) + 19398656*pow(k,6)*pow(l0,2)*pow(l1,4)*pow(r1,8) + 7471104*pow(k,4)*pow(l0,4)*pow(l1,4)*pow(r1,8) + 786432*pow(k,2)*pow(l0,6)*pow(l1,4)*pow(r1,8) + 65536*pow(l0,8)*pow(l1,4)*pow(r1,8) + 6291456*pow(k,6)*pow(l0,4)*pow(r1,10) + 3145728*pow(k,4)*pow(l0,6)*pow(r1,10) + 393216*pow(k,2)*pow(l0,8)*pow(r1,10) + 2097152*pow(k,6)*pow(l0,2)*pow(l1,2)*pow(r1,10) + 524288*pow(k,4)*pow(l0,4)*pow(l1,2)*pow(r1,10) - 4194304*pow(k,6)*pow(l0,2)*pow(r1,12) - 2097152*pow(k,4)*pow(l0,4)*pow(r1,12) - 262144*pow(k,2)*pow(l0,6)*pow(r1,12),-1228800*pow(k,6)*pow(l0,13) - 614400*pow(k,4)*pow(l0,15) - 76800*pow(k,2)*pow(l0,17) + 1114112*pow(k,6)*pow(l0,11)*pow(l1,2) + 868352*pow(k,4)*pow(l0,13)*pow(l1,2) + 147456*pow(k,2)*pow(l0,15)*pow(l1,2) - 1310720*pow(k,6)*pow(l0,9)*pow(l1,4) - 360448*pow(k,4)*pow(l0,11)*pow(l1,4) - 16384*pow(k,2)*pow(l0,13)*pow(l1,4) + 262144*pow(k,6)*pow(l0,7)*pow(l1,6) + 131072*pow(k,4)*pow(l0,9)*pow(l1,6) - 131072*pow(k,6)*pow(l0,5)*pow(l1,8) - 4456448*pow(k,6)*pow(l0,11)*pow(r1,2) - 2228224*pow(k,4)*pow(l0,13)*pow(r1,2) - 278528*pow(k,2)*pow(l0,15)*pow(r1,2) + 23855104*pow(k,6)*pow(l0,9)*pow(l1,2)*pow(r1,2) + 11337728*pow(k,4)*pow(l0,11)*pow(l1,2)*pow(r1,2) + 1343488*pow(k,2)*pow(l0,13)*pow(l1,2)*pow(r1,2) - 11534336*pow(k,6)*pow(l0,7)*pow(l1,4)*pow(r1,2) - 6815744*pow(k,4)*pow(l0,9)*pow(l1,4)*pow(r1,2) - 917504*pow(k,2)*pow(l0,11)*pow(l1,4)*pow(r1,2) + 9437184*pow(k,6)*pow(l0,5)*pow(l1,6)*pow(r1,2) + 2097152*pow(k,4)*pow(l0,7)*pow(l1,6)*pow(r1,2) + 5767168*pow(k,6)*pow(l0,9)*pow(r1,4) + 2883584*pow(k,4)*pow(l0,11)*pow(r1,4) + 360448*pow(k,2)*pow(l0,13)*pow(r1,4) + 19922944*pow(k,6)*pow(l0,7)*pow(l1,2)*pow(r1,4) + 7077888*pow(k,4)*pow(l0,9)*pow(l1,2)*pow(r1,4) + 524288*pow(k,2)*pow(l0,11)*pow(l1,2)*pow(r1,4) - 62914560*pow(k,6)*pow(l0,5)*pow(l1,4)*pow(r1,4) - 24641536*pow(k,4)*pow(l0,7)*pow(l1,4)*pow(r1,4) - 2359296*pow(k,2)*pow(l0,9)*pow(l1,4)*pow(r1,4) + 12582912*pow(k,6)*pow(l0,7)*pow(r1,6) + 6291456*pow(k,4)*pow(l0,9)*pow(r1,6) + 786432*pow(k,2)*pow(l0,11)*pow(r1,6) - 29360128*pow(k,6)*pow(l0,5)*pow(l1,2)*pow(r1,6) - 13631488*pow(k,4)*pow(l0,7)*pow(l1,2)*pow(r1,6) - 1572864*pow(k,2)*pow(l0,9)*pow(l1,2)*pow(r1,6) - 12582912*pow(k,6)*pow(l0,5)*pow(r1,8) - 6291456*pow(k,4)*pow(l0,7)*pow(r1,8) - 786432*pow(k,2)*pow(l0,9)*pow(r1,8),3629056*pow(k,6)*pow(l0,12) + 1814528*pow(k,4)*pow(l0,14) + 226816*pow(k,2)*pow(l0,16) - 1867776*pow(k,6)*pow(l0,10)*pow(l1,2) - 1417216*pow(k,4)*pow(l0,12)*pow(l1,2) - 237568*pow(k,2)*pow(l0,14)*pow(l1,2) + 1966080*pow(k,6)*pow(l0,8)*pow(l1,4) + 507904*pow(k,4)*pow(l0,10)*pow(l1,4) + 8192*pow(k,2)*pow(l0,12)*pow(l1,4) - 131072*pow(k,6)*pow(l0,6)*pow(l1,6) - 65536*pow(k,4)*pow(l0,8)*pow(l1,6) + 65536*pow(k,6)*pow(l0,4)*pow(l1,8) + 7471104*pow(k,6)*pow(l0,10)*pow(r1,2) + 3735552*pow(k,4)*pow(l0,12)*pow(r1,2) + 466944*pow(k,2)*pow(l0,14)*pow(r1,2) - 38141952*pow(k,6)*pow(l0,8)*pow(l1,2)*pow(r1,2) - 17465344*pow(k,4)*pow(l0,10)*pow(l1,2)*pow(r1,2) - 1982464*pow(k,2)*pow(l0,12)*pow(l1,2)*pow(r1,2) + 5767168*pow(k,6)*pow(l0,6)*pow(l1,4)*pow(r1,2) + 3407872*pow(k,4)*pow(l0,8)*pow(l1,4)*pow(r1,2) + 458752*pow(k,2)*pow(l0,10)*pow(l1,4)*pow(r1,2) - 4718592*pow(k,6)*pow(l0,4)*pow(l1,6)*pow(r1,2) - 1048576*pow(k,4)*pow(l0,6)*pow(l1,6)*pow(r1,2) - 13369344*pow(k,6)*pow(l0,8)*pow(r1,4) - 6684672*pow(k,4)*pow(l0,10)*pow(r1,4) - 835584*pow(k,2)*pow(l0,12)*pow(r1,4) - 9961472*pow(k,6)*pow(l0,6)*pow(l1,2)*pow(r1,4) - 3538944*pow(k,4)*pow(l0,8)*pow(l1,2)*pow(r1,4) - 262144*pow(k,2)*pow(l0,10)*pow(l1,2)*pow(r1,4) + 31457280*pow(k,6)*pow(l0,4)*pow(l1,4)*pow(r1,4) + 12320768*pow(k,4)*pow(l0,6)*pow(l1,4)*pow(r1,4) + 1179648*pow(k,2)*pow(l0,8)*pow(l1,4)*pow(r1,4) - 6291456*pow(k,6)*pow(l0,6)*pow(r1,6) - 3145728*pow(k,4)*pow(l0,8)*pow(r1,6) - 393216*pow(k,2)*pow(l0,10)*pow(r1,6) + 14680064*pow(k,6)*pow(l0,4)*pow(l1,2)*pow(r1,6) + 6815744*pow(k,4)*pow(l0,6)*pow(l1,2)*pow(r1,6) + 786432*pow(k,2)*pow(l0,8)*pow(l1,2)*pow(r1,6) + 6291456*pow(k,6)*pow(l0,4)*pow(r1,8) + 3145728*pow(k,4)*pow(l0,6)*pow(r1,8) + 393216*pow(k,2)*pow(l0,8)*pow(r1,8),-6553600*pow(k,6)*pow(l0,11) - 3276800*pow(k,4)*pow(l0,13) - 409600*pow(k,2)*pow(l0,15) + 1572864*pow(k,6)*pow(l0,9)*pow(l1,2) + 1179648*pow(k,4)*pow(l0,11)*pow(l1,2) + 196608*pow(k,2)*pow(l0,13)*pow(l1,2) - 1572864*pow(k,6)*pow(l0,7)*pow(l1,4) - 393216*pow(k,4)*pow(l0,9)*pow(l1,4) - 6291456*pow(k,6)*pow(l0,9)*pow(r1,2) - 3145728*pow(k,4)*pow(l0,11)*pow(r1,2) - 393216*pow(k,2)*pow(l0,13)*pow(r1,2) + 31457280*pow(k,6)*pow(l0,7)*pow(l1,2)*pow(r1,2) + 14155776*pow(k,4)*pow(l0,9)*pow(l1,2)*pow(r1,2) + 1572864*pow(k,2)*pow(l0,11)*pow(l1,2)*pow(r1,2) + 12582912*pow(k,6)*pow(l0,7)*pow(r1,4) + 6291456*pow(k,4)*pow(l0,9)*pow(r1,4) + 786432*pow(k,2)*pow(l0,11)*pow(r1,4),7077888*pow(k,6)*pow(l0,10) + 3538944*pow(k,4)*pow(l0,12) + 442368*pow(k,2)*pow(l0,14) - 524288*pow(k,6)*pow(l0,8)*pow(l1,2) - 393216*pow(k,4)*pow(l0,10)*pow(l1,2) - 65536*pow(k,2)*pow(l0,12)*pow(l1,2) + 524288*pow(k,6)*pow(l0,6)*pow(l1,4) + 131072*pow(k,4)*pow(l0,8)*pow(l1,4) + 2097152*pow(k,6)*pow(l0,8)*pow(r1,2) + 1048576*pow(k,4)*pow(l0,10)*pow(r1,2) + 131072*pow(k,2)*pow(l0,12)*pow(r1,2) - 10485760*pow(k,6)*pow(l0,6)*pow(l1,2)*pow(r1,2) - 4718592*pow(k,4)*pow(l0,8)*pow(l1,2)*pow(r1,2) - 524288*pow(k,2)*pow(l0,10)*pow(l1,2)*pow(r1,2) - 4194304*pow(k,6)*pow(l0,6)*pow(r1,4) - 2097152*pow(k,4)*pow(l0,8)*pow(r1,4) - 262144*pow(k,2)*pow(l0,10)*pow(r1,4),-4194304*pow(k,6)*pow(l0,9) - 2097152*pow(k,4)*pow(l0,11) - 262144*pow(k,2)*pow(l0,13),1048576*pow(k,6)*pow(l0,8) + 524288*pow(k,4)*pow(l0,10) + 65536*pow(k,2)*pow(l0,12)};

  gsl_poly_complex_workspace * v
      = gsl_poly_complex_workspace_alloc (9);

  gsl_poly_complex_solve (a, 9, v, u);

  gsl_poly_complex_workspace_free (v);

return 0;
}

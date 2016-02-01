

#include <iostream>
#include "finite_difference.h"

using namespace lars;
using namespace std;


complex<double> F(double x,double y,double z){
  return (-(complex<double>(0.0,3040.638455120569)*((pow((((6000<z)|((0.0<2)&((5.*abs(x))<2)&(0<=x)))?(1):(true)?(complex<double>(0.999993994,-6.32e-07)):0),2.))+(-1.))));
}

complex<double> u_boundary(double x,double y,double z){
  return (exp(-(complex<double>(0.0,3040.638455120569)*z*((pow(((((5*abs(x))<2)&((5.*abs(y))<2)&(0<=x))?(1):(true)?(complex<double>(0.999993994,-6.32e-07)):0.),2.))+(-1.)))));
}


int main(){
  
  finite_difference_1D solver;
  
  solver.set_field(finite_difference_1D::vector::Ones(600));
  solver.F = [](double x,double z){ return F(x,0,z); };
  solver.u_boundary = [](double x,double z){ return u_boundary(x,0,z); };
  solver.A = finite_difference_1D::complex(0,-8.221957450383127e-05);
  solver.xmin = -1;
  solver.xmax = 1;
  solver.dz = 40./3;
  solver.z = 0;
  
  
  for(int i=0;i<1000;++i){
    std::cout << abs(solver.get_field().mean()) << "\t" << abs(solver.get_field()(0))  << "\t" << abs(solver.get_field()(1)) << std::endl;

    solver.step();
  }
  
}


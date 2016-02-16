

#include <iostream>
#include "crank_nicolson.h"
#include "finite_difference.h"

using namespace lars;
using namespace std;

auto A = complex<double>(0,-0.00008221957450383126835);
auto F = complex<double>(-0.00060812829915180474809, - 0.0060812799204732078091);
auto dz = 16.032064128256513026;
auto dx = 0.02020202020202020202;
auto dy = 0.0067114093959731543624;

void init(crank_nicolson_2D &C){
  
  auto dt = dz/2;
  
  std::cout << "ra: " << A * dt/(dy*dy) << std::endl;
  std::cout << "rb: " << A * dt/(dx*dx) << std::endl;
  

  C.ra.fill(A * dt/(dy*dy));
  C.rc.fill(A * dt/(dx*dx));
  C.rd.fill(0);
  C.re.fill(0);
  C.rf.fill(F * dt/2.);
  C.u.fill(1);
  
}

std::complex<double> u_boundary(double x, double y, double z){
  return exp(-z*std::complex<double>(0.0019357315193375983642, - 0.00019357324968815317454j)* std::complex<double>(0.0 , 3.1415926535897932385j));
}

int main(){
  
  crank_nicolson_2D C;
  finite_difference_2D FD;
  
  C.resize(100, 150);
  
  init(C);
  C.update();
  init(C);

  FD.dz = dz;
  FD.A = A;
  FD.F = [=](double x,double y,double z){ return F; };
  FD.u_boundary = u_boundary;
  FD.set_field(C.u.transpose());
  
  FD.xmin = -1;
  FD.xmax = 1;
  FD.ymin = -0.5;
  FD.ymax = 0.5;
  
  for(int i = 0;i<100;++i){
    FD.step();
    
    C.step_1();
    C.update();
    C.step_2();
    C.update();
    
    std::cout << pow(abs(FD.get_field()),2).mean() - pow(abs(FD.get_field()).mean(),2) << std::endl;
    //std::cout << abs(FD.get_field().transpose() - C.u).mean() << std::endl;
  }
  
}


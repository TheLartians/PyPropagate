

#include <iostream>
#include "crank_nicolson.h"
#include "finite_difference.h"

using namespace lars;
using namespace std;


void init(crank_nicolson_2D &C){
  C.ra.fill(crank_nicolson_2D::scalar(0,1));
  C.rb.fill(crank_nicolson_2D::scalar(0,1));
  C.rc.fill(0);
  C.rd.fill(0);
  C.re.fill(0);
  C.rf.fill(0);
  C.u.fill(1);
}

int main(){
  
  crank_nicolson_2D C;
  finite_difference_2D F;
  
  C.resize(5, 10);
  
  init(C);
  C.update();
  init(C);
  

  F.dz = 2;
  F.A = finite_difference_2D::complex(0,1);
  F.F = [](double x,double y,double z){ return 0; };
  F.u_boundary = [](double x,double y,double z){ return 1; };
  F.set_field(C.u.transpose());
  F.xmin = 0;
  F.ymin = 0;
  F.ymax = C.u.rows()+1;
  F.xmax = C.u.cols()+1;
  

  for(int i = 0;i<100;++i){
    F.step();
    
    C.step_1();
    C.update();
    C.step_2();
    C.update();

    std::cout << F.get_field()(5,2) << "\t" << C.u(2,5) << std::endl;
  }
  
}


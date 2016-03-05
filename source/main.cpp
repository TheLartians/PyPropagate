

#include <chrono>
#include <iostream>
#include "finite_difference.h"

using namespace lars;

int main(){
  
  finite_difference_aF propagator;
  
  auto boundary = 0;
  
  propagator.resize(10);
  propagator.u.fill(1);
  propagator.u[0] = boundary;
  propagator.u[9] = boundary;
  propagator.rf.fill(finite_differences::complex(0,0));
  
  propagator.update();
  propagator.u.fill(1);
  propagator.u[0] = boundary;
  propagator.u[9] = boundary;
  propagator.rf.fill(finite_differences::complex(0,0));
  
  propagator.ra = finite_differences::complex(1,0);
  
  for(auto i UNUSED:range(10)){
    propagator.update();
    propagator.step();
    ndarray<double, 1> intensity(propagator.u.size());
    intensity.element_wise([&](dynamic_index_tuple<1> idx){ return abs(propagator.u(idx))*abs(propagator.u(idx)) ; });
    //std::cout << propagator.u << std::endl;
    std::cout << intensity << std::endl;
  }
  
  return 0;
}
  

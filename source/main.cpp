

#include <chrono>
#include <iostream>
#include "finite_difference.h"

using namespace lars;

void test_1D(){
  
  finite_difference_aF propagator;
  
  propagator.resize(10);
  propagator.u.fill(1);
  propagator.rf.fill(finite_differences::complex(0,0));
  
  propagator.update();
  propagator.u.fill(1);
  propagator.rf.fill(finite_differences::complex(0,0));
  
  propagator.ra = finite_differences::complex(0,2);
  
  for(auto i UNUSED:range(10)){
    propagator.update();
    propagator.step();
    ndarray<double, 1> intensity(propagator.u.size());
    intensity.element_wise([&](dynamic_index_tuple<1> idx){ return abs(propagator.u(idx))*abs(propagator.u(idx)) ; });
    std::cout << intensity << std::endl;
  }

}

void test_2D(){
  

  finite_difference_acF propagator;
  
  propagator.resize(50,10);
  propagator.u.fill(1);
  propagator.rf.fill(finite_differences::complex(0,0));
  
  propagator.update();
  propagator.u.fill(1);
  propagator.rf.fill(finite_differences::complex(0,0));
  
  propagator.ra = finite_differences::complex(0,1);
  propagator.rc = finite_differences::complex(0,1);
  
  for(auto i UNUSED:range(10)){
    propagator.update();
    propagator.step_1();
    propagator.update();
    propagator.step_2();
    
    using real_array = ndarray<double, propagator.u.ndim()>;
    real_array intensity(propagator.u.shape());
    intensity.element_wise([&](real_array::index_type idx){ return abs(propagator.u(idx))*abs(propagator.u(idx)) ; });
    std::cout << intensity << std::endl << std::endl;
  }
  
}

int main(){
  
  test_1D();
  std::cout << "\nTest 2D:" << std::endl;
  test_2D();
  
  return 0;
}
  

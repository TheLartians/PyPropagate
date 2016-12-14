

#include <chrono>
#include <iostream>

#include <lars/timeit.h>
#include "finite_difference.h"

using namespace lars;

int main(){
  

  finite_difference_ABC propagator;  
  propagator.resize(10);

  auto dz = 0.01;

  propagator.u.fill(1);
  propagator.ra.fill(1);
  propagator.rb.fill(0);
  propagator.rc.fill(0);
  propagator.rz.fill(1);

  propagator.update();

  propagator.u.fill(1);
  propagator.ra.fill(1);
  propagator.rb.fill(0);
  propagator.rc.fill(0);
  propagator.rz.fill(1);

  std::cout << "initial: " << propagator.u.transpose() << std::endl;

  for(auto i:range(5)){
    propagator.update();
    propagator.step();
    std::cout << "step " << i << ": " << propagator.u.transpose() << std::endl;
  }

  return 0;
  
  timeit<100>("AF", [](){
    finite_difference_AF propagator;
    propagator.resize(1000);
  
    for(auto i UNUSED:range(1000)){
      propagator.update();
      propagator.step();
    }
    
  });
  
  timeit<2>("A0F", [](){
    finite_difference_A0F propagator;
    propagator.resize(1000, 100);
    for(auto i UNUSED:range(1000)){
      propagator.update();
      propagator.step();
    }
  });

  timeit<2>("ACF", [](){
    finite_difference_ACF propagator;
    propagator.resize(1000, 100);
    for(auto i UNUSED:range(1000)){
      propagator.update();
      propagator.step_1();
      propagator.update();
      propagator.step_2();
    }
  });
  
  return 0;
}
  

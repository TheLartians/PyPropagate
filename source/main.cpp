

#include <chrono>
#include <iostream>

#include <lars/timeit.h>
#include "finite_difference.h"

using namespace lars;

int main(){
  
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
  

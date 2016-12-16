

#include <chrono>
#include <iostream>

#include <lars/timeit.h>
#include "finite_difference.h"

using namespace lars;

int main(){
  
  timeit<2>("AF", [](){
    finite_difference_AF propagator;
    propagator.resize(10);
  
    for(auto i UNUSED:range(10)){
      propagator.update();
      propagator.step();
    }
    
  });
  
  timeit<2>("ABC", [](){
    finite_difference_ABC propagator;
    propagator.resize(10);
  
    for(auto i UNUSED:range(10)){
      propagator.update();
      propagator.step();
    }
    
  });
  
  timeit<2>("A0F", [](){
    finite_difference_A0F propagator;
    propagator.resize(10, 5);
    for(auto i UNUSED:range(10)){
      propagator.update();
      propagator.step();
    }
  });

  timeit<2>("ACF", [](){
    finite_difference_ACF propagator;
    propagator.resize(10, 5);
    for(auto i UNUSED:range(10)){
      propagator.update();
      propagator.step_1();
      propagator.update();
      propagator.step_2();
    }
  });
  
  return 0;
}
  

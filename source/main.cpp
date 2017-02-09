
#include <chrono>
#include <iostream>

#include <lars/timeit.h>
#include "finite_difference.h"
#include "ring_derivative.h"

using namespace lars;

int main(){

  HeapNDArray<double,DynamicIndexTuple<2>> test;
  test.resize(5,10);
  
  test.for_all_indices([&](DynamicIndexTuple<2> idx){ test(idx) = fmod(1*idx.get<0>() - 0.5 * idx.get<1>(),M_PI) - 1; });

  std::cout << test << std::endl << std::endl;

  std::cout << ring_derivative<0>(test,-1,M_PI-1) << std::endl;
  std::cout << ring_derivative<1>(test,-1,M_PI-1) << std::endl;

  return 0;
  
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
  

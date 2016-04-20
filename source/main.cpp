

#include <chrono>
#include <iostream>
#include "finite_difference.h"

using namespace lars;

int main(){
  
  finite_difference_a0F propagator;
  propagator.resize(1000, 100);
  propagator.update();
  propagator.step();
  propagator.update();
  propagator.step();

  return 0;
}
  

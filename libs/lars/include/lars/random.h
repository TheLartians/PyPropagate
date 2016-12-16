#pragma once

#include <random>
#include <initializer_list>

namespace lars{

  inline std::default_random_engine & random_number_engine(){
    static std::default_random_engine engine;
    static bool initialized = false;
    if(initialized) return engine;
    std::random_device device;
    auto seed = device();
    engine.seed(seed);
    initialized = true;
    return engine;
  }
  
  template <class I = int> I uniform_random_int(I min, I max){
    std::uniform_int_distribution<I> dis(min,max);
    return dis(random_number_engine());
  }
  
  inline bool random_bool(float probability = 0.5){
    std::bernoulli_distribution dis(probability);
    return dis(random_number_engine());
  }
  
  template <class C> size_t random_index(const C &c){
    return uniform_random_int<size_t>(0,c.size()-1);
  }
  
  template <class C> typename C::value_type random_element(const C &c){
    return c[random_index(c)];
  }

  template <class T> T random_element(const std::initializer_list<T> &i){
    return random_element(std::vector<T>(i));
  }
  
}

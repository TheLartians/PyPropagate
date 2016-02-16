#pragma once

#include "algebra.h"
#include <memory>

namespace lars {
  
  namespace crank_nicolson{
    using real = double;
    using scalar = algebra::complex<real>;
    using vector = algebra::vector<scalar>;
    using array_1D = algebra::array<scalar,algebra::dynamic_size,1>;
    using array_2D = algebra::array<scalar,algebra::dynamic_size,algebra::dynamic_size>;
  }
  
  class crank_nicolson_2D{
    
  public:
    using real = crank_nicolson::real;
    using scalar = crank_nicolson::scalar;
    using field2D = crank_nicolson::array_2D;
    using field1D = crank_nicolson::array_1D;
    
  private:
    
    struct parallel_data{
      field1D A,B,C,R,U,tmp;
      parallel_data(unsigned s):A(s),B(s),C(s),R(s),U(s),tmp(s){}
    };
    
    field2D rap,rcp,rdp,rep,rfp,up;
    scalar v1(unsigned i,unsigned j);
    scalar v2(unsigned i,unsigned j);
    
  public:
    
    field2D ra,rc,rd,re,rf,u;
    
    void resize(unsigned Nx,unsigned Ny);
    
    void step_1();
    void update();
    void step_2();
    
  };
  
  
}

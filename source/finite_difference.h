/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//    Finite differences c++ solver for parabolic differential equations
//
//
//    File created in 2012 by Lars Melchior for the Institute for X-Ray Physics, Göttingen
//
//

#pragma once

#include "algebra.h"
#include <memory>

namespace lars {
  
  namespace finite_differences{
    using real = double;
    using complex = algebra::complex<real>;

    using array_1D   = ndarray<complex, 1>;
    using array_2D   = ndarray<complex, 2>;
    
    class pseudo_field{
      private:
      complex _value;
      public:
      pseudo_field(const complex & value):_value(value){}
      const complex &operator[](size_t i)const{ return _value; }
    };
    
    struct tridiagonal_data{
      array_1D A,B,C,R,U,tmp;
      tridiagonal_data(unsigned s):A(s),B(s),C(s),R(s),U(s),tmp(s){}
    };
    
  }
  
  class finite_difference_aF{
  public:
    
    using scalar = finite_differences::complex;
    using complex = finite_differences::complex;
    using field = finite_differences::array_1D;
    
  private:
    
    field up,rfp;
    field B,R,tmp;
    
  public:
    
    scalar ra;
    field u,rf;
    
    void step();
    void update();
    
    void resize(size_t N);
  };
  
  class finite_difference_acF{
  public:
    
    using scalar = finite_differences::complex;
    using complex = finite_differences::complex;
    using field = finite_differences::array_2D;
    
  private:
    
    field up,rfp;
    
  public:
    
    scalar ra,rc;
    field u,rf;
    
    void step_1();
    void step_2();
    void update();
    
    void resize(size_t nx,size_t ny);
  };
  
  
}








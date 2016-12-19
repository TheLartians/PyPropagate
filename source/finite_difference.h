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
#include <lars/parallel.h>

namespace lars {
  
  namespace finite_differences{
    
    using real = double;
    using complex = algebra::complex<real>;
    
    using scalar = complex;

    using array_1D   = HeapNDArray<scalar,DynamicIndexTuple<1>>;
    using array_2D   = HeapNDArray<scalar,DynamicIndexTuple<2>>;
    
    class pseudo_field{
      private: 
      scalar _value;
      public:
      pseudo_field(const scalar & value):_value(value){}
      const scalar &operator[](size_t i)const{ return _value; }
    };
    
    struct tridiagonal_data{
      array_1D A,B,C,R,U,tmp;
      tridiagonal_data(unsigned s):A(s),B(s),C(s),R(s),U(s),tmp(s){}
    };

  }
  
  /*
   
   ra = A dz/dx^2
   rc = C dz/dy^2
   rf = F dz/2
   
  */
  
  class finite_difference_AF{
  public:
    
    using scalar = finite_differences::complex;
    using field = finite_differences::array_1D;
    
  private:
    
    field up,rfp,rap;
    field A,B,R,tmp;
    
  public:
    
    field u,rf,ra;
    
    void step();
    void update();
    
    void resize(size_t N);
  };
  
  class finite_difference_ACF{
  public:
    
    using scalar = finite_differences::complex;
    using field = finite_differences::array_2D;
    
  private:
    
    field up,rfp,rap,rcp;
    
  public:
    
    field u,rf,ra,rc;
    unsigned thread_count = hardware_thread_count();
    
    void step_1();
    void step_2();
    void update();
    
    void resize(size_t nx,size_t ny);
  };
  
  
  class finite_difference_A0F{
  public:
    
    using scalar = finite_differences::complex;
    using field = finite_differences::array_2D;
    
  private:
    
    field up,rfp,rap;
    
  public:
    
    field u,rf,ra;
    unsigned thread_count = hardware_thread_count();
    
    void step();
    void update();
    
    void resize(size_t nx,size_t ny);
  };
  
  
  class finite_difference_ABC{
  public:
    
    using scalar = finite_differences::complex;
    using field = finite_differences::array_2D;
    
  private:
    
    field up,rap,rbp,rcp,rzp; 
    
  public:
    // alpha dz u = A dr^2 u + B d_r u + C u 
    // ra = A/2Dr^2, rb = B/4Dr, rc = C/2, rz = alpha/Dz
    field u,ra,rb,rc,rz; 
    unsigned thread_count = hardware_thread_count();
    
    void step();
    void update();
    
    void resize(size_t nx,size_t ny=1);
  };
  

}








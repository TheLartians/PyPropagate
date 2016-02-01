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
    using vector  = algebra::vector<complex>;
    using array_1D   = algebra::array<complex,algebra::dynamic_size,1>;
    using array_2D   = algebra::array<complex,algebra::dynamic_size,algebra::dynamic_size>;
  }
  
  class finite_difference_1D{
  public:
    using real = finite_differences::real;
    using complex = finite_differences::complex;
    using vector = finite_differences::vector;
    
  private:
    vector Ax,Bx,Dx,tmp;
    complex rx;
    bool ready;
    unsigned s;
    vector field;
    real dx;

  public:
    real xmin,xmax;
    real dz,z;
    
    complex A;
    
    std::function<complex(real,real)> F;
    std::function<complex(real,real)> u_boundary;
    
    const vector & get_field()const{ return field; }
    void set_field(const vector::Base &f){ ready = false; field = f; }
    void resize(unsigned nx){ field.resize(nx); }
    
    finite_difference_1D();
    void init();
    
    void step();
  };
  
  class finite_difference_2D{
  public:
    using real = finite_differences::real;
    using complex = finite_differences::complex;
    using vector = finite_differences::vector;
    using field = finite_differences::array_2D;
    
  private:
    field CField1,CField2;
    field field1,field2;
    vector Ax,Ay,&Cx,&Cy;
    complex rx,ry;
    bool ready;
    int sx,sy;
    real dx,dy;
    
    void update();
    
    complex C(int xi,int yi ,int zi);
    complex &u(int xi,int yi,int zi);
    
  public:
    real xmin,xmax,ymin,ymax,dz,z;
    complex A;
    
    std::function<complex(real,real,real)> F;
    std::function<complex(real,real,real)> u_boundary;
    
    finite_difference_2D();
    
    void init();
    void step();
    
    field & get_field();
    
    void set_field(const field &field);
  };
  
}


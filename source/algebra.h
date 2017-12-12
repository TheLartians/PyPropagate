#pragma once

#include <ostream>
#include <lars/ndarray.h>

namespace lars {
  namespace algebra{
    
    
    // we create our own complex class since std::complex seems to be extremely slow (3x slower than the class below on our hardware)
    
    template <class Real> struct Complex{
      using Scalar = Real;
      Scalar real,imag;
      
      Complex(Scalar r=0,Scalar i=0):real(r),imag(i){ }
      
      Complex(const Complex &) = default;
      Complex(Complex &&) = default;
      Complex &operator=(const Complex &) = default;
      Complex &operator=(Complex &&) = default;
      
      bool operator==(const Complex &other)const{ return real == other.real && imag == other.imag; }
      bool operator!=(const Complex &other)const{ return !(*this == other); }
      
      Complex &operator=(const Scalar &s){ real = s; imag = 0; return *this; }
      Complex &operator+=(const Complex &other){ real += other.real; imag += other.imag; return *this; }
      Complex &operator-=(const Complex &other){ real -= other.real; imag -= other.imag; return *this; }
      
      Complex operator+(const Complex &other)const{ return Complex(real + other.real, imag + other.imag ); }
      Complex operator-(const Complex &other)const{ return Complex(real - other.real, imag - other.imag ); }
      
      Complex operator*(const Scalar &other)const{ return Complex(real * other, imag * other ); }
      Complex operator/(const Scalar &other)const{ return Complex(real/other, imag/other ); }
      Complex operator*(const Complex &other)const{ return Complex(real * other.real - imag * other.imag,real * other.imag + imag * other.real ); }
      Complex operator/(const Complex &other)const{
        auto d = other.real * other.real + other.imag * other.imag;
        return Complex((real * other.real + imag * other.imag)/d,(-real * other.imag + imag * other.real)/d);
      }
      
      Complex operator-()const{ return Complex(-real,-imag); }
    };
    
    template <class Real> Complex<Real> operator+(const Real &s,const Complex<Real> &c){ return c+s; }
    template <class Real> Complex<Real>  operator-(const Real &s,const Complex<Real> &c){ return Complex<Real>(s,0)-c; }
    template <class Real> Complex<Real>  operator*(const Real &s,const Complex<Real> &c){ return c*s; }
    template <class Real> Complex<Real>  operator/(const Real &s,const Complex<Real> &c){ return Complex<Real>(s,0)/c; }

    template <class Real> std::ostream & operator<<(std::ostream &stream,const Complex<Real> &c){
      stream << '(' << c.real << ',' << c.imag << ')';
      return stream;
    }

    
    template <class scalar> using complex = Complex<scalar>;
    template <class scalar> using vector = HeapNDArray<scalar, DynamicIndexTuple<1>>;
    template <class scalar> using matrix = HeapNDArray<scalar,DynamicIndexTuple<2>>;
    template <class scalar> using array = HeapNDArray<scalar,DynamicIndexTuple<2>>;
    
    template <class v1,class v2,class v3, class v4,class v5,class v6> void tridiagonal(const v1 &a,const v2 &b,const v3 &c,const v4 &r,v5 &u,v6 &gam){
      int j,n=u.size();
      using scalar = typename v5::Scalar;
      scalar bet;
      // assert(b[0]!=0.);
      u[0]=r[0]/(bet=b[0]);
      for(j=1;j<n;++j){
        gam[j]=c[j-1]/bet;
        bet=b[j]-a[j]*gam[j];
        /// assert(bet!=0.);
        u[j]=(r[j]-a[j]*u[j-1])/bet;
      }
      for(j=n-2;j>=0;--j){
        u[j]-=gam[j+1]*u[j+1];
      }
    }
    
    template <class v1,class v2,class v3, class v4,class v5> void tridiagonal(const v1 &a,const v2 &b,const v3 &c,const v4 &r,v5 &u){
      v5 gam(a.size());
      tridiagonal(a,b,c,r,u,gam);
    }

    
    
    template <class A,class B,class R> void square_matrix_multiplication(const A &a,const B &b,R &r,size_t n){
      for(size_t j = 0; j < n; ++j){
        for(size_t i = 0; i < n; ++i){
          r(i,j) = 0;
          for(size_t k = 0; k < n; ++k){
            r(i,j) += a(i,k) * b(k,j);
          }
        }
      }
    };
    
  }
}







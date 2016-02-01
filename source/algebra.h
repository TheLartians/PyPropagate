#pragma once

#include <complex>
#include <Eigen/Dense>

namespace lars {
  namespace algebra{
  
    const int dynamic_size = Eigen::Dynamic;
    
    template <class real> using complex = std::complex<real>;
    template <class scalar,int dimension = dynamic_size> using vector = Eigen::Matrix< scalar, dimension, 1 >;
    template <class scalar,int m = dynamic_size,int n = dynamic_size> using matrix = Eigen::Matrix< scalar, m, n >;
    template <class scalar,int m = dynamic_size,int n = dynamic_size> using array = Eigen::Array< scalar, m, n >;
    
    template <class v1,class v2,class v3, class v4,class v5,class v6> void tridiagonal(const v1 &a,const v2 &b,const v3 &c,const v4 &r,v5 &u,v6 &gam){
      int j,n=a.size();
      using scalar = typename v5::Scalar;
      scalar bet;
      //ASSERT(b[0]!=0.);
      u[0]=r[0]/(bet=b[0]);
      for(j=1;j<n;++j){
        gam[j]=c[j-1]/bet;
        bet=b[j]-a[j]*gam[j];
        //ASSERT(bet!=0.);
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

    
  }
}

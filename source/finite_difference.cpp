/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//    Finite differences c++ solver for parabolic differential equations
//
//
//    File created in 2012 by Lars Melchior for the Institute for X-Ray Physics, Göttingen
//    
//

#include "finite_difference.h"
#include <lars/parallel.h>

namespace lars {
  
  //
  // ----------------------- Finite differences AF -----------------------
  //
  
  void finite_difference_aF::resize(size_t N){
    rf.resize(N);
    rfp.resize(N);
    u.resize(N);
    up.resize(N);
    
    B.resize(N-2);
    R.resize(N-2);
    tmp.resize(N-2);
  }
  
  void finite_difference_aF::step(){
    unsigned n = u.size()-2;
    finite_differences::pseudo_field A(-ra/2.);
    
    for (unsigned i=1; i<=n; ++i) {
      B[i-1] = 1.+ra-rf[i];
      R[i-1] = (up[1+i]+up[i-1])*ra/2.+up[i]*(1.+rfp[i]-ra);
    }
    
    R[0]   += u[0] * ra/2.;
    R[n-1] += u[n+1] * ra/2.;
        
    auto us = u.slice(static_index_tuple<1>(),make_dynamic_index_tuple(n));
    algebra::tridiagonal(A,B,A,R,us,tmp);
  }
  
  void finite_difference_aF::update(){
    std::swap(rf, rfp);
    std::swap(u, up);
  }
  
  //
  // ----------------------- Finite differences ACF -----------------------
  //
  
  void finite_difference_acF::resize(size_t Nx,size_t Ny){
    rf.resize(Nx,Ny);
    rfp.resize(Nx,Ny);
    u.resize(Nx,Ny);
    up.resize(Nx,Ny);
  }
  
  using namespace finite_differences;
  
  template<typename field> void acF_step(const complex &ra, const complex &rc, const field &rf, const field &rfp, field &u, const field &up){
    
    unsigned nx = u.size()-2;
    unsigned ny = u[0].size()-2;
    
    finite_differences::pseudo_field A(-rc);
    
    struct parallel_data{
      finite_differences::array_1D B,R,tmp;
      parallel_data(unsigned s):B(s),R(s),tmp(s){}
    };
    
    unique_parallel_for(1, nx+1, [&](unsigned i,parallel_data &d){
      
      for (unsigned j=1; j<=ny; ++j) {
        d.B[j-1] = 1.+rc*2.-rf(i,j);
        d.R[j-1] = (up(1.+i,j)+up(i-1.,j))*ra+up(i,j)*(1.+rfp(i,j)-ra*2.);
      }
      
      d.R[0]    += u(i,0) * rc;
      d.R[ny-1] += u(i,ny+1) * rc;
      
      auto us = u[i].slice(static_index_tuple<1>(), make_dynamic_index_tuple(ny));
      algebra::tridiagonal(A,d.B,A,d.R,us,d.tmp);
    },parallel_data(ny));

  }
  
  void finite_difference_acF::step_1(){
    acF_step(ra,rc,rf,rfp,u,up);
  }
  
  void finite_difference_acF::update(){
    std::swap(rf, rfp);
    std::swap(u, up);
  }

  void finite_difference_acF::step_2(){
    auto u_transposed = u.transpose();
    acF_step(rc,ra,rf.transpose(),rfp.transpose(),u_transposed,up.transpose());
  }
    
}



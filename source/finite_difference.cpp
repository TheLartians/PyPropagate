/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//    Finite differences c++ solver for parabolic differential equations
//
//
//    File created in 2012 by Lars Melchior for the Institute for X-Ray Physics, Göttingen
//    
//

#include "finite_difference.h"

namespace lars {
  
  //
  // ----------------------- Finite differences AF -----------------------
  //
 
  void finite_difference_AF::resize(size_t N){
    rf.resize(N);
    ra.resize(N);
    rfp.resize(N);
    rap.resize(N);
    u.resize(N);
    up.resize(N);
    
    A.resize(N-2);
    B.resize(N-2);
    R.resize(N-2);
    tmp.resize(N-2);
  }
  
  void finite_difference_AF::step(){
    unsigned n = u.size()-2;
    
    for (unsigned i=1; i<=n; ++i) {
      A[i-1] = -ra[i]/2.;
      B[i-1] = 1.+ra[i]-rf[i];
      R[i-1] = (up[1+i]+up[i-1])*rap[i]/2.+up[i]*(1.+rfp[i]-rap[i]);
    }
    
    R[0]   += u[0] * ra[0]/2.;
    R[n-1] += u[n+1] * ra[n+1]/2.;
    
    auto us = u.slice(StaticIndexTuple<1>(),make_dynamic_index_tuple(n));
    algebra::tridiagonal(A,B,A,R,us,tmp);
  }
  
  void finite_difference_AF::update(){
    std::swap(ra, rap);
    std::swap(rf, rfp);
    std::swap(u, up);
  }
  
  //
  // ----------------------- Finite differences ACF -----------------------
  //
  
  void finite_difference_ACF::resize(size_t Nx,size_t Ny){
    ra.resize(Ny,Nx);
    rap.resize(Ny,Nx);
    rc.resize(Ny,Nx);
    rcp.resize(Ny,Nx);
    rf.resize(Ny,Nx);
    rfp.resize(Ny,Nx);
    u.resize(Ny,Nx);
    up.resize(Ny,Nx);
  }
  
  using namespace finite_differences;
  
  struct trig_parallel_data{
    finite_differences::array_1D A,B,R,tmp;
    trig_parallel_data(unsigned s):A(s),B(s),R(s),tmp(s){ }
  };
  

  template<typename field> void ACF_step(const field &ra, const field &rc, const field &rf, const field &rap, const field &rcp, const field &rfp, field &u, const field &up,unsigned thread_count){
    
    unsigned nx = u.size()-2;
    unsigned ny = u[0].size()-2;
    
    unique_parallel_for(1, nx+1, [&](unsigned i,trig_parallel_data &d){
      
      for (unsigned j=1; j<=ny; ++j) {
        d.A[j-1] = -ra[i][j];
        d.B[j-1] = 1.+ra[i][j]*2.-rf[i][j];
        d.R[j-1] = (up[1.+i][j]+up[i-1.][j])*rcp[i][j]+up[i][j]*(1.+rfp[i][j]-rcp[i][j]*2.);
      }
      
      d.R[0]    += u[i][0] * ra[i][0];
      d.R[ny-1] += u[i][ny+1] * ra[i][ny+1];
      
      auto us = u[i].slice(StaticIndexTuple<1>(), make_dynamic_index_tuple(ny));
      algebra::tridiagonal(d.A,d.B,d.A,d.R,us,d.tmp);
    },trig_parallel_data(ny),thread_count);

  }

/*
  template<typename field> void ACF_step(const field &ra, const field &rc, const field &rf, const field &rap, const field &rcp, const field &rfp, field &u, const field &up){
    
    unsigned nx = u[0].size()-2;
    unsigned ny = u.size()-2;
    
    unique_parallel_for(1, ny+1, [&](unsigned i,trig_parallel_data &d){
      
      for (unsigned j=1; j<=nx; ++j) {
        d.A[j-1] = -ra[i][j];
        d.B[j-1] = 1.+ra[i][j]*2.-rf[i][j];
        d.R[j-1] = (up[i][j+1]+up[i][j-1])*rcp[i][j]+up[i][j]*(1.+rfp[i][j]-rcp[i][j]*2.);
      }
      
      d.R[0]    += u[i][0] * ra[i][0];
      d.R[nx-1] += u[i][nx+1] * ra[i][nx+1];
      
      auto us = u[i].slice(StaticIndexTuple<1>(), make_dynamic_index_tuple(nx));
      algebra::tridiagonal(d.A,d.B,d.A,d.R,us,d.tmp);
    },trig_parallel_data(nx));

  }
*/  

  void finite_difference_ACF::step_1(){
    ACF_step(ra,rc,rf,rap,rcp,rfp,u,up,thread_count);
  }
  
  void finite_difference_ACF::update(){
    std::swap(rf, rfp);
    std::swap(ra, rap);
    std::swap(rc, rcp);
    std::swap(u, up);
  }

  void finite_difference_ACF::step_2(){
    auto u_transposed = u.transpose();
    ACF_step(rc.transpose(),ra.transpose(),rf.transpose(),rcp.transpose(),rap.transpose(),rfp.transpose(),u_transposed,up.transpose(),thread_count);
  }
  
  void finite_difference_A0F::resize(size_t Nx,size_t Ny){
    rf.resize(Ny,Nx);
    ra.resize(Ny,Nx);
    rfp.resize(Ny,Nx);
    rap.resize(Ny,Nx);
    u.resize(Ny,Nx);
    up.resize(Ny,Nx);
  }
  
  void finite_difference_A0F::update(){
    std::swap(ra, rap);
    std::swap(rf, rfp);
    std::swap(u, up);
  }

  void finite_difference_A0F::step(){
    unsigned nx = u[0].size()-2;
    unsigned ny = u.size();
    
    unique_parallel_for(0, ny, [&](unsigned i,trig_parallel_data &d){
      
      for (unsigned j=1; j<=nx; ++j) {
        d.A[j-1] = -ra[i][j]/2.;
        d.B[j-1] = 1.+ra[i][j]-rf[i][j];
        d.R[j-1] = (up[i][j+1]+up[i][j-1])*rap[i][j]/2.+up[i][j]*(1.+rfp[i][j]-rap[i][j]);
      }
    
      d.R[0]    += u[i][0]    * ra[i][0]/2.;
      d.R[nx-1] += u[i][nx+1] * ra[i][nx+1]/2.;
      
      auto us = u[i].slice(StaticIndexTuple<1>(),make_dynamic_index_tuple(nx));
      algebra::tridiagonal(d.A,d.B,d.A,d.R,us,d.tmp);
      
    },trig_parallel_data(nx),thread_count);

  }

  
  void finite_difference_ABC::resize(size_t Nx,size_t Ny){
    ra.resize(Ny,Nx);
    rb.resize(Ny,Nx);
    rc.resize(Ny,Nx);
    rz.resize(Ny,Nx);
    rap.resize(Ny,Nx);
    rbp.resize(Ny,Nx);
    rcp.resize(Ny,Nx);
    rzp.resize(Ny,Nx);
    u.resize(Ny,Nx);
    up.resize(Ny,Nx);
  }
  
  void finite_difference_ABC::update(){
    std::swap(ra, rap);
    std::swap(rb, rbp);
    std::swap(rc, rcp);
    std::swap(rz, rzp);
    std::swap(u, up);
  }

  void finite_difference_ABC::step(){
    unsigned nx = u[0].size()-2;
    unsigned ny = u.size();
  
    struct trig_parallel_data{
      finite_differences::array_1D A,B,C,R,tmp;
      trig_parallel_data(unsigned s):A(s),B(s),C(s),R(s),tmp(s){ }
    };
    
    unique_parallel_for(0, ny, [&](unsigned i,trig_parallel_data &d){
      
      for (unsigned j=1; j<=nx; ++j) {
        d.A[j-1] = rb[i][j]   - ra[i][j];
        d.B[j-1] = 2.*ra[i][j] + rz[i][j] - rc[i][j];
        d.C[j-1] = - rb[i][j] - ra[i][j];
        d.R[j-1] = (rap[i][j] + rbp[i][j]) * up[i][j+1] + (rap[i][j] - rbp[i][j]) * up[i][j-1] + (rcp[i][j] + rzp[i][j] - 2.*rap[i][j]) * up[i][j];
      }
    
      d.R[0]    += (-rb[i][1]   + ra[i][1]) * u[i][0]; //TODO u[0][j]    * ra[0][j]/2.;
      d.R[nx-1] += ( rb[i][nx] + ra[i][nx]) * u[i][nx+1]; //TODO u[nx+1][j] * ra[nx+1][j]/2.;
 
      auto us = u[i].slice(StaticIndexTuple<1>(),make_dynamic_index_tuple(nx));
      algebra::tridiagonal(d.A,d.B,d.C,d.R,us,d.tmp);
      
    },trig_parallel_data(nx),thread_count);

  }
}



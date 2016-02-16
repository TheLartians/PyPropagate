#include "crank_nicolson.h"
#include <lars/parallel.h>
#include <iostream>

namespace lars {
  
  void crank_nicolson_2D::resize(unsigned Nx,unsigned Ny){
    auto fields = std::vector<field2D*>{&ra ,&rc ,&rd ,&re ,&rf ,&u , &rap, &rcp ,&rdp ,&rep ,&rfp ,&up };
    for(auto * f:fields){ f->resize(Nx, Ny); f->fill(0); }
  }
  
  
  crank_nicolson_2D::scalar crank_nicolson_2D::v1(unsigned i,unsigned j){
    return (up(i,1.+j)-up(i,j-1.))*rep(i,j)/4.+up(1.+i,j)*(rap(i,j)+rdp(i,j)/2.)+up(i,j)*(1.+rfp(i,j)-rap(i,j)*2.)+up(i-1.,j)*(rap(i,j)-rdp(i,j)/2.);
  }
  
  void crank_nicolson_2D::step_1(){
    
    unsigned nx = u.rows()-2;
    unsigned ny = u.cols()-2;
    
    unique_parallel_for(1, nx+1, [&](unsigned i,parallel_data &d){
      
      for (unsigned j=1; j<=ny; ++j) {
        d.A[j-1] = re(i,j)/4.-rc(i,j);
        d.B[j-1] = 1.+rc(i,j)*2.-rf(i,j);
        d.C[j-1] = -rc(i,j)-re(i,j)/4.;
        d.R[j-1] = v1(i,j);
      }
      
      d.R[0]    -= u(i,0) * (re(i,0)/4.-rc(i,0));
      d.R[ny-1] -= u(i,ny+1) * (-rc(i,ny+1)-re(i,ny+1)/4.);
      
      algebra::tridiagonal(d.A,d.B,d.C,d.R,d.U,d.tmp);
      for (unsigned j=1; j<=ny; ++j) u(i,j) = d.U[j-1];
      
    },parallel_data(ny));
    
  }
  
  void crank_nicolson_2D::update(){
    auto c = std::vector<field2D*>{&ra ,&rc ,&rd ,&re ,&rf ,&u };
    auto p = std::vector<field2D*>{&rap,&rcp,&rdp,&rep,&rfp,&up};
    for ( auto i : range(c.size()) ) c[i]->swap(*p[i]);
  }
  
  crank_nicolson_2D::scalar crank_nicolson_2D::v2(unsigned i,unsigned j){
    return (up(i,1.+j)-up(i,j-1.))*rdp(i,j)/4.+up(i,1.+j)*(rcp(i,j)+rep(i,j)/2.)+up(i,j)*(1.+rfp(i,j)-rcp(i,j)*2.)+up(i,j-1.)*(rcp(i,j)-rep(i,j)/2.);
  }
  
  void crank_nicolson_2D::step_2(){
    
    unsigned nx = u.rows()-2;
    unsigned ny = u.cols()-2;
    
    unique_parallel_for(1, ny+1, [&](unsigned j,parallel_data &d){
      
      for (unsigned i=1; i<=nx; ++i) {
        d.A[i-1] = rd(i,j)/4.-ra(i,j);
        d.B[i-1] = 1.+ra(i,j)*2.-rf(i,j);
        d.C[i-1] = -ra(i,j)-rd(i,j)/4.;
        d.R[i-1] = v2(i,j);
      }
      
      d.R[0]    -= u(0,j) * (rd(0,j)/4.-ra(0,j));
      d.R[nx-1] -= u(nx+1,j) * (-ra(nx+1,j)-rd(nx+1,j)/4.);
      
      algebra::tridiagonal(d.A,d.B,d.C,d.R,d.U,d.tmp);
      for (unsigned i=1; i<=nx; ++i) u(i,j) = d.U[i-1];
      
    },parallel_data(nx));
    
  }

  /*
  void crank_nicolson_2D::step_2(){
    
    bool lazy_step_2 = true;
    
    if(lazy_step_2){
      auto tp = std::vector<field2D*>{&ra ,&rc ,&rd ,&re ,&rf ,&u, &rap ,&rcp,&rdp,&rep,&rfp,&up};

      auto s1 = std::vector<field2D*>{&ra , &rd , &rap , &rdp };
      auto s2 = std::vector<field2D*>{&rc , &re , &rcp , &rep };
      
      for ( auto i : range(s1.size()) ) s1[i]->swap(*s2[i]);
      for ( auto i : range(tp.size()) ) tp[i]->transposeInPlace();
      
      step_1();
      
      for ( auto i : range(tp.size()) ) tp[i]->transposeInPlace();
      for ( auto i : range(s1.size()) ) s1[i]->swap(*s2[i]);

      return;
    }
   }
   //*/

}













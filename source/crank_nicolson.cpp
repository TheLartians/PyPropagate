#include "crank_nicolson.h"
#include <lars/parallel.h>
#include <iostream>

namespace lars {
  
  void crank_nicolson_2D::resize(unsigned Nx,unsigned Ny){
    auto fields = std::vector<field2D*>{&ra ,&rb ,&rc ,&rd ,&re ,&rf ,&u , &rap ,&rbp ,&rcp ,&rdp ,&rep ,&rfp ,&up };
    for(auto * f:fields){ f->resize(Nx, Ny); }
  }
  
  crank_nicolson_2D::scalar crank_nicolson_2D::v1(unsigned i,unsigned j){
    return ((up(i,1+j)-up(i,j-1))*rep(i,j)+rcp(i,j)*(up(1+i,1+j)+up(i-1,j-1)-up(1+i,j-1)-up(i-1,1+j)))/4.+up(1+i,j)*(rap(i,j)+rdp(i,j)/2.)+up(i,j)*(1.+rfp(i,j)-rap(i,j)*2.)+up(i-1,j)*(rap(i,j)-rdp(i,j)/2.);
  }
  
  void crank_nicolson_2D::step_1(){
    
    struct parallel_data{
      field1D A,B,C,R,U,tmp;
      parallel_data(unsigned s):A(s),B(s),C(s),R(s),U(s),tmp(s){}
    };
    
    unsigned nx = u.rows()-2;
    unsigned ny = u.cols()-2;
    
    unique_parallel_for(1, nx+1, [&](unsigned i,parallel_data &d){
      
      for (unsigned j=1; j<=ny; ++j) {
        d.A[j-1] = re(i,j)/4.-rb(i,j);
        d.B[j-1] = 1.+2.*rb(i,j)-rf(i,j);
        d.C[j-1] = -rb(i,j)-re(i,j)/4.;
        d.R[j-1] = v1(i,j);
      }
      
      d.R[0]    += u(i,0) * rb(i,0);
      d.R[ny-1] += u(i,ny+1) * rb(i,ny+1);
      
      algebra::tridiagonal(d.A,d.B,d.C,d.R,d.U,d.tmp);
      for (unsigned j=1; j<=ny; ++j) u(i,j) = d.U[j-1];
      
    },parallel_data(ny));
    
  }
  
  void crank_nicolson_2D::update(){
    auto c = std::vector<field2D*>{&ra ,&rb ,&rc ,&rd ,&re ,&rf ,&u };
    auto p = std::vector<field2D*>{&rap,&rbp,&rcp,&rdp,&rep,&rfp,&up};
    for ( auto i : range(c.size()) ) c[i]->swap(*p[i]);
  }
  
  void crank_nicolson_2D::step_2(){
    
    // TODO: Implement step 2 the right way!
    bool lazy_step_2 = true;
    
    if(lazy_step_2){
      auto tp = std::vector<field2D*>{&ra ,&rb ,&rc ,&rd ,&re ,&rf ,&u, &rap,&rbp,&rcp,&rdp,&rep,&rfp,&up};

      auto s1 = std::vector<field2D*>{&ra , &rd , &rap , &rdp };
      auto s2 = std::vector<field2D*>{&rb , &re , &rbp , &rep };
      
      for ( auto i : range(s1.size()) ) s1[i]->swap(*s2[i]);
      for ( auto i : range(tp.size()) ) tp[i]->transposeInPlace();
      
      step_1();
      
      for ( auto i : range(tp.size()) ) tp[i]->transposeInPlace();
      for ( auto i : range(s1.size()) ) s1[i]->swap(*s2[i]);

      return;
    }
    
  }
  

}













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
#include <iostream>

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
  
  void finite_difference_acF::step_1(){
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
  
  void finite_difference_acF::update(){
    std::swap(rf, rfp);
    std::swap(u, up);
  }

  void finite_difference_acF::step_2(){    
    
    auto transpose = [](field &f){ field transposed = f.transpose(); f = std::move(transposed); };
    
    std::swap(ra, rc);
    transpose(rf);
    transpose(rfp);
    transpose(u);
    transpose(up);
    step_1();
    std::swap(ra, rc);
    transpose(rf);
    transpose(rfp);
    transpose(u);
    transpose(up);
  }
  
  //
  // ----------------------- Finite differences 2D -----------------------
  //
  
  
  finite_difference_2D::finite_difference_2D():Cx(Ax),Cy(Ay){
    ready=false;
    sx=sy=0;
  }
  
  void finite_difference_2D::init(){
    ready=false;

    sx = field1.size();
    sy = field1[0].size();
    
    dx=(xmax-xmin)/(sx-1);
    dy=(ymax-ymin)/(sy-1);
    
    if(field1.size() < 9) throw std::runtime_error("field has size smaller than 9");
    
    rx=A*dz/(2*dx*dx);
    ry=A*dz/(2*dy*dy);
    
    Ax.resize(sx-2);
    Ay.resize(sy-2);
    for(unsigned i=0;i<Ax.size();++i)Ax[i]=-rx;
    for(unsigned i=0;i<Ay.size();++i)Ay[i]=-ry;
    
    field2.resize(sx, sy);
    CField1.resize(sx-2,sy-2);
    CField2.resize(sx-2,sy-2);
    
    z-=dz;
    update();
    update();
    
    ready=true;
  }
  
  finite_differences::complex finite_difference_2D::C(int xi,int yi,int zi){
    return (zi==0) ? CField1(xi-1,yi-1) : CField2(xi-1,yi-1);
  }
  
  finite_differences::complex &finite_difference_2D::u(int xi,int yi,int zi){
    return (zi==0) ? field1(xi,yi) : field2(xi,yi);
  }
  
  void finite_difference_2D::step(){
    if(!ready) init();
    
    struct parallel_data{
      array_1D B,D,U,tmp;
      parallel_data(unsigned s):B(s),D(s),U(s),tmp(s){}
    };
    
    unique_parallel_for(1, sy-1, [&](unsigned yi,parallel_data &d){
      for (unsigned xi=1; xi+1<sx; ++xi) {
        d.B[xi-1]=1.+2.*rx-C(xi,yi,1);
        d.D[xi-1]=(1.-2.*ry+C(xi,yi,0))*u(xi,yi,0)+ry*(u(xi,yi-1,0)+u(xi,yi+1,0));
      }
      
      d.D[0]+=rx*u(0,yi,1);
      d.D[sx-3]+=rx*u(sx-1,yi,1);
      
      algebra::tridiagonal(Ax,d.B,Cx,d.D,d.U,d.tmp);
      for (unsigned xi=1; xi<=sx-2; ++xi) u(xi,yi,1)=d.U[xi-1];
      
    },parallel_data(sx-2));
    
    update();
    
    unique_parallel_for(1, sx-1, [&](unsigned xi,parallel_data &d){
      for (unsigned yi=1; yi+1<sy; ++yi) {
        d.B[yi-1]=1.+2.*ry-C(xi,yi,1);
        d.D[yi-1]=(1.-2.*rx+C(xi,yi,0))*u(xi,yi,0)+rx*(u(xi-1,yi,0)+u(xi+1,yi,0));
      }
      
      d.D[0]+=ry*u(xi,0,1);
      d.D[sy-3]+=ry*u(xi,sy-1,1);
      
      algebra::tridiagonal(Ay,d.B,Cy,d.D,d.U,d.tmp);
      for (unsigned yi=1; yi<=sy-2; ++yi) u(xi,yi,1)=d.U[yi-1];
    },parallel_data(sy-2));
    
    update();
    
  }
  
  void finite_difference_2D::update(){
    real zn=z+dz;
    z+=dz/2.;
    
    if(!constant_F || !ready){
      CField1.swap(CField2);
      parallel_for(1, sx-1, [&](unsigned i){ for(unsigned j=1;j+1<sy;++j){ CField2(i-1,j-1)=F(xmin+i*dx,ymin+j*dy,zn)*dz/4.; } });
    }
    
    field1.swap(field2);
    
    auto f1 = std::async(std::launch::async,[&](){ for(unsigned i=1;i+1<sx;++i)field2(i , 0  )=u_boundary(xmin+i*dx ,ymin,zn); });
    auto f2 = std::async(std::launch::async,[&](){ for(unsigned i=1;i+1<sx;++i)field2(i ,sy-1)=u_boundary(xmin+i*dx ,ymax,zn); });
    auto f3 = std::async(std::launch::async,[&](){ for(unsigned i=1;i+1<sy;++i)field2( 0  ,i )=u_boundary(xmin, ymin+i*dy,zn); });
    auto f4 = std::async(std::launch::async,[&](){ for(unsigned i=1;i+1<sy;++i)field2(sx-1,i )=u_boundary(xmax, ymin+i*dy,zn); });
    
    f1.wait(); f2.wait(); f3.wait(); f4.wait();
  }
  
  finite_difference_2D::field & finite_difference_2D::get_field(){
    return field1;
  }
  
  void finite_difference_2D::set_field(const field &field){
    field1 = field;
    ready = false;
  }
  
}



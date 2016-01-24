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
    // ----------------------- Finite differences 1D -----------------------
    //
    
    finite_difference_1D::finite_difference_1D(){
      ready=false;
    }
    
    void finite_difference_1D::init(){
      if(field.size()==0) throw std::runtime_error("field has size 0");
      rx=A*dz/(dx*dx);
      s=field.size();
      xb=-1.*s*dx/2+dx/2;
      Ax.resize(s);
      for(int i=0;i<s;++i)Ax[i]=-rx/2.;
      Bx.resize(s);
      Dx.resize(s);
      tmp.resize(s);
      
      for(int i=0;i<s;++i){
        real x=xb+dx*i;
        field[i] = u0(x,z);
      }

      ready=true;
    }
    
    void finite_difference_1D::step(){
      if(!ready) init();
      real zh=z+dz/2;
      real zn=z+dz;
      
      for(int i=0;i<s;++i){
        real x=xb+dx*i;
        complex C=F(x,zh)*dz/2.;
        Bx[i]=1.+rx-C;
        Dx[i]=(1.-rx+C)*field[i];
        if(i>0)Dx[i]+=rx/2.*field[i-1];
        else Dx[i]+=rx/2.*(u0_boundary(x-dx,z)+u0_boundary(x-dx,zn));
        if(i<s-1)Dx[i]+=rx/2.*field[i+1];
        else Dx[i]+=rx/2.*(u0_boundary(x+dx,z)+u0_boundary(x+dx,zn));
      }
      
      algebra::tridiagonal(Ax,Bx,Ax,Dx,field,tmp);
      
      z+=dz;
    }
    
    //
    // ----------------------- Finite differences 2D -----------------------
    //
  
  
    finite_difference_2D::finite_difference_2D():Cx(Ax),Cy(Ay){
      ready=false;
      sx=sy=0;
    }
  
    void finite_difference_2D::init(){
      rx=A*dz/(dx*dx);
      ry=A*dz/(dy*dy);
      xb=-sx*dx/2+dx*0.5;
      yb=-sy*dy/2+dy*0.5;
      
      Ax.resize(sx);
      Ay.resize(sy);
      for(int i=0;i<sx;++i)Ax[i]=-rx/2.;
      for(int i=0;i<sy;++i)Ay[i]=-ry/2.;
      
      CField1.resize(sx,sy);
      CField2.resize(sx,sy);
      
      z-=dz;
      update();
      update();
      
      for(int i=0;i<sx;++i)for(int j=0;j<sy;++j){ field1(i+1,j+1)=u0(xb+i*dx,yb+j*dy,z); }
      
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
        vector B,D,U,tmp;
        parallel_data(unsigned s):B(s),D(s),U(s),tmp(s){}
      };
      
      unique_parallel_for(1, sy+1, [&](unsigned yi,parallel_data &d){
        for (int xi=1; xi<=sx; ++xi) {
          d.B[xi-1]=1.+rx-C(xi,yi,1);
          d.D[xi-1]=(1.-ry+C(xi,yi,0))*u(xi,yi,0)+ry/2.*(u(xi,yi-1,0)+u(xi,yi+1,0));
        }

        d.D[0]+=rx/2.*u(0,yi,1);
        d.D[sx-1]+=rx/2.*u(sx+1,yi,1);
        
        algebra::tridiagonal(Ax,d.B,Cx,d.D,d.U,d.tmp);
        for (int xi=1; xi<=sx; ++xi) u(xi,yi,1)=d.U[xi-1];
      },parallel_data(sx));
      
      update();
      
      unique_parallel_for(1, sx+1, [&](unsigned xi,parallel_data &d){
        for (int yi=1; yi<=sy; ++yi) {
          d.B[yi-1]=1.+ry-C(xi,yi,1);
          d.D[yi-1]=(1.-rx+C(xi,yi,0))*u(xi,yi,0)+rx/2.*(u(xi-1,yi,0)+u(xi+1,yi,0));
        }
        
        d.D[0]+=ry/2.*u(xi,0,1);
        d.D[sy-1]+=ry/2.*u(xi,sy+1,1);
        
        algebra::tridiagonal(Ay,d.B,Cy,d.D,d.U,d.tmp);
        for (int yi=1; yi<=sy; ++yi) u(xi,yi,1)=d.U[yi-1];
      },parallel_data(sy));
      
      update();
    
    }
        
    void finite_difference_2D::update(){
      real zn=z+dz;
      z+=dz/2.;
      
      CField1.swap(CField2);
        
      parallel_for(0, sx, [&](unsigned i){ for(int j=0;j<sy;++j){ CField2(i,j)=F(xb+i*dx,yb+j*dy,zn)*dz/4.; } });
      
      field1.swap(field2);

      auto f1 = std::async(std::launch::async,[&](){ for(int i=0;i<sx;++i)field2(i+1 , 0  )=u0_boundary(xb+i*dx ,yb-dy,zn   ); });
      auto f2 = std::async(std::launch::async,[&](){ for(int i=0;i<sx;++i)field2(i+1 ,sy+1)=u0_boundary(xb+i*dx ,yb+dy*sy,zn); });
      auto f3 = std::async(std::launch::async,[&](){ for(int i=0;i<sy;++i)field2( 0  ,i+1 )=u0_boundary(xb-dx   ,yb+i*dy,zn ); });
      auto f4 = std::async(std::launch::async,[&](){ for(int i=0;i<sy;++i)field2(sx+1,i+1 )=u0_boundary(xb+dx*sx,yb+i*dy,zn ); });
      
      f1.wait(); f2.wait(); f3.wait(); f4.wait();
    }
  
    Eigen::Block<finite_difference_2D::field> finite_difference_2D::get_field(){
      return field1.block(1, 1, sx, sy);
    }
  
    finite_difference_2D::field & finite_difference_2D::get_full_field(){
      return field1;
    }

  
    void finite_difference_2D::set_field(const field &field){
      field1 = field;
    }
  
    void finite_difference_2D::resize(int nsx,int nsy){
      sx=nsx;
      sy=nsy;
      field1.resize(sx+2,sy+2);
      field2.resize(sx+2,sy+2);
      field1(0,0)=field1(sx+1,sy+1)=field1(0,sy+1)=field1(sx+1,0)=0;
      field2(0,0)=field2(sx+1,sy+1)=field2(0,sy+1)=field2(sx+1,0)=0;
    }
}

 

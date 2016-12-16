#pragma once

#include <array>
#include <stdexcept>
#include <limits>
#include <cmath>

namespace lars {
    
    template <class V> class Rectangle{
    public:
        
        using Vector = V;
        using Scalar = typename Vector::Scalar;
        
        struct const_iterator{
            const Rectangle * parent;
            unsigned char idx;
            const_iterator(const Rectangle * p,unsigned char i):parent(p),idx(i){}
            void operator++(){ idx++; }
            Vector operator*()const{
                if(idx==0)return parent->center() + ( - parent->extent()[0] - parent->extent()[1] ) / 2 ;
                if(idx==1)return parent->center() + (   parent->extent()[0] - parent->extent()[1] ) / 2 ;
                if(idx==2)return parent->center() + (   parent->extent()[0] + parent->extent()[1] ) / 2 ;
                if(idx==3)return parent->center() + ( - parent->extent()[0] + parent->extent()[1] ) / 2 ;
                throw std::runtime_error("iterator overflow");
            }
            bool operator!=(const const_iterator &other)const{ return other.idx != idx || other.parent != parent; }
        };
        
    private:
        Vector pcenter;
        std::array<Vector,2> pextent;
        
    public:
        
        Rectangle():pcenter(Vector::create_zeros()),pextent({{Vector::create_zeros(),Vector::create_zeros()}}){}
        Rectangle(Scalar width,Scalar height):pcenter(Vector::create(0,0)),pextent({{Vector::create(width,0),Vector::create(0,height)}}){ }
        Rectangle(const Vector &c,Scalar width,Scalar height):pcenter(c),pextent({{Vector::create(width,0),Vector::create(0,height)}}){ }
        Rectangle(const Vector &ll,const Vector &ur):pcenter((ll + ur)/2),pextent({{Vector::create((ur.x()-ll.x()),0),Vector::create(0,(ur.y()-ll.y()))}}){ }
        Rectangle(Scalar x1,Scalar y1,Scalar x2,Scalar y2):pcenter(Vector::create((x1+x2)/2,(y1+y2)/2)),pextent({{Vector::create((x2-x1),0),Vector::create(0,(y2-y1))}}){ }
        Rectangle(Vector c,Vector extent_1,Vector extent_2):pcenter(c),pextent({{extent_1,extent_2}}){ }
        
        const Vector &center()const{ return pcenter; }
        Vector &center(){ return pcenter; }
        const std::array<Vector,2> &extent()const{ return pextent; }
        std::array<Vector,2> &extent(){ return pextent; }
      
        Scalar angle()const{ return extent()[0].angle(); }
        bool is_aligned()const{ return extent()[0](1) == 0 && extent()[1](0) == 0; }
      
        Scalar width() const{ return extent()[0].norm(); }
        Scalar height()const{ return extent()[1].norm(); }
        
        Vector lower_left()const { return center() + ( - extent()[0] - extent()[1] ) / 2 ; }
        Vector upper_right()const{ return center() + ( extent()[0] + extent()[1] ) / 2 ; }
        Vector upper_left()const { return center() + ( - extent()[0] + extent()[1] ) / 2 ; }
        Vector lower_right()const{ return center() + (  extent()[0] - extent()[1] ) / 2 ; }
      
        const_iterator begin()const{ return const_iterator(this,0); }
        const_iterator end()const{ return const_iterator(this,4); }
        
        template <class F> void for_all_edges(F f)const{
            f(upper_left());
            f(upper_right());
            f(lower_left());
            f(lower_right());
        }
        
    };
    
    template <class V> class AlignedRectangle{
    public:
        
        using Vector = V;
        using Scalar = typename Vector::Scalar;
        
        struct const_iterator{
            const AlignedRectangle & parent;
            unsigned char idx;
            const_iterator(const AlignedRectangle & p,unsigned char i):parent(p),idx(i){}
            const_iterator operator++(){ return const_iterator(parent,idx++); }
            Vector operator*()const{
                if(idx == 0) return parent.lower_left();
                if(idx == 1) return parent.lower_right();
                if(idx == 2) return parent.upper_right();
                if(idx == 3) return parent.upper_left();
                throw std::runtime_error("iterator overflow");
            }
            bool operator!=(const const_iterator &other)const{ return other.idx != idx || other.parent != parent; }
        };
        
    private:
        Vector ll,ur;
        
    public:
        
        AlignedRectangle():ll(Vector::create(std::numeric_limits<Scalar>::max(),std::numeric_limits<Scalar>::max())),ur(Vector::create(std::numeric_limits<Scalar>::lowest(),std::numeric_limits<Scalar>::lowest())){}
        AlignedRectangle(Scalar width,Scalar height):ll(Vector::create(-width/2,-height/2)),ur(Vector::create(width/2,height/2)){ }
        AlignedRectangle(const Vector &c,Scalar width,Scalar height):ll(c - Vector::create(width/2,height/2)),ur(c+Vector::create(width/2,height/2)){ }
        AlignedRectangle(const Vector &ll,const Vector &ur):ll(ll),ur(ur){ }
        AlignedRectangle(Scalar x1,Scalar y1,Scalar x2,Scalar y2):ll(Vector::create(x1,y1)),ur(Vector::create(x2,y2)){ }
        
        bool is_valid()const{ return xmax() >= xmin(); }
        
        Scalar &xmin(){ return ll(0);  }
        Scalar &xmax(){ return ur(0); }
        Scalar &ymin(){ return ll(1);  }
        Scalar &ymax(){ return ur(1); }
        
        const Scalar &xmin()const{ return ll(0);  }
        const Scalar &xmax()const{ return ur(0); }
        const Scalar &ymin()const{ return ll(1);  }
        const Scalar &ymax()const{ return ur(1); }
        
        Scalar angle() const{ return 0; }
        Scalar width() const{ return fabs(ur(0) - ll(0)); }
        Scalar height()const{ return fabs(ur(1) - ll(1)); }
        Vector center()const{ return (ll + ur)/2; }
        
        AlignedRectangle operator+(const Vector &v)const  { return AlignedRectangle(ll+v,ur+v); }
        AlignedRectangle operator*(const Scalar &v)const  { return AlignedRectangle(ll*v,ur*v); }
        AlignedRectangle operator*(const Vector &v)const  { return AlignedRectangle(ll.as_array()*v,ur.as_array()*v); }
        AlignedRectangle &operator+=(const Vector &v){ ll+=v; ur+=v; return *this; }
        AlignedRectangle &operator*=(const Scalar &v){ ll*=v; ur*=v; return *this; }
        AlignedRectangle &operator*=(const Vector &v){ ll.as_array()*=v; ur.as_array()*=v; return *this; }
      
        const Vector &lower_left()const{ return ll; }
        const Vector &upper_right()const{ return ur; }
        Vector &lower_left(){ return ll; }
        Vector &upper_right(){ return ur; }
        Vector upper_left()const{ return Vector::create(xmin(),ymax()); }
        Vector lower_right()const{ return Vector::create(xmax(),ymin()); }
        
        bool operator==(const AlignedRectangle &other)const{ return ll == other.ll && ur == other.ur; }
        bool operator!=(const AlignedRectangle &other)const{ return ll != other.ll || ur != other.ur; }
        
        const_iterator begin()const{ return const_iterator(*this,0); }
        const_iterator end()const{ return const_iterator(*this,4); }
        
        operator Rectangle<Vector>()const{ return Rectangle<Vector>(lower_left(), upper_right()); }
        
        template <class F> void for_all_edges(F f)const{
            f(upper_left());
            f(upper_right());
            f(lower_left());
            f(lower_right());
        }
        
    };
    
    template<class Vector> bool convex_polygons_intersect(const Vector &p,const AlignedRectangle<Vector> &rect){
      return p(0) >= rect.xmin() && p(0) <= rect.xmax() && p(1) >= rect.ymin() && p(1) <= rect.ymax();
    }
    
    template<class Vector> bool convex_polygons_intersect(const AlignedRectangle<Vector> &rect,const Vector &p){
        return convex_polygons_intersect(p,rect);
    }
    
    template<class Vector> class BoundingBoxCreator{
      AlignedRectangle<Vector> bb;
      using Scalar = typename Vector::Scalar;
        
    public:
      
      BoundingBoxCreator() = default;
      
      void add_point(Scalar x, Scalar y){
          if(x<bb.xmin()) bb.xmin() = x;
          if(y<bb.ymin()) bb.ymin() = y;
          if(x>bb.xmax()) bb.xmax() = x;
          if(y>bb.ymax()) bb.ymax() = y;
      }
        
      void reset(){ bb = AlignedRectangle<Vector>();  }
      void add_point(const Vector &v){ add_point(v.x(), v.y()); }
      template <class Polygon> void add_polygon(const Polygon &v){ for(auto p:v) add_point(p); }
      const AlignedRectangle<Vector> &get_bounding_box()const{ return bb; }
      
      BoundingBoxCreator & operator+=(const Vector &v){ bb += v; }
      BoundingBoxCreator & operator-=(const Vector &v){ bb -= v; }
      
      template <class Polygon> static AlignedRectangle<Vector> get_bounding_box(Polygon p){
        BoundingBoxCreator bbc;
        bbc.add_polygon(p);
        return bbc.get_bounding_box();
      }

    };
  
  template<typename Char, typename Traits,class Vector> std::basic_ostream<Char,Traits> &operator<<(std::basic_ostream<Char,Traits> &stream, const Rectangle<Vector> &rect){
    stream << "Rectange("
     << rect.lower_left() << ','
     << rect.lower_right() << ','
     << rect.upper_left() << ','
     << rect.upper_right() << ')';
    return stream;
  }
  
  template<typename Char, typename Traits,class Vector> std::basic_ostream<Char,Traits> &operator<<(std::basic_ostream<Char,Traits> &stream, const AlignedRectangle<Vector> &rect){
    stream << "AlignedRectangle(" << rect.xmin() << ',' << rect.ymin() << ',' << rect.xmax() << ',' << rect.ymax() << ')';
    return stream;
  }

  
}













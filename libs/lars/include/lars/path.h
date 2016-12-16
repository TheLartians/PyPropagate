#pragma once

#include <memory>
#include <vector>
#include <functional>
#include <stdexcept>
#include <unordered_map>
#include <math.h> // std::nextafter
#include <limits>
#include <stdexcept>

#include <lars/hashers.h>
#include <lars/iterators.h>
#include <lars/rectangle.h>
#include <lars/math.h>

namespace lars {
  
  namespace curves{
    template <typename Vector> struct Line;
    template <typename Vector> struct Conic;
    template <typename Vector> struct Cubic;
  }
  
  template <typename Vector> class Path;
  
  template <typename Vector> class PathVisitor{
  public:
    virtual void visit(Path<Vector> *) = 0;
    virtual void visit(curves::Line<Vector> *) = 0;
    virtual void visit(curves::Conic<Vector> *) = 0;
    virtual void visit(curves::Cubic<Vector> *) = 0;
    virtual ~PathVisitor(){}
  };

  template <typename Vector> class ConstPathVisitor{
  public:
    virtual void visit(const Path<Vector> *) = 0;
    virtual void visit(const curves::Line<Vector> *) = 0;
    virtual void visit(const curves::Conic<Vector> *) = 0;
    virtual void visit(const curves::Cubic<Vector> *) = 0;
    virtual ~ConstPathVisitor(){}
  };

  template <typename Vector> struct Curve{
    virtual void accept(PathVisitor<Vector> * v) = 0;
    virtual void accept(ConstPathVisitor<Vector> * v) const = 0;
    virtual ~Curve(){}
  };
  
  namespace curves{
    template <typename Vector> struct Line:public Curve<Vector>{
      Vector to;
      Line(const Line &other) = default;
      Line(const Vector &to):to(to){}
      void accept(PathVisitor<Vector> * v) { v->visit(this); }
      void accept(ConstPathVisitor<Vector> * v)const{ v->visit(this); }
    };
    template <typename Vector> struct Conic:public Curve<Vector>{
      Vector to,control;
      Conic(const Conic &other) = default;
      Conic(const Vector &control,const Vector &to):to(to),control(control){}
      void accept(PathVisitor<Vector> * v) { v->visit(this); }
      void accept(ConstPathVisitor<Vector> * v)const{ v->visit(this); }
    };
    template <typename Vector> struct Cubic:public Curve<Vector>{
      Vector to,control_1,control_2;
      Cubic(const Cubic &other) = default;
      Cubic(const Vector &control_1,const Vector &control_2,const Vector &to):to(to),control_1(control_1),control_2(control_2){}
      void accept(PathVisitor<Vector> * v) { v->visit(this); }
      void accept(ConstPathVisitor<Vector> * v)const{ v->visit(this); }
    };
  }
  
  namespace PathVisitors{
    template <typename Vector> class CopyVisitor;
    template <typename Vector> class IncrementalSegmenter;
    template <typename Vector> class RecursiveSegmenter;
    template <typename Vector> class BoundingBoxVisitor;
    template <typename Vector> class MultiplyVisitor;
    template <typename Vector> class TranslateVisitor;
  }
    
  template <typename Vector> class Path:public std::vector<std::unique_ptr<Curve<Vector>>>{
  public:
    Vector start;
    
    Path(){}
    Path(const Vector &start):start(start){ }
    Path(const Path &other){ PathVisitors::CopyVisitor<Vector> v(this); v.visit(&other); }
    Path(Path &&other) = default;
    Path &operator=(const Path &other){ PathVisitors::CopyVisitor<Vector> v(this); v.visit(&other); return *this; }
    Path &operator=(Path &&other) = default;
    
    void accept(PathVisitor<Vector> *v){ v->visit(this); }
    void accept(ConstPathVisitor<Vector> *v)const{ v->visit(this); }
    
    template < class Curve,typename...Args> void add_point(Args... args){
      this->emplace_back(new Curve(args...));
    }
    
    std::vector<Vector> get_segmented_path_with_angle_tolerance(double angle_tolerance,double length_tolerance = 0)const{
      PathVisitors::RecursiveSegmenter<Vector> v(angle_tolerance,length_tolerance);
      v.visit(this);
      return v.get_points();
    }
    
    std::vector<Vector> get_segmented_path_with_subdivisions(double steps)const{
      PathVisitors::IncrementalSegmenter<Vector> v(steps);
      v.visit(this);
      return v.get_points();
    }
    
    AlignedRectangle<Vector> get_bounding_box()const{
      PathVisitors::BoundingBoxVisitor<Vector> v;
      v.visit(this);
      return v.get_bounding_box();
    }
    
    Path &operator*=(typename Vector::Scalar x){
      PathVisitors::MultiplyVisitor<Vector> v(x);
      v.visit(this);
      return *this;
    }
    
    Path &operator+=(Vector x){
      PathVisitors::TranslateVisitor<Vector> v(x);
      v.visit(this);
      return *this;
    }
    
  };
  
  
  namespace PathVisitors{
    
    template<typename Vector> class CopyVisitor:public ConstPathVisitor<Vector>{
      public:
      Path<Vector> * copy;
      CopyVisitor(Path<Vector> *c):copy(c){}
      virtual void visit(const Path<Vector> *p){ copy->start = p->start; for(auto &c:*p) c->accept(this); }
      virtual void visit(const curves::Line<Vector> *c){ copy->template add_point<curves::Line<Vector>>(*c); }
      virtual void visit(const curves::Conic<Vector> *c){ copy->template add_point<curves::Conic<Vector>>(*c); }
      virtual void visit(const curves::Cubic<Vector> *c){ copy->template add_point<curves::Cubic<Vector>>(*c); }
    };
    
    template <typename Vector> class IncrementalSegmenter:public ConstPathVisitor<Vector>{
      
      using real = typename Vector::Scalar;
      
      Vector intermediate_Line(Vector p0,Vector p1,real t){
        return (1-t)*p0+t*p1;
      }
      
      Vector intermediate_Conic(Vector p0,Vector p1,Vector p2,real t){
        return (1-t)*(1-t)*p0 + 2*t*(1-t)*p1 + t*t*p2;
      }
      
      Vector intermediate_Cubic(Vector p0,Vector p1,Vector p2,Vector p3,real t){
        return (1-t)*(1-t)*(1-t)*p0 + 3*(1-t)*(1-t)*t*p1 + 3*(1-t)*t*t*p2 + t*t*t*p3;
      }
      
      unsigned steps;
      
      std::vector<Vector> points;
      
      const Vector & current_point(){
        return points.back();
      }
      
      virtual void add_point(const Vector &v){
        points.push_back(v);
      }
      
    public:
      
      IncrementalSegmenter(unsigned steps):steps(steps){}
      
      void visit(const Path<Vector> * e){
        points.clear();
        add_point(e->start);
        for(auto &p:*e)p->accept(this);
      }
      
      void visit(const curves::Line<Vector> * e){
        add_point(e->to);
      }
      
      void visit(const curves::Conic<Vector> * e){
        auto from = current_point();
        for(auto i:range<unsigned>(1,steps+1)){
          real t = i/double(steps+1);
          add_point(intermediate_Conic(from,e->control,e->to,t));
        }
        add_point(e->to);
      }
      
      void visit(const curves::Cubic<Vector> * e){
        auto from = current_point();
        for(auto i:range<unsigned>(1,steps+1)){
          real t = i/double(steps+1);
          add_point(intermediate_Cubic(from,e->control_1,e->control_2,e->to,t));
        }
        add_point(e->to);
      }
      
      std::vector<Vector> get_points(){
        return points;
      }
      
    };
    
    
    template <typename Vector> class RecursiveSegmenter:public ConstPathVisitor<Vector>{
      
      std::vector<Vector> points;
      using scalar = typename Vector::Scalar;
      
      scalar angle_tolerance,cos_angle_tolerance_squared,cos_half_angle_tolerance_squared;
      scalar length_tolerance,length_tolerance_squared;
      
      const Vector & current_point(){ return points.back(); }
      virtual void add_point(const Vector &v){ points.push_back(v); }
      
      void draw_Conic_Curve(const Vector &A,const Vector &C,const Vector &B){
        Vector AC = (C - A);
        Vector BC = (C - B);
        
        scalar l1_squared = AC.norm_squared(), l2_squared = BC.norm_squared();
        scalar dot = (-AC).dot(BC);
        
        scalar cos_squared = dot*dot / (l1_squared * l2_squared);
        
        auto length_test = [&](){ return l1_squared < length_tolerance_squared && l2_squared < length_tolerance_squared; };
        auto angle_test  = [&](){ return sign(dot) * cos_squared > cos_angle_tolerance_squared; };
        
        if( length_test() || angle_test() ) add_point(B);
        else {
          Vector C1 = A + AC/2;
          Vector C2 = B + BC/2;
          Vector M  = C1 + (C2 - C1)/2;
          draw_Conic_Curve(A, C1, M);
          draw_Conic_Curve(M, C2, B);
        }
      }
      
      void draw_Cubic_Curve(const Vector &A,const Vector &C1,const Vector &C2,const Vector &B){
        
        Vector AC1 = C1 - A;
        Vector BC2 = C2 - B;
        Vector C1C2 = C2-C1;
        
        scalar l1_squared = AC1.norm_squared(), l2_squared = BC2.norm_squared(), l3_squared = C1C2.norm_squared();
        
        scalar dot1 = AC1.dot(C1C2);
        scalar dot2 = (-C1C2).dot(BC2);
        
        scalar cos1_squared = dot1*dot1 / (l1_squared * l3_squared);
        scalar cos2_squared = dot2*dot2 / (l2_squared * l3_squared);
        
        auto length_test = [&](){ return l1_squared < length_tolerance_squared && l2_squared < length_tolerance_squared && l3_squared < length_tolerance_squared; };
        auto angle_test  = [&](){ return sign(dot1) * cos1_squared > cos_half_angle_tolerance_squared && sign(dot2) * cos2_squared > cos_half_angle_tolerance_squared; };
        
        if( length_test() || angle_test() ) add_point(B);
        else {
          Vector C11 = A + AC1/2;
          Vector C22 = B + BC2/2;
          Vector CM = C1 + C1C2/2;
          Vector C12 = C11 + (CM - C11)/2;
          Vector C21 = C22 + (CM - C22)/2;
          
          Vector M = C12 + (C21 - C12)/2;

          draw_Cubic_Curve(A, C11, C12, M);
          draw_Cubic_Curve(M, C21, C22, B);
        }
      }
      
    public:
      
      RecursiveSegmenter(scalar angle_tolerance,scalar length_tolerance):angle_tolerance(angle_tolerance),length_tolerance(length_tolerance){
        if(fabs(angle_tolerance)  >= M_PI/2) throw std::domain_error("path segmentation angle tolerance must be smaller than pi/2");
        if(length_tolerance <  0) throw std::domain_error("path segmentation length tolerance must be larger than zero");
        cos_angle_tolerance_squared = cos(angle_tolerance)*cos(angle_tolerance);
        cos_half_angle_tolerance_squared = cos(angle_tolerance/2)*cos(angle_tolerance/2);
        length_tolerance_squared = length_tolerance * length_tolerance;
      }
      
      void visit(const Path<Vector> * e){
        points.clear();
        add_point(e->start);
        for(auto &p:*e)p->accept(this);
      }
      
      void visit(const curves::Line<Vector> * e){
        add_point(e->to);
      }
      
      void visit(const curves::Conic<Vector> * e){
        draw_Conic_Curve(current_point(),e->control,e->to);
      }
      
      void visit(const curves::Cubic<Vector> * e){
        draw_Cubic_Curve(current_point(),e->control_1,e->control_2,e->to);
      }
      
      std::vector<Vector> get_points(){
        return points;
      }
      
    };
    
    template <typename Vector> class BoundingBoxVisitor:public ConstPathVisitor<Vector>,public BoundingBoxCreator<Vector>{
    public:
      
      void visit(const Path<Vector> * e){
        BoundingBoxCreator<Vector>::add_point(e->start);
        for(auto &p:*e)p->accept(this);
      }
      
      void visit(const curves::Line<Vector> * e){
        BoundingBoxCreator<Vector>::add_point(e->to);
      }
      
      void visit(const curves::Conic<Vector> * e){
        BoundingBoxCreator<Vector>::add_point(e->to);
        BoundingBoxCreator<Vector>::add_point(e->control);
      }
      
      void visit(const curves::Cubic<Vector> * e){
        BoundingBoxCreator<Vector>::add_point(e->to);
        BoundingBoxCreator<Vector>::add_point(e->control_1);
        BoundingBoxCreator<Vector>::add_point(e->control_2);
      }
      
    };
    
    template <typename Vector> class MultiplyVisitor:public PathVisitor<Vector>{
      typename Vector::Scalar value;
      
    public:
      
      MultiplyVisitor(typename Vector::Scalar value):value(value){}
      
      void visit(Path<Vector> * e){
        e->start *= value;
        for(auto &p:*e)p->accept(this);
      }
      
      void visit(curves::Line<Vector> * e){
        e->to *= value;
      }
      
      void visit(curves::Conic<Vector> * e){
        e->to *= value;
        e->control *= value;
      }
      
      void visit(curves::Cubic<Vector> * e){
        e->to *= value;
        e->control_1 *= value;
        e->control_2 *= value;
      }
      
    };
    
    template <typename Vector> class TranslateVisitor:public PathVisitor<Vector>{
      Vector value;
      
    public:
      
      TranslateVisitor(Vector value):value(value){}
      
      void visit(Path<Vector> * e){
        e->start += value;
        for(auto &p:*e)p->accept(this);
      }
      
      void visit(curves::Line<Vector> * e){
        e->to += value;
      }
      
      void visit(curves::Conic<Vector> * e){
        e->to += value;
        e->control += value;
      }
      
      void visit(curves::Cubic<Vector> * e){
        e->to += value;
        e->control_1 += value;
        e->control_2 += value;
      }
      
    };
    
    template <typename Vector> class TouchingPointResolver:public PathVisitor<Vector>{
      std::unordered_map<Vector,const Vector *,ArrayHasher> points;
      
      const Vector * last_point = nullptr;
      
      void add_point(Vector &p){
        if(points.find(p) != points.end()){
          Vector old = p;
          auto & p2 = *points.find(p)->second;
          for(auto i:range(p.size())){
            if( (*last_point)[i] > p2[i] ) p[i] = nextafter(p[i],std::numeric_limits<typename Vector::Scalar>::max());
            else if((*last_point)[i] < p2[i]) p[i] = nextafter(p[i],std::numeric_limits<typename Vector::Scalar>::min());
          }
          if(p == old) throw std::runtime_error("cannot move point") ;
          add_point(p);
        }
        else {
          points[p] = last_point;
          last_point = &p;
        }
      }
      
    public:
      
      void visit(Path<Vector> * e){
        last_point = &e->start;
        for(auto &p:*e)p->accept(this);
      }
      
      void visit(curves::Line<Vector> * e){
        add_point(e->to);
      }
      
      void visit(curves::Conic<Vector> * e){
        add_point(e->to);
      }
      
      void visit(curves::Cubic<Vector> * e){
        add_point(e->to);
      }
      
    };
    
  }
  
}

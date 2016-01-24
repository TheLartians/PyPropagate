


#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <boost/python.hpp>
#include <numpy/arrayobject.h>

#include <stdexcept>
#include <array>
#include <complex>
#include <iostream>
using namespace std;

#include <Eigen/Core>
using namespace Eigen;

#include <lars/iterators.h>
#include <lars/parallel.h>
using namespace lars;

#include "finite_difference.h"

namespace python{
  using scalar = std::complex<double>;
  
  using vector = Matrix<scalar,Dynamic,1>;
  using vector_map = Map< Matrix<scalar,Dynamic,1> >;
  
  using matrix = Matrix<scalar,Dynamic,Dynamic>;
  using matrix_map = Map< Matrix<scalar,Dynamic,Dynamic> >;
  
  const auto py_scalar = NPY_COMPLEX128;
  
  bool initialized = false;

  void init(){
    if(initialized) return;
    import_array();
    initialized = true;
  }
  
  vector_map convert_to_vector(PyObject * Arr){
    init();
    PyArrayObject * converted = (PyArrayObject *)PyArray_FromAny(Arr,PyArray_DescrFromType(py_scalar),1,1,NPY_ARRAY_CARRAY,NULL);
    if(!converted) throw std::invalid_argument("cannot convert python object to vector");
    vector_map M((scalar *) PyArray_DATA(converted),PyArray_DIM(converted,0));
    return M;
  }
  
  template <class vector> PyObject * convert_vector_to_python(const vector & b,bool writeable,bool owns_data = false){
    init();
    long size = b.size();
    PyArrayObject * converted = (PyArrayObject *) PyArray_SimpleNewFromData(1,&size,py_scalar,(void*)b.data());
    if(!writeable) PyArray_CLEARFLAGS(converted,NPY_ARRAY_WRITEABLE);
    if(owns_data)  PyArray_ENABLEFLAGS(converted,NPY_ARRAY_OWNDATA);
    return (PyObject *)converted;
  }
  
  matrix_map convert_to_matrix(PyObject * Arr){
    init();
    auto flags = NPY_ARRAY_CARRAY;
    PyArrayObject * converted = (PyArrayObject *)PyArray_FromAny(Arr,PyArray_DescrFromType(py_scalar),2,2,flags,NULL);
    if(!converted) throw std::invalid_argument("cannot convert python object to matrix");
    matrix_map M((scalar *) PyArray_DATA(converted),PyArray_DIM(converted,0),PyArray_DIM(converted,1));
    return M;
  }
  
  template <class matrix> PyObject * convert_matrix_to_python(const matrix & b,bool writeable = false,bool owns_data = false){
    init();
    std::array<long,2> size = {{b.cols(),b.rows()}};
    PyArrayObject * converted = (PyArrayObject *) PyArray_SimpleNewFromData(2,size.data(),py_scalar,(void*)b.data());
    
    if( true ){
      PyArray_STRIDES(converted)[0] = (size_t)&b(0,1) - (size_t)&b(0,0);
      PyArray_STRIDES(converted)[1] = (size_t)&b(1,0) - (size_t)&b(0,0);
    
      PyArray_CLEARFLAGS(converted,NPY_ARRAY_C_CONTIGUOUS);
      PyArray_ENABLEFLAGS(converted, NPY_ARRAY_ALIGNED);
    }
    
    if(!writeable) PyArray_CLEARFLAGS(converted,NPY_ARRAY_WRITEABLE);
    if(owns_data)  PyArray_ENABLEFLAGS(converted,NPY_ARRAY_OWNDATA);
    return (PyObject *)converted;
  }
  
  struct c_complex{ double real; double imag; };
  inline scalar convert(c_complex z){ return scalar(z.real,z.imag); }

}

std::function<python::scalar(double)> get_1D_function_from_pointer(size_t f){
  using f_type = python::c_complex (*)(double);
  return [f](double x){ return python::convert(f_type(f)(x)); };
}

std::function<python::scalar(double,double)> get_2D_function_from_pointer(size_t f){
  using f_type = python::c_complex (*)(double,double);
  return [f](double x,double z){ return python::convert(f_type(f)(x,z)); };
}

std::function<python::scalar(double,double,double)> get_3D_function_from_pointer(size_t f){
  using f_type = python::c_complex (*)(double,double,double);
  return [f](double x,double y,double z){ return python::convert(f_type(f)(x,y,z)); };
}


class py_finite_difference_1D:public lars::finite_difference_1D{
  public:
  PyObject * get_field() const{ return python::convert_vector_to_python(finite_difference_1D::get_field(),false); }
  
  void set_F(size_t f){ F = get_2D_function_from_pointer(f); }
  void set_u0(size_t f){ u0 = get_2D_function_from_pointer(f); }
  void set_u0_boundary(size_t f){ u0_boundary = get_2D_function_from_pointer(f); }
};


class py_finite_difference_2D:public lars::finite_difference_2D{
public:
  PyObject * get_field(){ return python::convert_matrix_to_python(finite_difference_2D::get_field(),false); }
  PyObject * get_full_field(){ return python::convert_matrix_to_python(finite_difference_2D::get_full_field(),false); }

  void set_F(size_t f){ F = get_3D_function_from_pointer(f); }
  void set_u0(size_t f){ u0 = get_3D_function_from_pointer(f); }
  void set_u0_boundary(size_t f){ u0_boundary = get_3D_function_from_pointer(f); }
};

PyObject * create_2D_field_from_function(size_t f,unsigned nx,unsigned ny,double xmin,double xmax,double ymin,double ymax){
  python::matrix_map map((python::scalar*)malloc( sizeof(python::scalar) * nx * ny ),ny,nx);
  auto g = get_2D_function_from_pointer(f);
  
  auto dx = (xmax-xmin)/(nx-1);
  auto dy = (ymax-ymin)/(ny-1);
  
  parallel_for<unsigned>(0, nx, [&](unsigned xi){
    for(auto yi:range(ny)){
      auto x = xmin + xi * dx;
      auto y = ymin + yi * dy;
      map(yi,xi) = g(x,y);
    }
  });
  
  auto r = python::convert_matrix_to_python(map, true, true);
  return r;
}

PyObject * create_1D_field_from_function(size_t f,unsigned nx,double xmin,double xmax){
  python::vector_map map((python::scalar*)malloc( sizeof(python::scalar) * nx ),nx);
  auto g = get_1D_function_from_pointer(f);
  
  auto dx = (xmax-xmin)/(nx-1);
  
    for(auto xi:range(nx)){
      auto x = xmin + xi * dx;
      map(xi) = g(x);
    }
  
  auto r = python::convert_vector_to_python(map, true, true);
  return r;
}

using namespace boost::python;

BOOST_PYTHON_MODULE(_finitedifferences){
  
  def("create_2D_field_from_function", create_2D_field_from_function);
  def("create_1D_field_from_function", create_1D_field_from_function);
  
  class_<py_finite_difference_1D>("finite_difference_1D")
    .def("get_field",&py_finite_difference_1D::get_field)
    .def("set_F",&py_finite_difference_1D::set_F)
    .def("set_u0",&py_finite_difference_1D::set_u0)
    .def("set_u0_boundary",&py_finite_difference_1D::set_u0_boundary)
    .def_readwrite("dx",&py_finite_difference_1D::dx)
    .def_readwrite("dz",&py_finite_difference_1D::dz)
    .def_readwrite("z",&py_finite_difference_1D::z)
    .def_readwrite("A",&py_finite_difference_1D::A)
    .def("step",&py_finite_difference_1D::step)
    .def("resize",&py_finite_difference_1D::resize)
    .def("init",&py_finite_difference_1D::init)
  ;
  
  class_<py_finite_difference_2D>("finite_difference_2D")
  .def("get_field",&py_finite_difference_2D::get_field)
  .def("get_full_field",&py_finite_difference_2D::get_full_field)
  .def("set_F",&py_finite_difference_2D::set_F)
  .def("set_u0",&py_finite_difference_2D::set_u0)
  .def("set_u0_boundary",&py_finite_difference_2D::set_u0_boundary)
  .def_readwrite("dx",&py_finite_difference_2D::dx)
  .def_readwrite("dy",&py_finite_difference_2D::dy)
  .def_readwrite("dz",&py_finite_difference_2D::dz)
  .def_readwrite("z",&py_finite_difference_2D::z)
  .def_readwrite("A",&py_finite_difference_2D::A)
  .def("step",&py_finite_difference_2D::step)
  .def("resize",&py_finite_difference_2D::resize)
  .def("init",&py_finite_difference_2D::init)
  ;

}


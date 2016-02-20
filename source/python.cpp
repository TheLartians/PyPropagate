


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
    matrix_map M((scalar *) PyArray_DATA(converted),PyArray_DIM(converted,1),PyArray_DIM(converted,0));
    return M;
  }
  
  PyObject * convert_matrix_to_python(const finite_difference_2D::field & b,bool writeable = false,bool owns_data = false){
    init();
    
    std::array<long,2> size = {{b.cols(),b.rows()}};
    PyArrayObject * converted = (PyArrayObject *) PyArray_SimpleNewFromData(2,size.data(),py_scalar,(void*)b.data());
    
    if( false ){ // For non contignous types
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
  void set_field(PyObject *arr){ lars::finite_difference_1D::set_field(python::convert_to_vector(arr)); }
  
  void set_F(size_t f){ F = get_2D_function_from_pointer(f); }
  void set_u_boundary(size_t f){ u_boundary = get_2D_function_from_pointer(f); }
};

class py_finite_difference_2D:public lars::finite_difference_2D{
public:
  PyObject * get_field(){ return python::convert_matrix_to_python(finite_difference_2D::get_field(),false); }
  void set_field(PyObject *arr){ lars::finite_difference_2D::set_field( python::convert_to_matrix(arr) ); }

  void set_F(size_t f){ F = get_3D_function_from_pointer(f); }
  void set_u_boundary(size_t f){ u_boundary = get_3D_function_from_pointer(f); }
};

void tridiagonal(PyObject * v1,PyObject * v2,PyObject * v3,PyObject * v4,PyObject * v5,PyObject * v6){
  auto a1 = python::convert_to_vector(v1),a2 = python::convert_to_vector(v2),a3 = python::convert_to_vector(v3),a4 = python::convert_to_vector(v4),a5 = python::convert_to_vector(v5),a6 = python::convert_to_vector(v6);
  lars::algebra::tridiagonal(a1, a2, a3, a4, a5, a6);
}

using namespace boost::python;

BOOST_PYTHON_MODULE(_pypropagate){
  
  class_<py_finite_difference_1D>("finite_difference_1D")
    .def("get_field",&py_finite_difference_1D::get_field)
    .def("set_field",&py_finite_difference_1D::set_field)
    .def("set_F",&py_finite_difference_1D::set_F)
    .def("set_u_boundary",&py_finite_difference_1D::set_u_boundary)
    .def_readwrite("xmin",&py_finite_difference_1D::xmin)
    .def_readwrite("xmax",&py_finite_difference_1D::xmax)
    .def_readwrite("dz",&py_finite_difference_1D::dz)
    .def_readwrite("z",&py_finite_difference_1D::z)
    .def_readwrite("A",&py_finite_difference_1D::A)
    .def_readwrite("constant_F",&py_finite_difference_1D::constant_F)
    .def("step",&py_finite_difference_1D::step)
    .def("init",&py_finite_difference_1D::init)
  ;
  
  class_<py_finite_difference_2D>("finite_difference_2D")
  .def("get_field",&py_finite_difference_2D::get_field)
  .def("set_field",&py_finite_difference_2D::set_field)
  .def("set_F",&py_finite_difference_2D::set_F)
  .def("set_u_boundary",&py_finite_difference_2D::set_u_boundary)
  .def_readwrite("xmin",&py_finite_difference_2D::xmin)
  .def_readwrite("xmax",&py_finite_difference_2D::xmax)
  .def_readwrite("ymin",&py_finite_difference_2D::ymin)
  .def_readwrite("ymax",&py_finite_difference_2D::ymax)
  .def_readwrite("dz",&py_finite_difference_2D::dz)
  .def_readwrite("z",&py_finite_difference_2D::z)
  .def_readwrite("A",&py_finite_difference_2D::A)
  .def_readwrite("constant_F",&py_finite_difference_2D::constant_F)
  .def("step",&py_finite_difference_2D::step)
  .def("init",&py_finite_difference_2D::init)
  ;

}


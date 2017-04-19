
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <boost/python.hpp>
#include <numpy/arrayobject.h>

#include <stdexcept>
#include <sstream>

#include <unordered_map>
#include <functional>
#include <typeinfo>

#include "finite_difference.h"
#include "ring_derivative.h"


namespace python_converters{
  
  const auto py_complex = NPY_COMPLEX128;
 
  using TypeInfoRef = std::reference_wrapper<const std::type_info>;
 
  struct TypeInfoHasher {
    std::size_t operator()(TypeInfoRef code) const{
        return code.get().hash_code();
    }
  };
 
  struct TypeInfoEqualTo {
    bool operator()(TypeInfoRef lhs, TypeInfoRef rhs) const {
        return lhs.get() == rhs.get();
    }
  };
  
  const std::unordered_map<TypeInfoRef, decltype(NPY_BOOL), TypeInfoHasher, TypeInfoEqualTo> numpy_types = {{typeid(bool),NPY_BOOL},
                                                                                                             {typeid(int),NPY_INT},
                                                                                                             {typeid(double),NPY_DOUBLE}};

  struct init_numpy{ init_numpy(){ import_array(); } } initializer;
  
  using namespace boost::python;
  using namespace lars;
  
  template <class T,size_t N> using MappedNDArray = NDArrayBase<T, DynamicIndexTuple<N>, DynamicIndexTuple<N>, DynamicIndex, BorrowedData<T>, BasicNDArrayCreator<HeapNDArray>>;

  template <class T,size_t N> MappedNDArray<T,N> numpy_to_ndarray(object data_obj){
    PyArrayObject * arr = reinterpret_cast<PyArrayObject*>(data_obj.ptr());

    if(PyArray_NDIM(arr) != N) throw std::runtime_error("array is not two-dimensional");
    
    auto it = numpy_types.find(typeid(T));
    if(it == numpy_types.end() || PyArray_TYPE(arr) != it->second) throw std::runtime_error("array has wrong type");   
 
    auto npyshape = PyArray_SHAPE(arr);
    auto npystrides = PyArray_STRIDES(arr);
    auto data = PyArray_BYTES(arr);
    
    using Index = DynamicIndexTuple<N>;
    Index strides,shape;
    
    shape.apply([&](size_t idx,size_t &val){ val = npyshape[idx]; }); 
    strides.apply([&](size_t idx,size_t &val){ val = npystrides[idx] / sizeof(T); }); 
    
    return MappedNDArray<double,2>(shape,strides,0,(T*)data); 
  }
 
  PyObject * array_1D_as_numpy(finite_differences::array_1D & array){
    long size = array.size();
    PyArrayObject * converted = (PyArrayObject *) PyArray_SimpleNewFromData(1,&size,py_complex,(void*)array.data());
    return (PyObject *)converted;
  }
  
  PyObject * array_2D_as_numpy(finite_differences::array_2D & array){
    long size[2] = {static_cast<long>(array.size()),static_cast<long>(array[0].size())};
    PyArrayObject * converted = (PyArrayObject *) PyArray_SimpleNewFromData(2,size,py_complex,(void*)array.data());
    return (PyObject *)PyArray_Transpose(converted,nullptr);
  }
  
  template <typename T> std::string to_string(const T &o){
    std::stringstream stream; stream << o; return stream.str();
  }
  
}

BOOST_PYTHON_MODULE(_pypropagate){
  using namespace boost::python;
  using namespace lars;
  
  def("ring_derivative_2D",+[](object pyarr,object pydx,object pydy,double s){ 
    auto arr = python_converters::numpy_to_ndarray<double,2>(pyarr); 
    auto dx = python_converters::numpy_to_ndarray<double,2>(pydx); 
    auto dy = python_converters::numpy_to_ndarray<double,2>(pydy);
    ring_derivative<0>(arr,dx,s); 
    ring_derivative<1>(arr,dy,s); 
  },args("array","dx","dy","s"));

  class_<finite_differences::array_1D>("Array1D")
  .def("as_numpy",python_converters::array_1D_as_numpy)
  .def("resize",+[](finite_differences::array_1D &arr,size_t size){ arr.resize(size); })
  .def("__str__",python_converters::to_string<finite_differences::array_1D>)
  ;
  
  class_<finite_differences::array_2D>("Array2D")
  .def("as_numpy",python_converters::array_2D_as_numpy)
  .def("resize",+[](finite_differences::array_2D &arr,size_t a,size_t b){ arr.resize(a,b); })
  .def("__str__",python_converters::to_string<finite_differences::array_2D>)
  ;
  
  class_<finite_difference_AF>("finite_difference_AF")
  .def_readwrite("ra",&finite_difference_AF::ra)
  .def_readwrite("rf",&finite_difference_AF::rf)
  .def_readwrite("u",&finite_difference_AF::u)
  .def("step", &finite_difference_AF::step)
  .def("update", &finite_difference_AF::update)
  .def("resize", &finite_difference_AF::resize)
  ;
  
  class_<finite_difference_ACF>("finite_difference_ACF")
  .def_readwrite("ra",&finite_difference_ACF::ra)
  .def_readwrite("rc",&finite_difference_ACF::rc)
  .def_readwrite("rf",&finite_difference_ACF::rf)
  .def_readwrite("u",&finite_difference_ACF::u)
  .def_readwrite("thread_count",&finite_difference_ACF::thread_count)
  .def("step_1", &finite_difference_ACF::step_1)
  .def("step_2", &finite_difference_ACF::step_2)
  .def("update", &finite_difference_ACF::update)
  .def("resize", &finite_difference_ACF::resize)
  ;
  
  class_<finite_difference_A0F>("finite_difference_A0F")
  .def_readwrite("ra",&finite_difference_A0F::ra)
  .def_readwrite("rf",&finite_difference_A0F::rf)
  .def_readwrite("u",&finite_difference_A0F::u)
  .def_readwrite("thread_count",&finite_difference_A0F::thread_count)
  .def("step", &finite_difference_A0F::step)
  .def("update", &finite_difference_A0F::update)
  .def("resize", &finite_difference_A0F::resize)
  ;

  class_<finite_difference_ABC>("finite_difference_ABC")
  .def_readwrite("ra",&finite_difference_ABC::ra)
  .def_readwrite("rb",&finite_difference_ABC::rb)
  .def_readwrite("rc",&finite_difference_ABC::rc)
  .def_readwrite("rz",&finite_difference_ABC::rz)
  .def_readwrite("thread_count",&finite_difference_ABC::thread_count)
  .def_readwrite("u",&finite_difference_ABC::u)
  .def("step", &finite_difference_ABC::step)
  .def("update", &finite_difference_ABC::update)
  .def("resize", &finite_difference_ABC::resize)
  ;
}

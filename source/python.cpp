
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <boost/python.hpp>
#include <numpy/arrayobject.h>

#include <stdexcept>
#include <sstream>

#include "finite_difference.h"

namespace python_converters{
  
  const auto py_scalar = NPY_COMPLEX128;
  
  struct init_numpy{ init_numpy(){ import_array(); } } initializer;
  
  using namespace boost::python;
  using namespace lars;

  PyObject * array_1D_as_numpy(finite_differences::array_1D & array){
    long size = array.size();
    PyArrayObject * converted = (PyArrayObject *) PyArray_SimpleNewFromData(1,&size,py_scalar,(void*)array.data());
    return (PyObject *)converted;
  }
  
  PyObject * array_2D_as_numpy(finite_differences::array_2D & array){
    long size[2] = {static_cast<long>(array.size()),static_cast<long>(array[0].size())};
    PyArrayObject * converted = (PyArrayObject *) PyArray_SimpleNewFromData(2,size,py_scalar,(void*)array.data());
    return (PyObject *)converted;
  }
  
  template <typename T> std::string to_string(const T &o){
    std::stringstream stream; stream << o; return stream.str();
  }
  
}

class ScopedGILRelease {
public:
  inline ScopedGILRelease() { m_thread_state = PyEval_SaveThread(); }
  inline ~ScopedGILRelease() {
  PyEval_RestoreThread(m_thread_state);m_thread_state = NULL;}
private:
  PyThreadState* m_thread_state;
};

BOOST_PYTHON_MODULE(_pypropagate){
  using namespace boost::python;
  using namespace lars;
  PyEval_InitThreads();
  
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
  .def("step", +[](finite_difference_AF &fd){ ScopedGILRelease gil; fd.step(); })
  .def("update", &finite_difference_AF::update)
  .def("resize", &finite_difference_AF::resize)
  ;
  
  class_<finite_difference_ACF>("finite_difference_ACF")
  .def_readwrite("ra",&finite_difference_ACF::ra)
  .def_readwrite("rc",&finite_difference_ACF::rc)
  .def_readwrite("rf",&finite_difference_ACF::rf)
  .def_readwrite("u",&finite_difference_ACF::u)
  .def("step_1", +[](finite_difference_ACF &fd){ ScopedGILRelease gil; fd.step_1(); })
  .def("step_2", +[](finite_difference_ACF &fd){ ScopedGILRelease gil; fd.step_2(); })
  .def("update", &finite_difference_ACF::update)
  .def("resize", &finite_difference_ACF::resize)
  ;
  
  class_<finite_difference_A0F>("finite_difference_A0F")
  .def_readwrite("ra",&finite_difference_A0F::ra)
  .def_readwrite("rf",&finite_difference_A0F::rf)
  .def_readwrite("u",&finite_difference_A0F::u)
  .def("step", +[](finite_difference_A0F &fd){ ScopedGILRelease gil; fd.step(); })
  .def("update", &finite_difference_A0F::update)
  .def("resize", &finite_difference_A0F::resize)
  ;

  class_<finite_difference_ABC>("finite_difference_ABC")
  .def_readwrite("ra",&finite_difference_ABC::ra)
  .def_readwrite("rb",&finite_difference_ABC::rb)
  .def_readwrite("rc",&finite_difference_ABC::rc)
  .def_readwrite("rz",&finite_difference_ABC::rz)
  .def_readwrite("u",&finite_difference_ABC::u)
  .def("step", +[](finite_difference_ABC &fd){ ScopedGILRelease gil; fd.step(); })
  .def("update", &finite_difference_ABC::update)
  .def("resize", &finite_difference_ABC::resize)
  ;
}





from .propagator import Propagator
from .._pypropagate import finite_difference_1D,finite_difference_2D
from ..coordinate_ndarray import CoordinateNDArray
from ..compile_sympy import compile_sympy_expression,compile_sympy_expressions,ctypes
import numpy

class FiniteDifferencesPropagator1D(Propagator):

    ndim = 1
    dtype = numpy.complex128
    
    def __init__(self,settings,initial=None):
        super(FiniteDifferencesPropagator1D,self).__init__(settings)
    
        s = settings.simulation_box
        fd = settings.finitedifferences
        
        y = (s.ymin + s.ymax)/2
        
        F,u0,u0_boundary = settings.get_unitless((fd.F.subs(s.y,y),fd.u0.subs(s.y,y),fd.u0_boundary.subs(s.y,y)))
        A = settings.get(fd.A,complex)
        dx,dz,self.__zmin = settings.get((s.dx,s.dz,s.zmin),float)
        nx = settings.get(s.nx,int)
    
        lib = compile_sympy_expressions([("F",F,(s.x,s.z),ctypes.c_complex),
                                   ("u0_boundary",u0_boundary,(s.x,s.z),ctypes.c_complex),
                                   ("u0",u0,(s.x,s.z),ctypes.c_complex)])
        
        self._solver = finite_difference_1D()
        
        self._solver.set_F(lib.F.address)
        self._solver.set_u0(lib.u0.address)
        self._solver.set_u0_boundary(lib.u0_boundary.address)
        self._solver.A = complex(A)
        self._solver.dz = dz
        self._solver.dx = dx
        self._solver.z = self.__zmin
        self._solver.resize(nx)
        self._solver.init()
        
        if initial != None:
            self.set_field(initial)
        
    def reset(self):
        self._solver.z = self.__zmin
        self._solver.init()
        
    def _step(self):
        self._solver.step()
    
    def _get_field(self):
        return self._solver.get_field()
    
    def _set_field(self,field):
        self._solver.set_field(field)

    
class FiniteDifferencesPropagator2D(Propagator):

    ndim = 2
    dtype = numpy.complex128
    
    def __init__(self,settings,initial = None):
        super(FiniteDifferencesPropagator2D,self).__init__(settings)
    
        s = settings.simulation_box
        fd = settings.finitedifferences
        
        F,u0,u0_boundary = settings.get_unitless((fd.F,fd.u0,fd.u0_boundary))
        A = settings.get(fd.A,complex)
        dx,dy,dz,self.__zmin = settings.get((s.dx,s.dy,s.dz,s.zmin),float)
        nx,ny = settings.get((s.nx,s.ny),int)
    
        lib = compile_sympy_expressions([("F",F,(s.x,s.y,s.z),ctypes.c_complex),
                                   ("u0_boundary",u0_boundary,(s.x,s.y,s.z),ctypes.c_complex),
                                   ("u0",u0,(s.x,s.y,s.z),ctypes.c_complex)])
        
        self._solver = finite_difference_2D()
        
        self._solver.set_F(lib.F.address)
        self._solver.set_u0(lib.u0.address)
        self._solver.set_u0_boundary(lib.u0_boundary.address)
        self._solver.A = complex(A)
        self._solver.dz = dz
        self._solver.dy = dy
        self._solver.dx = dx
        self._solver.z = self.__zmin
        self._solver.resize(nx,ny)
        self._solver.init()
        
        if initial != None:
            self.set_field(initial)
        
    def _reset(self):
        self._solver.z = self.__zmin
        self._solver.init()
        
    def _step(self):
        self._solver.step()
    
    def _get_field(self):
        return self._solver.get_field().transpose()
    
    def _set_field(self,field):
        self._solver.set_field(field)



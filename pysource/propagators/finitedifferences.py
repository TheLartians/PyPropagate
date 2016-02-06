
from .propagator import Propagator
from .._pypropagate import finite_difference_1D,finite_difference_2D
from ..coordinate_ndarray import CoordinateNDArray
import pycas as pc
import numpy as np

class FiniteDifferencesPropagator1D(Propagator):

    ndim = 1
    dtype = np.complex128
    
    def __init__(self,settings):
        super(FiniteDifferencesPropagator1D,self).__init__(settings)
    
        sb = settings.simulation_box
        pe = settings.paraxial_equation
        
        F,A,u0,u_boundary = settings.get_unitless(pc.Tuple(pe.F,pe.A, pe.u0, pe.u_boundary).subs(sb.y,sb.fy))
        xmin,xmax,zmin,dz = settings.get_as((sb.xmin,sb.xmax,sb.zmin,sb.dz),float)

        lib = pc.ccompile(
                pc.FunctionDefinition("F",(sb.x,sb.z),F,return_type=pc.Types.Complex),
                pc.FunctionDefinition("u_boundary",(sb.x,sb.z),u_boundary,return_type=pc.Types.Complex),
        )

        self._solver = finite_difference_1D()

        self._solver.set_F(lib.F.address())
        self._solver.set_u_boundary(lib.u_boundary.address())
        self._solver.A = complex(A)
        self._solver.xmin = xmin
        self._solver.xmax = xmax
        self._solver.z = zmin
        self._solver.dz = dz
        self._solver.constant_F = self._F_is_constant

        self._set_initial_field(settings)

    def _set_z(self,z):
        self._solver.z = z

    def _step(self):
        self._solver.step()
    
    def _get_field(self):
        return self._solver.get_field()
    
    def _set_field(self,field):
        self._solver.set_field(field)
        self._solver.init()

    
class FiniteDifferencesPropagator2D(Propagator):

    ndim = 2
    dtype = np.complex128
    
    def __init__(self,settings):
        super(FiniteDifferencesPropagator2D,self).__init__(settings)
    
        sb = settings.simulation_box
        pe = settings.paraxial_equation

        F,A,u0,u_boundary = settings.get_unitless(( pe.F, pe.A, pe.u0, pe.u_boundary, ))
        xmin,xmax,ymin,ymax,zmin,dz = settings.get_as((sb.xmin,sb.xmax,sb.ymin,sb.ymax,sb.zmin,sb.dz),float)

        lib = pc.ccompile(
                pc.FunctionDefinition("F",(sb.x,sb.y,sb.z),F,return_type=pc.Types.Complex),
                pc.FunctionDefinition("u_boundary",(sb.x,sb.y,sb.z),u_boundary,return_type=pc.Types.Complex),
        )

        self._solver = finite_difference_2D()
        self._solver.set_F(lib.F.address())
        self._solver.set_u_boundary(lib.u_boundary.address())
        self._solver.A = complex(A)
        self._solver.xmin = xmin
        self._solver.xmax = xmax
        self._solver.ymin = ymin
        self._solver.ymax = ymax
        self._solver.dz = dz
        self._solver.z = zmin
        self._solver.constant_F = self._F_is_constant

        self._set_initial_field(settings)

    def _set_z(self,z):
        self._solver.z = z

    def _step(self):
        self._solver.step()
    
    def _get_field(self):
        return self._solver.get_field().transpose()
    
    def _set_field(self,field):
        self._solver.set_field(field)
        self._solver.init()



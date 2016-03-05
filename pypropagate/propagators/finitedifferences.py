
from .propagator import Propagator
import expresso.pycas as pc
import numpy as np

class FiniteDifferencesPropagator1D(Propagator):

    ndim = 1
    dtype = np.complex128
    
    def __init__(self,settings):        
        super(FiniteDifferencesPropagator1D,self).__init__(settings)
        from _pypropagate import finite_difference_aF

        pde = settings.partial_differential_equation        
        ra = settings.get_as(pde.ra,complex)
        self.__u_boundary,self.__rf = self._get_evaluators([ pde.u_boundary, pde.rf ],settings,return_type=pc.Types.Complex,compile_to_c = True,parallel=False)

        self._solver = finite_difference_aF()
        self._solver.resize(self._nx)
        self._solver.ra = ra       

        self._set_initial_field(settings)
        self.__boundary_values = np.array([self._nxmin,self._nxmax],dtype=float)
        
        self._update()
        
    def _update(self):
        self._solver.u.as_numpy()[[0,-1]] = self.__u_boundary(self.__boundary_values,[self._current_nz]*2)
        self._solver.update()
        self.__rf(*self._get_coordinates(),res=self._solver.rf.as_numpy())
        
    def _step(self):
        self._update()
        self._solver.step()
    
    def _get_field(self):
        return self._solver.u.as_numpy()
    
    def _set_field(self,field):
        self._solver.u.as_numpy()[:] = field


    
class FiniteDifferencesPropagator2D(Propagator):

    ndim = 2
    dtype = np.complex128
    
    def __init__(self,settings):
        from _pypropagate import finite_difference_2D

        super(FiniteDifferencesPropagator2D,self).__init__(settings)
    
        sb = settings.simulation_box
        pe = settings.paraxial_equation

        F,A,u_boundary = settings.get_unitless(( pe.F, pe.A, pe.u_boundary, ))
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
        self._solver.set_field(field.transpose())
        self._solver.init()


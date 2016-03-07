
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
        self.__u_boundary,self.__rf = self._get_evaluators([ pde.u_boundary, pde.rf ], settings, return_type=pc.Types.Complex, compile_to_c = not self._F_is_constant_in_z, parallel=False)

        self._solver = finite_difference_aF()
        self._solver.resize(self._nx)
        self._solver.ra = ra       

        self._set_initial_field(settings)
        self.__boundary_values = np.array([0,self._nx-1],dtype=np.uint)
        self.__z_values = np.array([0,0],dtype=np.uint)

        self._reset()

    def _reset(self):
        self.__rf(*self._get_coordinates(),res=self._solver.rf.as_numpy())
        self._solver.update()
        super(FiniteDifferencesPropagator1D,self)._reset()
        self.__rf(*self._get_coordinates(),res=self._solver.rf.as_numpy())

    def _update(self):
        self._solver.update()
        self.__z_values.fill(self._i)
        self._solver.u.as_numpy()[self.__boundary_values] = self.__u_boundary(self.__boundary_values,self.__z_values)
        if not self._F_is_constant_in_z:
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
        super(FiniteDifferencesPropagator2D,self).__init__(settings)
        from _pypropagate import finite_difference_acF

        pde = settings.partial_differential_equation        
        ra = settings.get_as(pde.ra/2,complex)
        rc = settings.get_as(pde.rc/2,complex)

        self.__rf, = self._get_evaluators([ pde.rf/2 ],settings,return_type=pc.Types.Complex,compile_to_c = True ,parallel=True)
        self.__u_boundary, =  self._get_evaluators([ pde.u_boundary ],settings,return_type=pc.Types.Complex,compile_to_c = True ,parallel=False)

        self._solver = finite_difference_acF()
        self._solver.resize(self._nx,self._ny)
        self._solver.ra = ra
        self._solver.rc = rc
        
        d,u,l,r = [(self._get_x_coordinates(),np.zeros(self._nx,dtype = np.uint)),
                   (self._get_x_coordinates(),np.ones(self._nx,dtype = np.uint)*self._ny),
                   (np.zeros(self._ny,dtype = np.uint),
                    self._get_y_coordinates()),
                   (np.ones(self._ny,dtype = np.uint)*self._nx,self._get_y_coordinates())]

        self.__boundary_values = [np.concatenate([v[0] for v in d,u,l,r]),
                                  np.concatenate([v[1] for v in d,u,l,r]),
                                  np.zeros(2*self._nx+2*self._ny,dtype=np.uint)]

        self._set_initial_field(settings)
        self._reset()
        
    def _reset(self):
        super(FiniteDifferencesPropagator2D,self)._reset()
        self._update(True)
        super(FiniteDifferencesPropagator2D,self)._reset()
        self._update(True)

    def _update_boundary(self):
        self.__boundary_values[2].fill(self._i)
        boundary = self.__u_boundary(*self.__boundary_values)

        self._solver.u.as_numpy()[:,0] = boundary[0:self._nx]
        self._solver.u.as_numpy()[:,-1] = boundary[self._nx:2*self._nx]
        self._solver.u.as_numpy()[0,:] = boundary[2*self._nx:2*self._nx+self._ny]
        self._solver.u.as_numpy()[-1,:] = boundary[2*self._nx+self._ny:]
        
    def _update(self,force_update_F = False):
        self._solver.update()
        self._update_boundary()
        if (not self._F_is_constant_in_z) or force_update_F:
            self.__rf(*self._get_coordinates(),res=self._solver.rf.as_numpy())
        
    def _step(self):
        self._update()
        self._solver.step_1()
        self._update()
        self._solver.step_2()

    def _get_field(self):
        return self._solver.u.as_numpy()
    
    def _set_field(self,field):
        self._solver.u.as_numpy()[:] = field



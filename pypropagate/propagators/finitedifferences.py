
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
        sb = settings.simulation_box

        self._2step = settings.get_numeric( pc.equal(pde.C, 0)  ) != pc.S(True)
        sf = 0.5 if self._2step else 1

        ra = settings.get_as(pde.ra*sf,complex)
        rc = settings.get_as(pde.rc*sf,complex)

        evaluators = self._get_evaluators([ (pde.rf*sf),
                                            (pde.rf*sf).subs(sb.z,sb.z-sb.dz*sf),
                                            pde.u_boundary,
                                            pde.u_boundary.subs(sb.z,sb.z-sb.dz*sf) ],
                                          settings,return_type=pc.Types.Complex,compile_to_c = True,parallel=True)

        self.__rf = evaluators[:2]
        self.__u_boundary = evaluators[2:]

        self._solver = finite_difference_acF()
        self._solver.resize(self._nx,self._ny)

        self._solver.ra = ra
        self._solver.rc = rc
        
        d,u,l,r = [(self._get_x_coordinates(),np.zeros(self._nx,dtype = np.uint)),
                   (self._get_x_coordinates(),np.ones(self._nx,dtype = np.uint)*(self._ny-1)),
                   (np.zeros(self._ny,dtype = np.uint),
                    self._get_y_coordinates()),
                   (np.ones(self._ny,dtype = np.uint)*(self._nx-1),self._get_y_coordinates())]

        self.__boundary_values = [np.concatenate([v[0] for v in d,u,l,r]),
                                  np.concatenate([v[1] for v in d,u,l,r]),
                                  np.zeros(2*self._nx+2*self._ny,dtype=np.uint)]

        self._set_initial_field(settings)
        self._reset()
        
    def _reset(self):
        self.__rf[1](*self._get_coordinates(),res=self._solver.rf.as_numpy())
        self._solver.u.as_numpy().fill(0)
        self._solver.update()
        self.__rf[0](*self._get_coordinates(),res=self._solver.rf.as_numpy())
        super(FiniteDifferencesPropagator2D,self)._reset()

    def _update_boundary(self,half_step):
        self.__boundary_values[2].fill(self._i)
        boundary = self.__u_boundary[half_step](*self.__boundary_values)
        u  = self._solver.u.as_numpy()

        u[:,0] = boundary[0:self._nx]
        u[:,-1] = boundary[self._nx:2*self._nx]
        u[0,:] = boundary[2*self._nx:2*self._nx+self._ny]
        u[-1,:] = boundary[2*self._nx+self._ny:]
        
    def _update(self,half_step):
        self._solver.update()
        self._update_boundary(half_step)
        if (not self._F_is_constant_in_z):
            self.__rf[half_step](*self._get_coordinates(),res=self._solver.rf.as_numpy())
        
    def _step(self):
        if self._2step:
            self._update(True)
            self._solver.step_1()
        self._update(False)
        self._solver.step_2()

    def _get_field(self):
        return self._solver.u.as_numpy()
    
    def _set_field(self,field):
        self._solver.u.as_numpy()[:] = field



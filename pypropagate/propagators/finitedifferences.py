
# TODO: update A in 1D FD propagation
# TODO: transpose arrays in 1D FD propagation

from .propagator import Propagator
import expresso.pycas as pc
import numpy as np

class FiniteDifferencesPropagator1D(Propagator):

    ndim = 1
    dtype = np.complex128
    
    def __init__(self,settings):        
        super(FiniteDifferencesPropagator1D,self).__init__(settings)
        from _pypropagate import finite_difference_AF

        pde = settings.partial_differential_equation        

        self.__u_boundary,self.__rf,self.__ra = self._get_evaluators([ pde.u_boundary, pde.rf, pde.ra ], settings, return_type=pc.Types.Complex, compile_to_c = not self._F_is_constant_in_z, parallel=False)

        self._solver = finite_difference_AF()
        self._solver.resize(self._nx)

        self._set_initial_field(settings)
        self.__boundary_values = np.array([0,self._nx-1],dtype=np.uint)
        self.__z_values = np.array([0,0],dtype=np.uint)

        self._reset()

    def _reset(self):
        self.__ra(*self._get_indices(), res=self._solver.ra.as_numpy())
        self.__rf(*self._get_indices(), res=self._solver.rf.as_numpy())
        self._solver.update()
        super(FiniteDifferencesPropagator1D,self)._reset()
        self.__ra(*self._get_indices(), res=self._solver.ra.as_numpy())
        self.__rf(*self._get_indices(), res=self._solver.rf.as_numpy())

    def _update(self):
        self._solver.update()
        self.__z_values.fill(self._i)
        self._solver.u.as_numpy()[self.__boundary_values] = self.__u_boundary(self.__boundary_values,self.__z_values)
        if not self._F_is_constant_in_z:
            self.__rf(*self._get_indices(), res=self._solver.rf.as_numpy())

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
        from _pypropagate import finite_difference_ACF,finite_difference_A0F

        pde = settings.partial_differential_equation
        sb = settings.simulation_box

        self.__C_is_zero = settings.get_numeric(pde.C) == pc.S(0)

        sf = 0.5 if not self.__C_is_zero else 1

        z,dz = sb.coordinates[2].symbol,sb.coordinates[2].step

        evaluators = self._get_evaluators([ (pde.ra*sf),
                                            (pde.ra*sf).subs(z,z-dz*sf),
                                            (pde.rc*sf),
                                            (pde.rc*sf).subs(z,z-dz*sf),
                                            (pde.rf*sf),
                                            (pde.rf*sf).subs(z,z-dz*sf),
                                            pde.u_boundary,
                                            pde.u_boundary.subs(z,z-dz*sf) ],
                                          settings,return_type=pc.Types.Complex,compile_to_c = True,parallel=True)

        self.__ra = evaluators[0:2]
        self.__rc = evaluators[2:4]
        self.__rf = evaluators[4:6]
        self.__u_boundary = evaluators[6:8]

        self._solver = finite_difference_ACF() if not self.__C_is_zero else finite_difference_A0F()
        self._solver.resize(self._nx,self._ny)

        d,u,l,r = [(self._get_x_indices(), np.zeros(self._nx, dtype = np.uint)),
                   (self._get_x_indices(), np.ones(self._nx, dtype = np.uint) * (self._ny - 1)),
                   (np.zeros(self._ny,dtype = np.uint),
                    self._get_y_indices()),
                   (np.ones(self._ny,dtype = np.uint)*(self._nx-1),self._get_y_indices())]

        self.__boundary_values = [np.concatenate([v[0] for v in d,u,l,r]),
                                  np.concatenate([v[1] for v in d,u,l,r]),
                                  np.zeros(2*self._nx+2*self._ny,dtype=np.uint)]

        self._set_initial_field(settings)
        self._reset()
        
    def _reset(self):
        self.__ra[1](*self._get_transposed_indices(), res=self._solver.ra.as_numpy().transpose())
        if not self.__C_is_zero: self.__rc[1](*self._get_transposed_indices(), res=self._solver.rc.as_numpy().transpose())
        self.__rf[1](*self._get_transposed_indices(), res=self._solver.rf.as_numpy().transpose())
        self._solver.u.as_numpy().fill(0)
        self._solver.update()
        self.__ra[0](*self._get_transposed_indices(), res=self._solver.ra.as_numpy().transpose())
        if not self.__C_is_zero: self.__rc[0](*self._get_transposed_indices(), res=self._solver.rc.as_numpy().transpose())
        self.__rf[0](*self._get_transposed_indices(), res=self._solver.rf.as_numpy().transpose())
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
            self.__ra[half_step](*self._get_transposed_indices(), res=self._solver.ra.as_numpy().transpose())
            if not self.__C_is_zero: self.__rc[half_step](*self._get_transposed_indices(), res=self._solver.rc.as_numpy().transpose())
            self.__rf[half_step](*self._get_transposed_indices(), res=self._solver.rf.as_numpy().transpose())

    def _step(self):
        if not self.__C_is_zero:
            self._update(True)
            self._solver.step_1()
            self._update(False)
            self._solver.step_2()
        else:
            self._update(False)
            self._solver.step()

    def _get_field(self):
        return self._solver.u.as_numpy()
    
    def _set_field(self,field):
        self._solver.u.as_numpy()[:] = field


class FiniteDifferencesPropagatorRS(Propagator):
    ndim = 1
    dtype = np.complex128

    def __init__(self, settings):
        super(FiniteDifferencesPropagatorRS, self).__init__(settings)
        from _pypropagate import finite_difference_ABC

        pde = settings.partial_differential_equation
        x, y, z = settings.simulation_box.coordinates

        from ..units import m
        evaluators = self._get_evaluators([ pde.u_boundary,
                                           x.symbol * pde.A / (2 * x.step ** 2),
                                           pde.A / (4 * x.step),
                                           x.symbol * pde.F / 2,
                                           x.symbol / z.step
                                           ], settings, return_type=pc.Types.Complex,
                                          compile_to_c=not self._F_is_constant_in_z, parallel=False)

        self.__u_boundary, self.__ra, self.__rb, self.__rc, self.__rz = evaluators

        self._solver = finite_difference_ABC()
        self._solver.resize(self._nx, 1)

        self._set_initial_field(settings)
        self.__boundary_values = np.array([0, self._nx - 1], dtype=np.uint)
        self.__z_values = np.array([0, 0], dtype=np.uint)

        self._reset()

    def _reset(self):
        self.__ra(*self._get_transposed_indices(), res=self._solver.ra.as_numpy().transpose()[:])
        self.__rb(*self._get_transposed_indices(), res=self._solver.rb.as_numpy().transpose()[:])
        self.__rc(*self._get_transposed_indices(), res=self._solver.rc.as_numpy().transpose()[:])
        self.__rz(*self._get_transposed_indices(), res=self._solver.rz.as_numpy().transpose()[:])
        self._solver.update()
        super(FiniteDifferencesPropagatorRS, self)._reset()
        self.__ra(*self._get_transposed_indices(), res=self._solver.ra.as_numpy().transpose()[:])
        self.__rb(*self._get_transposed_indices(), res=self._solver.rb.as_numpy().transpose()[:])
        self.__rc(*self._get_transposed_indices(), res=self._solver.rc.as_numpy().transpose()[:])
        self.__rz(*self._get_transposed_indices(), res=self._solver.rz.as_numpy().transpose()[:])

    def _update(self):
        self._solver.update()
        self.__z_values.fill(self._i)

        self._get_field()[self.get_boundary_indices()] = self.__u_boundary(self.get_boundary_indices(),self.__z_values)
        if not self._F_is_constant_in_z:
            self.__ra(*self._get_transposed_indices(), res=self._solver.ra.as_numpy().transpose()[:])
            self.__rb(*self._get_transposed_indices(), res=self._solver.rb.as_numpy().transpose()[:])
            self.__rc(*self._get_transposed_indices(), res=self._solver.rc.as_numpy().transpose()[:])

    def _step(self):
        self._update()
        self._solver.step()

    def _get_field(self):
        return self._solver.u.as_numpy().transpose()[0]

    def _set_field(self, field):
        self._solver.u.as_numpy().transpose()[0][:] = field


from .propagator import Propagator
from _pypropagate import finite_difference_1D,finite_difference_2D,crank_nicolson_2D
from ..coordinate_ndarray import CoordinateNDArray
import expresso.pycas as pc
import numpy as np

class FiniteDifferencesPropagator1D(Propagator):

    ndim = 1
    dtype = np.complex128
    
    def __init__(self,settings):
        super(FiniteDifferencesPropagator1D,self).__init__(settings)
    
        sb = settings.simulation_box
        pe = settings.paraxial_equation
        
        F,A,u_boundary = settings.get_unitless(pc.Tuple(pe.F,pe.A, pe.u_boundary).subs(sb.y,sb.fy))
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


class CrankNicolsonPropagator2D(Propagator):

    ndim = 2
    dtype = np.complex128

    def __init__(self,settings):
        super(CrankNicolsonPropagator2D,self).__init__(settings)

        sb = settings.simulation_box
        pe = settings.paraxial_equation

        F,A,u0,u_boundary = settings.get_unitless(( pe.F, pe.A, pe.u0, pe.u_boundary, ))

        dt = sb.dz/2
        rf = settings.get_unitless(pe.F*dt/2)

        lib = pc.ccompile(
                pc.FunctionDefinition("rf",(sb.x,sb.y,sb.z),rf,return_type=pc.Types.Complex),
                pc.FunctionDefinition("u_boundary",(sb.x,sb.y,sb.z),u_boundary,return_type=pc.Types.Complex,parallel=True),
                pc.FunctionDefinition("u0",(sb.x,sb.y,sb.z),u0,return_type=pc.Types.Complex),
        )

        self._rf = lib.rf
        self._ub = lib.u_boundary

        self._solver = crank_nicolson_2D()

        self._solver.resize(self._ny,self._nx)

        ra = settings.get_as(pe.A*dt/(sb.dy**2),complex)
        rb = settings.get_as(pe.A*dt/(sb.dx**2),complex)

        def init_solver():
            self._solver.ra.as_numpy().fill(ra)
            self._solver.rb.as_numpy().fill(rb)
            self._solver.rc.as_numpy().fill(0)
            self._solver.rd.as_numpy().fill(0)
            self._solver.re.as_numpy().fill(0)
            self._rf(*self._get_coordinates(),res = self._solver.rf.as_numpy())
            lib.u0(*self._get_coordinates(),res = self._solver.u.as_numpy())

        init_solver()
        self._solver.update()
        init_solver()

        self._set_initial_field(settings)

    def _update_boundary(self,i):
        z = self._nzmin + i * self._ndz
        self._ub(self._nxmin,self._get_y_coordinates(),z,res=self._solver.u.as_numpy()[0, :])
        self._ub(self._nxmax,self._get_y_coordinates(),z,res=self._solver.u.as_numpy()[-1,:])
        self._solver.u.as_numpy()[:, 0] = self._ub(self._get_x_coordinates(),self._nymin,z)
        self._solver.u.as_numpy()[:,-1] = self._ub(self._get_x_coordinates(),self._nymax,z)

    def _update(self):
        self._solver.update()
        self._rf(*self._get_coordinates(),res = self._solver.rf.as_numpy())

    def _step(self):
        self._update()
        self._update_boundary(self._i-0.5)
        self._solver.step_1()
        self._update()
        self._update_boundary(self._i)
        self._solver.step_2()

    def _get_field(self):
        return self._solver.u.as_numpy()

    def _set_field(self,field):
        self._solver.u.as_numpy()[:] = field



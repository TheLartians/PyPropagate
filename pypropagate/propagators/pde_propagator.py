
from .time_domain_propagator import TimeDomainPropagator
from _pypropagate import crank_nicolson_2D
import expresso.pycas as pc
import numpy as np

class PDEPropagator2D(TimeDomainPropagator):

    ndim = 2
    dtype = np.complex128

    def __init__(self,settings):
        super(PDEPropagator2D,self).__init__(settings)

        sb = settings.simulation_box
        pde = settings.PDE

        u_boundary,u0 = settings.get_unitless((pde.u_boundary,pde.u0))
        ra,rb,rc,rd,re,rf = settings.get_unitless(pc.Tuple(pde.ra,pde.rb,pde.rc,pde.rd,pde.re,pde.rf))

        lib = pc.ccompile(
                pc.FunctionDefinition("ra",(sb.x,sb.z,sb.t),ra,return_type=pc.Types.Complex),
                pc.FunctionDefinition("rb",(sb.x,sb.z,sb.t),rb,return_type=pc.Types.Complex),
                pc.FunctionDefinition("rc",(sb.x,sb.z,sb.t),rc,return_type=pc.Types.Complex),
                pc.FunctionDefinition("rd",(sb.x,sb.z,sb.t),rd,return_type=pc.Types.Complex),
                pc.FunctionDefinition("re",(sb.x,sb.z,sb.t),re,return_type=pc.Types.Complex),
                pc.FunctionDefinition("rf",(sb.x,sb.z,sb.t),rf,return_type=pc.Types.Complex),
                pc.FunctionDefinition("u_boundary",(sb.x,sb.z,sb.t),u_boundary,return_type=pc.Types.Complex,parallel=True),
                pc.FunctionDefinition("u0",(sb.x,sb.z,sb.t),u0,return_type=pc.Types.Complex,parallel=True),
        )

        self._lib = lib

        self._solver = crank_nicolson_2D()

        self._solver.resize(self._ny,self._nx)

        lib.u0(*self._get_coordinates(),res = self._solver.u.as_numpy())
        self._update()
        lib.u0(*self._get_coordinates(),res = self._solver.u.as_numpy())
        self._update()

        self._set_initial_field(settings)

    def _update(self):
        self._solver.update()
        self._lib.ra(*self._get_coordinates(),res = self._solver.ra.as_numpy())
        self._lib.rb(*self._get_coordinates(),res = self._solver.rb.as_numpy())
        self._lib.rc(*self._get_coordinates(),res = self._solver.rc.as_numpy())
        self._lib.rd(*self._get_coordinates(),res = self._solver.rd.as_numpy())
        self._lib.re(*self._get_coordinates(),res = self._solver.re.as_numpy())
        self._lib.rf(*self._get_coordinates(),res = self._solver.rf.as_numpy())

    def _step(self):
        self._update()
        self._solver.step_1()
        self._update()
        self._solver.step_2()

    def _get_field(self):
        return self._solver.u.as_numpy()

    def _set_field(self,field):
        self._solver.u.as_numpy()[:] = field



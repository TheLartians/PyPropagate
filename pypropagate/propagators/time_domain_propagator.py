
from ..solver import Solver

class TimeDomainPropagator(Solver):

    def __init__(self,settings):
        super(TimeDomainPropagator,self).__init__(settings)

        sb = settings.simulation_box

        self._x = sb.x
        self._y = sb.z
        self._z = sb.z
        self._t = sb.t
        self._nx = settings.get_as(sb.Nx,int)
        self._ny = settings.get_as(sb.Nz,int)
        self._nz = settings.get_as(sb.Nz,int)
        self._nt = settings.get_as(sb.Nt,int)
        self._xmin = settings.get_numeric(sb.xmin)
        self._xmax = settings.get_numeric(sb.xmax)
        self._ymin = settings.get_numeric(sb.zmin)
        self._ymax = settings.get_numeric(sb.zmax)
        self._zmin = settings.get_numeric(sb.zmin)
        self._zmax = settings.get_numeric(sb.zmax)
        self._tmin = settings.get_numeric(sb.tmin)
        self._tmax = settings.get_numeric(sb.tmax)

        self._nxmin = settings.get_as(sb.xmin,float)
        self._nxmax = settings.get_as(sb.xmax,float)
        self._nymin = settings.get_as(sb.zmin,float)
        self._nymax = settings.get_as(sb.zmax,float)
        self._nzmin = settings.get_as(sb.zmin,float)
        self._nzmax = settings.get_as(sb.zmax,float)
        self._ntmin = settings.get_as(sb.tmin,float)
        self._ntmax = settings.get_as(sb.tmax,float)

        self._ndx = settings.get_as(sb.dx,float)
        self._ndy = settings.get_as(sb.dz,float)
        self._ndz = settings.get_as(sb.dz,float)
        self._ndt = settings.get_as(sb.dt,float)


    def _set_initial_field(self,settings):
        sb = settings.simulation_box
        u0 = self._get_evaluators(settings.PDE.u0.subs(sb.t,sb.tmin),settings)
        initial = u0(*self._get_coordinates())
        self.__initial = initial
        self.set_field(initial)

    @property
    def _current_nt(self):
        return self._ntmin + self._i * self._ndt

    def _reset(self):
        self.set_field(self.__initial)

    def __get_x_coordinates(self):
        import numpy as np
        return np.linspace(self._nxmin,self._nxmax,self._nx)

    def _get_x_coordinates(self):
        try:
            return self.__x_coordinates
        except AttributeError:
            self.__x_coordinates = self.__get_x_coordinates()
            return self._get_x_coordinates()

    def __get_z_coordinates(self):
        import numpy as np
        return np.linspace(self._nzmin,self._nzmax,self._nz)

    def _get_z_coordinates(self):
        try:
            return self.__z_coordinates
        except AttributeError:
            self.__z_coordinates = self.__get_z_coordinates()
            return self._get_z_coordinates()

    def __get_xz_coordinates(self):
        import numpy as np
        npz,npx = np.meshgrid(np.linspace(self._nzmin,self._nzmax,self._nz),np.linspace(self._nxmin,self._nxmax,self._nx))
        return npx,npz

    def _get_coordinates(self):
        try:
            return self.__coordinates + [self._current_nt]
        except AttributeError:
            self.__coordinates = [self.__get_x_coordinates()] if self.ndim == 1 else list(self.__get_xz_coordinates())
            return self._get_coordinates()

    def _get_initial(self):
        return self.__initial

    def _get_evaluators(self,expressions,settings,**kwargs):
        import expresso.pycas as pc

        if not isinstance(expressions,(list,tuple)):
            return_single = True
            expressions = [expressions]
        else:
            return_single = False

        args = (self._x,self._z,self._t) if self.ndim == 2 else (self._x,self._t)
        sb = settings.simulation_box

        if self.ndim == 1:
            expressions = [settings.get_optimized(expr.subs(sb.y,sb.fy)) for expr in expressions]
        else:
            expressions = [settings.get_optimized(expr) for expr in expressions]

        definitions = [pc.FunctionDefinition('f%s' % i,args,expr,parallel=self.ndim>1,**kwargs)
                       for i,expr in enumerate(expressions)]

        if self.ndim == 1:
            lib = pc.ncompile(*definitions)
        else:
            lib = pc.ncompile(*definitions)

        res = [getattr(lib,'f%s' % i) for i in range(len(expressions))]

        if return_single:
            return res[0]
        return res


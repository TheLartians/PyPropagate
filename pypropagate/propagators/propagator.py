
from ..solver import Solver

class Propagator(Solver):
    
    def __init__(self,settings):
        super(Propagator,self).__init__(settings)

        sb = settings.simulation_box

        self._x = sb.x
        self._y = sb.y
        self._t = sb.z
        self._nx = settings.get_as(sb.Nx,int)
        if self.ndim > 1: self._ny = settings.get_as(sb.Ny,int)
        self._nt = settings.get_as(sb.Nz,int)
        self._xmin = settings.get_numeric(sb.xmin)
        self._xmax = settings.get_numeric(sb.xmax)
        if self.ndim > 1: self._ymin = settings.get_numeric(sb.ymin)
        if self.ndim > 1: self._ymax = settings.get_numeric(sb.ymax)
        self._tmin = settings.get_numeric(sb.zmin)
        self._tmax = settings.get_numeric(sb.zmax)

        self._nxmin = settings.get_as(sb.xmin,float)
        self._nxmax = settings.get_as(sb.xmax,float)
        if self.ndim > 1: self._nymin = settings.get_as(sb.ymin,float)
        if self.ndim > 1: self._nymax = settings.get_as(sb.ymax,float)
        self._nzmin = settings.get_as(sb.zmin,float)
        self._nzmax = settings.get_as(sb.zmax,float)

        self._ndx = settings.get_as(sb.dx,float)
        if self.ndim > 1: self._ndy = settings.get_as(sb.dy,float)
        self._ndz = settings.get_as(sb.dz,float)

        import expresso.pycas as pc
        pe = settings.paraxial_equation

        self._F_is_zero = settings.get_unitless( pe.F ) == pc.Zero
        self._F_is_constant = settings.get_unitless( pc.derivative(pe.F,sb.z) ) == pc.Zero

    def _set_initial_field(self,settings):
        sb = settings.simulation_box
        u0 = self._get_evaluators(settings.paraxial_equation.u0.subs(sb.z,sb.zmin),settings)
        initial = u0(*self._get_coordinates())
        self.__initial = initial
        self.set_field(initial)

    @property
    def _z(self):
        return self._t

    @property
    def _zmin(self):
        return self._tmin

    @property
    def _zmax(self):
        return self._tmax

    @property
    def _current_nz(self):
        return self._nzmin + self._i * self._ndz

    def _set_z(self,z):
        pass

    def _reset(self):
        self._set_z(self._nzmin)
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

    def __get_y_coordinates(self):
        import numpy as np
        return np.linspace(self._nymin,self._nymax,self._ny)

    def _get_y_coordinates(self):
        try:
            return self.__y_coordinates
        except AttributeError:
            self.__y_coordinates = self.__get_y_coordinates()
            return self._get_y_coordinates()

    def __get_xy_coordinates(self):
        import numpy as np
        npy,npx = np.meshgrid(np.linspace(self._nymin,self._nymax,self._ny),np.linspace(self._nxmin,self._nxmax,self._nx))
        return npx,npy

    def _get_coordinates(self):
        try:
            return self.__coordinates + [self._current_nz]
        except AttributeError:
            self.__coordinates = [self.__get_x_coordinates()] if self.ndim == 1 else list(self.__get_xy_coordinates())
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

        args = (self._x,self._y,self._z) if self.ndim == 2 else (self._x,self._z)
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

















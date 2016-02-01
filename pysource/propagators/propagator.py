
from ..solver import Solver

class Propagator(Solver):
    
    def __init__(self,settings):
        super(Propagator,self).__init__(settings)

        sb = settings.simulation_box

        self._x = sb.x
        self._y = sb.y
        self._t = sb.z
        self._nx = settings.get_as(sb.Nx,int)
        self._ny = settings.get_as(sb.Ny,int)
        self._nt = settings.get_as(sb.nz,int)
        self._xmin = settings.get_numeric(sb.xmin)
        self._xmax = settings.get_numeric(sb.xmax)
        self._ymin = settings.get_numeric(sb.ymin)
        self._ymax = settings.get_numeric(sb.ymax)
        self._tmin = settings.get_numeric(sb.zmin)
        self._tmax = settings.get_numeric(sb.zmax)
        self._downscale = settings.get_as(sb.downscale,int)

        self._nzmin = settings.get_as(sb.zmin,float)
        self._nzmax = settings.get_as(sb.zmax,float)

        self._ndx = settings.get_as(sb.dx,float)
        self._ndy = settings.get_as(sb.dy,float)
        self._ndz = settings.get_as(sb.dz,float)

    def _set_initial_field(self,initial,settings):

        if initial is None:
            from pycas import numpyfy
            pe = settings.paraxial_equation
            sb = settings.simulation_box

            if self.ndim == 1:
                u0 = settings.get_unitless(pe.u0.subs((sb.y,sb.fy),(sb.z,sb.zmin)))
                initial = numpyfy(u0)(x=self._get_x_coordinates(settings))

            if self.ndim == 2:
                u0 = settings.get_unitless(pe.u0.subs(sb.z,sb.zmin))
                initial = numpyfy(u0)(**{ n:c for n,c in zip(['x','y'],self._get_xy_coordinates(settings)) })

        self.__initial = initial
        self.set_field(initial)

    @property
    def _nz(self):
        return self._nzmin + self._i * self._ndz

    def _reset(self):
        self._set_z(self._nzmin)
        self.set_field(self.__initial)

    def _get_x_coordinates(self,settings):
        import numpy as np
        sb = settings.simulation_box
        xmin,xmax = settings.get_as((sb.xmin,sb.xmax),float)
        return np.linspace(xmin,xmax,self._nx)

    def _get_xy_coordinates(self,settings):
        import numpy as np
        sb = settings.simulation_box
        xmin,xmax,ymin,ymax = settings.get_as((sb.xmin,sb.xmax,sb.ymin,sb.ymax),float)
        npy,npx = np.meshgrid(np.linspace(ymin,ymax,self._ny),np.linspace(xmin,xmax,self._nx))
        return npx,npy

    def _get_initial(self):
        return self.__initial

    def get_z(self):
        return self.get_t()


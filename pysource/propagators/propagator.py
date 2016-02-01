
from ..solver import Solver

class Propagator(Solver):
    
    def __init__(self,settings):
        super(Propagator,self).__init__(settings)

        sb = settings.simulation_box

        self._x = sb.x
        self._y = sb.y
        self._t = sb.z
        self._nx = settings.get_as(sb.nx,int)
        self._ny = settings.get_as(sb.ny,int)
        self._nt = settings.get_as(sb.nz,int)
        self._xmin = settings.get_numeric(sb.xmin)
        self._xmax = settings.get_numeric(sb.xmax)
        self._ymin = settings.get_numeric(sb.ymin)
        self._ymax = settings.get_numeric(sb.ymax)
        self._tmin = settings.get_numeric(sb.zmin)
        self._tmax = settings.get_numeric(sb.zmax)
        self._downscale = settings.get_as(sb.downscale,int)

    def get_z(self):
        return self.get_t()


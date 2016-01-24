
from ..solver import Solver

class Propagator(Solver):
    
    def __init__(self,settings):
        super(Propagator,self).__init__(settings)

        self._x = settings.simulation_box.x
        self._y = settings.simulation_box.y
        self._t = settings.simulation_box.z
        self._dx = settings.get(settings.simulation_box.dx,float)
        self._dy = settings.get(settings.simulation_box.dy,float)
        self._dt = settings.get(settings.simulation_box.dz,float)
        self._nx = settings.get(settings.simulation_box.nx,int)
        self._ny = settings.get(settings.simulation_box.ny,int)
        self._nt = settings.get(settings.simulation_box.nz,int)
        self._xmin = settings.get_numeric(settings.simulation_box.xmin)
        self._xmax = settings.get_numeric(settings.simulation_box.xmax)
        self._ymin = settings.get_numeric(settings.simulation_box.ymin)
        self._ymax = settings.get_numeric(settings.simulation_box.ymax)
        self._tmin = settings.get_numeric(settings.simulation_box.zmin)
        self._tmax = settings.get_numeric(settings.simulation_box.zmax)
        self._downscale = settings.get(settings.simulation_box.downscale,int)
        
    def get_z(self):
        return self.get_t()


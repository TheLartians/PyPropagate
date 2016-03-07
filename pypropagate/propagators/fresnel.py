
from .propagator import Propagator
import numpy as np

class FresnelPropagator2D(Propagator):
    
    ndim = 2
    dtype = np.complex128
    
    def __init__(self,settings,thread_count = None):

        super(FresnelPropagator2D,self).__init__(settings)
        self._set_initial_field(settings)

        if thread_count == None:
            import multiprocessing
            thread_count = multiprocessing.cpu_count()
        self._thread_count = thread_count

        sb = settings.simulation_box
        pe = settings.partial_differential_equation

        import expresso.pycas as pc
        R, = self._get_evaluators([pc.exp(pe.F*sb.dz)],settings,return_type=pc.Types.Complex,parallel=not self._F_is_constant_in_z)

        if self._F_is_constant_in_z:
            self.__R_step = R(*self._get_coordinates())
        else:
            self.__R = R

        D = pc.numpyfy( settings.get_optimized( pc.exp(-pe.A*sb.dz*(pc.Symbol('kx')**2+pc.Symbol('ky')**2)) ) )

        import numpy as np

        fx = 2*np.pi/(self._nx*settings.get_as(sb.dx,float))
        fy = 2*np.pi/(self._ny*settings.get_as(sb.dy,float))

        kx,ky = np.meshgrid(
            fy*( self._ny/2.-np.abs(np.arange(self._ny)-self._ny/2.) ),
            fx*( self._nx/2.-np.abs(np.arange(self._nx)-self._nx/2.) )
        )

        self.__D_step = D(kx=kx,ky=ky)

    def _step(self):
        try:
            from pyfftw.interfaces.numpy_fft import fft2,ifft2
            kwargs = {'threads':self._thread_count}
        except ImportError:
            from numpy.fft import fft2,ifft2
            kwargs = {}

        if self._F_is_zero:
            freq = self.__initial_fft * self.__D_step**self._i
            self.__data = ifft2(freq,**kwargs)
        elif self._F_is_constant_in_z:
            freq = fft2(self.__data,**kwargs)
            freq *= self.__D_step
            self.__data = ifft2(freq,**kwargs)
            self.__data *= self.__R_step
        else:
            freq = fft2(self.__data,**kwargs)
            freq *= self.__D_step
            self.__data = ifft2(freq,**kwargs)
            self.__data *= self.__R(*self._get_coordinates())

    def _get_field(self):
        return self.__data
    
    def _set_field(self,field):
        import numpy as np
        from numpy.fft import fft2
        self.__data = field.astype(np.complex128)
        self.__initial_fft = fft2(field)

    
class FresnelPropagator1D(Propagator):

    ndim = 1
    dtype = np.complex128

    def __init__(self,settings):
        import expresso.pycas as pc

        super(FresnelPropagator1D,self).__init__(settings)
        self._set_initial_field(settings)

        sb = settings.simulation_box
        pe = settings.partial_differential_equation

        R, = self._get_evaluators([pc.exp(pe.F*sb.dz)],settings,return_type=pc.Types.Complex)

        D = pc.numpyfy( settings.get_optimized( pc.exp(-pe.A*sb.dz*(pc.Symbol('kx')**2)) ) )

        if self._F_is_constant_in_z:
            self.__R_step = R(*self._get_coordinates())
        else:
            self.__R = R

        import numpy as np
        fx = 2*np.pi/(self._nx*settings.get_as(sb.dx,float))
        kx = fx*( self._nx/2.-np.abs(np.arange(self._nx)- self._nx/2.) )
        self.__D_step = D(kx=kx)

    def _step(self):
        from numpy.fft import fft,ifft
        if self._F_is_zero:
            freq = self.__initial_fft * self.__D_step**self._i
            self.__data = ifft(freq)
        elif self._F_is_constant_in_z:
            freq = fft(self.__data)
            freq *= self.__D_step
            self.__data = ifft(freq)
            self.__data *= self.__R_step
        else:
            freq = fft(self.__data)
            freq *= self.__D_step
            self.__data = ifft(freq)
            self.__data *= self.__R(*self._get_coordinates())

    def _get_field(self):
        return self.__data

    def _set_field(self,field):
        from numpy.fft import fft
        import numpy as np
        self.__data = field.astype(np.complex128)
        self.__initial_fft = fft(field)

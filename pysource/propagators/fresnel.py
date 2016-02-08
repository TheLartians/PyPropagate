
from .propagator import Propagator
import numpy as np

class FresnelPropagator2D(Propagator):
    
    ndim = 2
    dtype = np.complex128
    
    def __init__(self,settings):

        super(FresnelPropagator2D,self).__init__(settings)
        self._set_initial_field(settings)

        sb = settings.simulation_box
        pe = settings.paraxial_equation

        import expresso.pycas as pc
        R,D = self._get_evaluators([pc.exp(pe.F*sb.dz),pc.exp(-pe.A*sb.dz*(sb.x**2+sb.y**2))],settings,return_type=pc.Types.Complex)

        self.__R = R
        if self._F_is_constant:
            self.__R_step = self.__R(*self._get_coordinates())

        import numpy as np

        fx = 2*np.pi/(self._nx*self._ndx)
        fy = 2*np.pi/(self._ny*self._ndx)

        kx,ky = np.meshgrid(
            fy*( self._ny/2.-np.abs(np.arange(self._ny,dtype = np.float64)-self._ny/2.) ),
            fx*( self._nx/2.-np.abs(np.arange(self._nx,dtype = np.float64)-self._nx/2.) )
        )

        self.__D_step = D(kx,ky,0)


    def _step(self):
        from numpy.fft import fft2,ifft2

        if self._F_is_zero:
            freq = fft2(self._get_initial())
            freq *= self.__D_step**self._i
            self.__data = ifft2(freq)
        else:
            freq = fft2(self.__data)
            freq *= self.__D_step
            self.__data = ifft2(freq)
            self.__data *= self.__R(*self._get_coordinates())

    def _get_field(self):
        return self.__data
    
    def _set_field(self,field):
        import numpy as np
        self.__data = field.astype(np.complex128)

    
class FresnelPropagator1D(Propagator):

    ndim = 1
    dtype = np.complex128

    def __init__(self,settings):

        super(FresnelPropagator1D,self).__init__(settings)
        self._set_initial_field(settings)

        sb = settings.simulation_box
        pe = settings.paraxial_equation

        import expresso.pycas as pc
        R,D = self._get_evaluators([pc.exp(pe.F*sb.dz),pc.exp(-pe.A*sb.dz*(sb.x**2))],settings,return_type=pc.Types.Complex)

        settings.get_unitless(pc.equal(pe.F,0))

        self.__R = R
        if self._F_is_constant:
            self.__R_step = self.__R(*self._get_coordinates())

        import numpy as np
        fx = 2*np.pi/(self._nx*self._ndx)
        kx = fx*( self._nx/2.-np.abs(np.arange(self._nx,dtype = np.float64)-self._nx/2.) )
        self.__D_step = D(kx,0)

    def _step(self):
        from numpy.fft import fft,ifft
        if self._F_is_zero:
            freq = fft(self._get_initial())
            freq *= self.__D_step**self._i
            self.__data = ifft(freq)
        else:
            freq = fft(self.__data)
            freq *= self.__D_step
            self.__data = ifft(freq)
            self.__data *= self.__R(*self._get_coordinates())

    def _get_field(self):
        return self.__data

    def _set_field(self,field):
        import numpy as np
        self.__data = field.astype(np.complex128)

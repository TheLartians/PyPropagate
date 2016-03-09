
from .propagator import Propagator
import numpy as np

class FresnelPropagator2D(Propagator):
    
    ndim = 2
    dtype = np.complex128
    
    def __init__(self,settings,thread_count = None):

        super(FresnelPropagator2D,self).__init__(settings)

        if thread_count == None:
            import multiprocessing
            thread_count = multiprocessing.cpu_count()
        self._thread_count = thread_count

        self._set_initial_field(settings)

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
    
    def get_fft(self):
        try:
	    import pyfftw
	    if not hasattr(self,'__fftw_ffts'):
	    	a = pyfftw.empty_aligned((self._nx, self._ny), dtype='complex128')
	    	b = pyfftw.empty_aligned((self._nx, self._ny), dtype='complex128')
	    	fft2 = pyfftw.FFTW(a, b, axes=(0,1), threads=self._thread_count, direction = 'FFTW_FORWARD')
	    	ifft2 = pyfftw.FFTW(b, a, axes=(0,1), threads=self._thread_count, direction = 'FFTW_BACKWARD')
		self.__fftw_ffts = fft2,ifft2
	    return self.__fftw_ffts
        except ImportError:
            from numpy.fft import fft2,ifft2
            return fft2,ifft2

    def _step(self):
        fft2,ifft2 = self.get_fft()        

        if self._F_is_constant:
            self.__freq_data *= self.__D_step * self.__R_step
            self.__data = ifft2(self.__freq_data)
        elif self._F_is_constant_in_z:
            self.__freq_data = fft2(self.__data)
            self.__freq_data *= self.__D_step
            self.__data = ifft2(self.__freq_data)
            self.__data *= self.__R_step
        else:
            self.__freq_data = fft2(self.__data)
            self.__freq_data *= self.__D_step
            self.__data = ifft2(self.__freq_data)
            self.__data *= self.__R(*self._get_coordinates())

    def _get_field(self):
        return self.__data
    
    def _set_field(self,field):
        import numpy as np
        fft2,ifft2 = self.get_fft()        
        self.__data = field.astype(np.complex128)
	self.__freq_data = fft2(self.__data)


class FresnelPropagator1D(Propagator):

    ndim = 1
    dtype = np.complex128

    def __init__(self,settings,thread_count=1):
        import expresso.pycas as pc

        super(FresnelPropagator1D,self).__init__(settings)
	self._thread_count = thread_count

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

    def get_fft(self):
        try:
	    import pyfftw
	    if not hasattr(self,'__fftw_ffts'):
	    	a = pyfftw.empty_aligned(self._nx, dtype='complex128')
	    	b = pyfftw.empty_aligned(self._nx, dtype='complex128')
	    	fft = pyfftw.FFTW(a, b, axes=(0,),  threads=self._thread_count, direction = 'FFTW_FORWARD')
	    	ifft = pyfftw.FFTW(b, a, axes=(0,), threads=self._thread_count, direction = 'FFTW_BACKWARD')
		self.__fftw_ffts = fft,ifft
	    return self.__fftw_ffts
        except ImportError:
            from numpy.fft import fft,ifft
            return fft,ifft

    def _step(self):
        fft,ifft = self.get_fft()        

        if self._F_is_constant:
            self.__freq_data *= self.__D_step * self.__R_step
            self.__data = ifft(self.__freq_data)
        elif self._F_is_constant_in_z:
            self.__freq_data = fft(self.__data)
            self.__freq_data *= self.__D_step
            self.__data = ifft(self.__freq_data)
            self.__data *= self.__R_step
        else:
            self.__freq_data = fft(self.__data)
            self.__freq_data *= self.__D_step
            self.__data = ifft(self.__freq_data)
            self.__data *= self.__R(*self._get_coordinates())

    def _get_field(self):
        return self.__data

    def _set_field(self,field):
        fft,ifft = self.get_fft()
        import numpy as np
        self.__data = field.astype(np.complex128)
        self.__freq_data = fft(self.__data)


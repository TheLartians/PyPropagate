
from .propagator import Propagator
#from ..create_data import create_front_data,create_1D_front_data
#from .._pypropagate import create_2D_field_from_function,create_1D_field_from_function
import numpy as np

def fresnel_propagate( bitmap, z, k, (dx, dy) ):
    from numpy.fft import fft2,ifft2
    import numpy as np

    shape = bitmap.shape

    fx = 2*np.pi/(bitmap.shape[0]*dx) 
    fy = 2*np.pi/(bitmap.shape[1]*dy)   

    qxmax,qymax = (shape[0]/2*fx, shape[1]/2*fy)

    if not ( k*k-qxmax*qxmax-qymax*qymax >= 0 ):
        raise ValueError('small angle approximation not valid')

    qx,qy = np.meshgrid( 
        fy*( shape[1]/2.-np.abs(np.arange(shape[1],dtype = np.float64)-shape[1]/2.) ),
        fx*( shape[0]/2.-np.abs(np.arange(shape[0],dtype = np.float64)-shape[0]/2.) )
    )
    
    kappa = ( k*k - qx*qx - qy*qy )**0.5

    I = fft2(bitmap)
    I *= np.exp(-1.j * z * kappa )
    I = ifft2(I)
    
    return I

class FresnelPropagator2D(Propagator):
    
    ndim = 2
    dtype = np.complex128
    
    def __init__(self,settings,initial = None):
        from sympy import exp
        
        super(FresnelPropagator2D,self).__init__(settings)
        
        psi0 = settings.wave_equation.psi0.subs(settings.simulation_box.z,settings.simulation_box.zmin)
        
        if initial != None:
            self.set_field(initial)
        else:
            self.__data = create_front_data(psi0,settings,z=settings.simulation_box.zmin).data

        self.__initial = self.__data
        
        self._transform = settings.get_transform()
        
        s = settings.symbols
        update_expr = settings.get_unitless(exp(-1j * s.k * s.dz * s.n))
        
        self._field_update = compile_sympy_expression(update_expr,(s.x,s.y),globals = s.z)
        self.__z_parameter = self._field_update.globals[s.z.name]
        self.__z_parameter.value = settings.get(s.zmin,float)
        self.__ndz = settings.get(s.dz,float)
        self.__nzmin = settings.get(s.zmin,float)

        self._k = settings.get(settings.wave_equation.k,float)
        
        xmin = settings.get(s.xmin,float)
        xmax = settings.get(s.xmax,float)
        ymin = settings.get(s.ymin,float)
        ymax = settings.get(s.ymax,float)
        
        mul = lambda: create_2D_field_from_function(self._field_update.address,self._nx,self._ny,xmin,xmax,ymin,ymax)
        self.__update_field = mul
        
    def reset(self):
        self.__z_parameter.value = self.__nzmin
        self.__data = self.__initial
        
    def _step(self):
        self.__z_parameter.value = self.__z_parameter.value + self.__ndz
        update = self.__update_field()
        self.__data *= update
        self.__data = fresnel_propagate(self.__data, self._dt ,  self._k ,  (self._dx,self._dy ))
    
    def _get_field(self):
        return self.__data
    
    def _set_field(self,field):
        self.__data = field

    
def fresnel_propagate_1D( bitmap, z, k, dx ):
    from numpy.fft import fft,ifft
    import numpy as np

    shape = bitmap.shape

    fx = 2*np.pi/(bitmap.shape[0]*dx) 

    qxmax = shape[0]/2*fx

    if not ( k*k-qxmax*qxmax >= 0 ):
        raise ValueError('small angle approximation not valid')

    qx = fx*( shape[0]/2.-np.abs(np.arange(shape[0],dtype = np.float64)-shape[0]/2.) )
    
    kappa = ( k*k - qx*qx )**0.5

    I = fft(bitmap)
    I *= np.exp(-1.j * z * kappa )
    I = ifft(I)
    
    return I

class FresnelPropagator1D(Propagator):
    
    ndim = 1
    dtype = np.complex128
    
    def __init__(self,settings,initial = None):
        from sympy import exp
        
        super(FresnelPropagator1D,self).__init__(settings)
        
        psi0 = settings.wave_equation.psi0.subs(settings.simulation_box.z,settings.simulation_box.zmin)
        
        if initial != None:
            self.set_field(initial)
        else:
            self.__data = create_1D_front_data(psi0,settings,z=settings.simulation_box.zmin).data
       
        self.__initial = self.__data

        self._transform = settings.get_transform()
        
        s = settings.symbols
        update_expr = settings.get_unitless(exp(-1j * s.k * s.dz * s.n).subs(s.y,(s.ymin+s.ymax)/2))
        
        self._field_update = compile_sympy_expression(update_expr,(s.x,),globals = s.z)
        self.__z_parameter = self._field_update.globals[s.z.name]
        self.__z_parameter.value = settings.get(s.zmin,float)
        self.__ndz = settings.get(s.dz,float)
        self.__nzmin = settings.get(s.zmin,float)

        self._k = settings.get(settings.wave_equation.k,float)
        
        xmin = settings.get(s.xmin,float)
        xmax = settings.get(s.xmax,float)
        
        mul = lambda: create_1D_field_from_function(self._field_update.address,self._nx,xmin,xmax)
        self.__update_field = mul
        
    def reset(self):
        self.__z_parameter.value = self.__nzmin
        self.__data = self.__initial
        
    def _step(self):
        self.__z_parameter.value = self.__z_parameter.value + self.__ndz
        update = self.__update_field()
        self.__data *= update
        self.__data = fresnel_propagate_1D(self.__data, self._dt ,  self._k , self._dx)
    
    def _get_field(self):
        return self.__data
    
    def _set_field(self,field):
        self.__data = field

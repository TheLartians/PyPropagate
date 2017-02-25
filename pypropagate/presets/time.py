
def fourier_transform(array,axis,new_axis,inverse=False):
    import numpy as np
    from numpy.fft import fftshift,ifftshift
    from expresso.pycas import pi
    from ..coordinate_ndarray import CoordinateNDArray
 
    try:
        from pyfftw.interfaces.numpy_fft import fft,ifft
    except ImportError:
        from numpy.fft import fft,ifft

    axi = array.axis.index(axis)

    if not inverse:
        new_data = fftshift(fft(array.data,axis=axi),axes=[axi])
        new_data *= 1/np.sqrt(2*np.pi)  
    else:
        new_data = ifft(ifftshift(array.data,axes=[axi]),axis=axi)
        new_data *= np.sqrt(2*np.pi) 

    sw = array.bounds[axi][1] - array.bounds[axi][0]
    tmin,tmax = array.evaluate((-(pi*array.shape[axi])/sw,
                                 (pi*array.shape[axi])/sw))

    new_bounds = [(b[0],b[1]) if i!=axi else (tmin,tmax) for i,b in enumerate(array.bounds)]
    new_axis = [a  if i!=axi else new_axis for i,a in enumerate(array.axis)]

    return CoordinateNDArray(new_data,new_bounds,new_axis,array.evaluate)

def inverse_fourier_transform(*args):
    return fourier_transform(*args,inverse=True)

def periodic_envelope_propagation(field,omega0,s=1,z=None,omega=None,t=None):

    import numpy as np
    import expresso.pycas as pc
    from .. import units
    
    if z is None:
        z = pc.Symbol('z')
    if t is None:
        t = pc.Symbol('t')
    if omega is None:
        omega = pc.Symbol('omega')

    zindex = field.get_axis_index(z)
    omegaindex = field.get_axis_index(omega)

    omegamin,omegamax = field.bounds[omegaindex]
    zmin, zmax = field.bounds[zindex]
    sz = zmax - zmin

    ukmin = float(field.evaluate( (omegamin-omega0)/units.c*sz ))
    ukmax = float(field.evaluate( (omegamax-omega0)/units.c*sz ))
    uzmin = 0
    uzmax = 1

    nz,ik = np.meshgrid(np.linspace(uzmin,uzmax,field.shape[zindex]),np.linspace(1j*ukmin,1j*ukmax,field.shape[omegaindex]))
    factor = np.exp(-ik*nz/float(s))
    
    transform = field * factor
    del factor,nz,ik

    tdfield = fourier_transform(transform, omega, t, inverse=True)
    tdfield.bounds[omegaindex] = [(b - tdfield.bounds[omegaindex][0]).evaluate() for b in tdfield.bounds[omegaindex]]
    
    return tdfield




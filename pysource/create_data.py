import _pypropagate
from .compile_sympy import compile_sympy_expression,ctypes
from .coordinate_ndarray import CoordinateNDArray

def create_side_data(expr,settings, y = None):
    settings.initialize()
    s = settings.simulation_box
    [xmin,xmax],[zmin,zmax] = settings.get(([s.xmin,s.xmax],[s.zmin,s.zmax]),float)
    size = settings.get((s.nx/s.downscale,s.nz/s.downscale+1),int)
    
    if y is None:
        y = (s.ymin+s.ymax)/2
    
    substituted = settings.get_unitless(expr.subs(s.y,y))
    compiled  = compile_sympy_expression(substituted,(s.x,s.z),return_type=ctypes.c_complex)
    
    data = _finitedifferences.create_2D_field_from_function(compiled.address,size[0],size[1],xmin,xmax,zmin,zmax)
        
    res = CoordinateNDArray(data,settings.get_numeric(([s.xmin,s.xmax],[s.zmin,s.zmax])),(s.x,s.z))
    return res

def create_front_data(expr,settings,z = None):
    settings.initialize()
    s = settings.simulation_box
    [xmin,xmax],[ymin,ymax] = settings.get(([s.xmin,s.xmax],[s.ymin,s.ymax]),float)
    size = settings.get((s.nx/s.downscale,s.ny/s.downscale),int)
    
    if z is None:
        z = (s.zmin)
    
    substituted = settings.get_unitless(expr.subs(s.z,z))
    
    compiled  = compile_sympy_expression(substituted,(s.x,s.y),return_type=ctypes.c_complex)
    
    data = _finitedifferences.create_2D_field_from_function(compiled.address,size[0],size[1],xmin,xmax,ymin,ymax)
        
    res = CoordinateNDArray(data,settings.get_numeric(([s.xmin,s.xmax],[s.ymin,s.ymax])),(s.x,s.y))
    return res
    
def create_1D_front_data(expr,settings,y = None, z = None):
    settings.initialize()
    s = settings.simulation_box
    [xmin,xmax],[ymin,ymax] = settings.get(([s.xmin,s.xmax],[s.ymin,s.ymax]),float)
    size = settings.get(s.nx/s.downscale,int)
    
    if z is None:
        z = (s.zmin)
    if y is None:
        y = (s.ymax + s.ymin)/2

    substituted = settings.get_unitless(expr.subs([(s.z,z),(s.y,y)]))
    
    compiled  = compile_sympy_expression(substituted,(s.x,s.y),return_type=ctypes.c_complex)
    
    data = _finitedifferences.create_1D_field_from_function(compiled.address,size,xmin,xmax)
        
    res = CoordinateNDArray(data,settings.get_numeric([[s.xmin,s.xmax]]),(s.x,))
    return res

def create_volume_data(expr,settings,size = None):
    settings.initialize()
    import numpy as np

    s = settings.symbols
    
    [xmin,xmax],[ymin,ymax],[zmin,zmax] = settings.get(([s.xmin,s.xmax],[s.ymin,s.ymax],[s.zmin,s.zmax]),float)
    dz = settings.get(s.dz,float)
    
    if size == None:
        size = settings.get((s.nx/s.downscale,s.ny/s.downscale,s.nz/s.downscale+1),int)

    data = np.zeros(size,dtype = np.complex128)

    substituted = settings.get_unitless(expr)
    compiled  = compile_sympy_expression(substituted,(s.x,s.y),return_type=ctypes.c_complex,globals=[s.z])
    
    zg = compiled.globals[s.z.name]
    
    for i in range(size[2]):
        zi = zmin + i*dz
        zg.value = zi
        data[:,:,i] = _finitedifferences.create_2D_field_from_function(compiled.address,size[0],size[1],xmin,xmax,ymin,ymax)

    res = CoordinateNDArray(data,settings.get_numeric(([s.xmin,s.xmax],[s.ymin,s.ymax],[s.zmin,s.zmax])),(s.x,s.y,s.z))
    return res



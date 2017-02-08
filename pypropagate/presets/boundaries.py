

def set_1D_boundary_condition(settings):
    from expresso.pycas import exp

    s = settings.simulation_box
    pe = settings.partial_differential_equation

    pe.u_boundary = pe.u0.subs(s.z,s.zmin) * exp(pe.F.subs(s.z,s.zmin)*s.z)


def set_plane_wave_initial_conditions(settings):
    """Sets the boundary conditions to a plane wave with intensity 1.
    The boundary are set to the index of refraction at z=0."""

    s = settings.simulation_box
    pe = settings.partial_differential_equation

    pe.u0 = 1
    set_1D_boundary_condition(settings)


def add_padding(array,factor,mode = 'edge',**kwargs):
    import numpy as np
    from ..coordinate_ndarray import CoordinateNDArray

    padding_points = [[int(x*factor)]*2 for x in array.data.shape]
    new_data = np.pad(array.data,padding_points,mode,**kwargs)

    extension = [d*p[0] for d,p in zip(array._dbounds,padding_points)]
    new_bounds = [(b-i,e+i) for i,(b,e) in zip(extension,array.bounds)]

    return CoordinateNDArray(new_data,new_bounds,array.axis,array.evaluate)

def set_initial(settings,initial_array):
    import expresso.pycas as pc
    from ..coordinate_ndarray import CoordinateNDArray

    if isinstance(initial_array,CoordinateNDArray):
        initial = pc.array("initial",initial_array.data)
    else:
        initial = pc.array("initial",initial_array)

    sb = settings.simulation_box

    if tuple(initial_array.axis) == (sb.x,):
        settings.partial_differential_equation.u0 = initial(sb.xi)
    elif tuple(initial_array.axis) == (sb.x,sb.y):
        settings.partial_differential_equation.u0 = initial(sb.yi,sb.xi)
        sb.Ny = initial_array.shape[1]
        if isinstance(initial_array,CoordinateNDArray):
            sb.unlock('ymin')
            sb.unlock('ymax')
            sb.unlock('sy')
            sb.ymin = initial_array.bounds[1][0]
            sb.ymax = initial_array.bounds[1][1]
            sb.sy = sb.ymax - sb.ymin
            sb.lock('ymin','defined by initial array')
            sb.lock('ymax','defined by initial array')
            sb.lock('sy','defined by ymin and ymax')

    else:
        raise ValueError('initial array axis must be (x,) or (x,y)')

    sb.Nx = initial_array.shape[0]

    if isinstance(initial_array,CoordinateNDArray):
        sb.unlock('xmin')
        sb.unlock('xmax')
        sb.unlock('sx')
        sb.xmin = initial_array.bounds[0][0]
        sb.xmax = initial_array.bounds[0][1]
        sb.sx = sb.xmax - sb.xmin
        sb.lock('xmin','defined by initial array')
        sb.lock('xmax','defined by initial array')
        sb.lock('sx','defined by xmin and xmax')


